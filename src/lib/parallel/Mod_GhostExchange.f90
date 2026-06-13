!>@brief Ghost cell MPI communication for inter-block boundaries.
!> Builds a communication schedule from the BC connectivity array,
!> then provides non-blocking exchange routines for P and R arrays.
!> Messages are aggregated per remote rank to minimize latency.
!> In serial mode (USE_MPI not defined), all routines are no-ops.
module ARES_Mod_GhostExchange
  use iso_fortran_env, only: R8 => real64, I4 => int32
  use ARES_Mod_MPI

  implicit none
  private

  !> Single ghost cell exchange entry
  type :: ghost_entry_type
    integer :: bc_idx          !< Index into domain%bc(:)
    integer :: local_block     !< Block on this rank
    integer :: remote_block    !< Block on the remote rank
    integer :: remote_rank     !< MPI rank of the remote block
    integer :: src_i, src_j, src_k, src_f  !< Source cell coords and face
    integer :: dst_i, dst_j, dst_k, dst_f  !< Destination cell coords and face
  end type ghost_entry_type

  !> Per-rank aggregation info
  type :: rank_group_type
    integer :: rank = -1       !< Remote MPI rank
    integer :: offset = 0      !< Start index in entry list (1-based)
    integer :: count  = 0      !< Number of entries for this rank
  end type rank_group_type

  !> Face group: contiguous run of entries sharing (block, face).
  !> Precomputes guide offsets; for rectangular faces enables direct coordinate sweep.
  type :: face_group_type
    integer :: block            !< Local block index
    integer :: face             !< Face number (1-6)
    integer :: entry_start      !< First entry index in send/recv list
    integer :: entry_count      !< Number of entries
    integer :: di, dj, dk       !< Precomputed guide(face,:)
  end type face_group_type

  !> Communication schedule: entries grouped by remote rank for aggregated messaging
  type :: ghost_schedule_type
    integer :: n_send = 0, n_recv = 0
    type(ghost_entry_type), allocatable :: send_list(:)
    type(ghost_entry_type), allocatable :: recv_list(:)
    logical :: built = .false.

    ! Per-rank aggregation for sends
    integer :: n_send_ranks = 0
    type(rank_group_type), allocatable :: send_groups(:)

    ! Per-rank aggregation for receives
    integer :: n_recv_ranks = 0
    type(rank_group_type), allocatable :: recv_groups(:)

    ! Face groups for aggregated pack/unpack (reduce per-cell overhead)
    integer :: n_send_faces = 0
    type(face_group_type), allocatable :: send_faces(:)
    integer :: n_recv_faces = 0
    type(face_group_type), allocatable :: recv_faces(:)

    ! Pre-allocated MPI buffers (one contiguous buffer per rank group)
    real(R8), allocatable :: P_send_buf(:)   !< (nprim*gc * n_send)
    real(R8), allocatable :: P_recv_buf(:)   !< (nprim*gc * n_recv)
    real(R8), allocatable :: Pg_send_buf(:)  !< (nprim*4  * n_send)
    real(R8), allocatable :: Pg_recv_buf(:)  !< (nprim*4  * n_recv)
    integer, allocatable  :: send_req(:)     !< (n_send_ranks)
    integer, allocatable  :: recv_req(:)     !< (n_recv_ranks)
    integer, allocatable  :: mpi_stat(:,:)   !< (MPI_STATUS_SIZE, max(n_send_ranks, n_recv_ranks))
    ! Persistent MPI requests for P exchange (initialized once, started each step)
    integer, allocatable  :: P_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable  :: P_recv_req_pers(:)  !< (n_recv_ranks)
    integer, allocatable  :: R_send_req_pers(:)  !< (n_send_ranks)
    integer, allocatable  :: R_recv_req_pers(:)  !< (n_recv_ranks)
    logical :: persistent_P_init = .false.
    logical :: persistent_R_init = .false.
  end type ghost_schedule_type

  type(ghost_schedule_type), public :: ghost_sched

  public :: build_ghost_schedule
  public :: build_local_bc_index
  public :: cleanup_ghost_schedule
  public :: exchange_ghost_P_post_recv, exchange_ghost_P_pack
  public :: exchange_ghost_P_post_send, exchange_ghost_P_wait_unpack
  public :: exchange_ghost_P_wait_send
  public :: exchange_ghost_R
  public :: exchange_ghost_Pg
  public :: gather_P_to_root
  public :: gather_diagnostic_to_root
  public :: scatter_P_from_root
  public :: mpi_io_barrier
  public :: Ghost_Interrank

contains


  !> Build the communication schedule by scanning all BC entries of type 1 (connection).
  !> Entries are sorted by remote rank for aggregated MPI messaging.
  !> Must be called after partition_blocks and domain setup.
  subroutine build_ghost_schedule(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, nprim

    implicit none
    type(ARES_domain_type), intent(in) :: domain
    ! Local
    integer :: i, ns, nr, bm, bs

    if (mpi_size_ <= 1) then
      ghost_sched%built = .true.
      return
    end if

#ifdef USE_MPI
    ! Count send and recv entries
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 101) cycle
      bm = domain%bc(i)%b   ! destination block
      bs = domain%bc(i)%bs  ! source block
      if (.not. allocated(block_owner)) cycle
      ! If I own the destination but not the source => I need to RECV
      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
      end if
      ! If I own the source but not the destination => I need to SEND
      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
      end if
    end do

    ghost_sched%n_send = ns
    ghost_sched%n_recv = nr
    if (allocated(ghost_sched%send_list)) deallocate(ghost_sched%send_list)
    if (allocated(ghost_sched%recv_list)) deallocate(ghost_sched%recv_list)
    allocate(ghost_sched%send_list(ns))
    allocate(ghost_sched%recv_list(nr))

    ! Fill send and recv lists
    ns = 0; nr = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 101) cycle
      bm = domain%bc(i)%b
      bs = domain%bc(i)%bs

      if (is_local_block(bm) .and. (.not. is_local_block(bs))) then
        nr = nr + 1
        ghost_sched%recv_list(nr)%bc_idx = i
        ghost_sched%recv_list(nr)%local_block  = bm
        ghost_sched%recv_list(nr)%remote_block = bs
        ghost_sched%recv_list(nr)%remote_rank  = block_owner(bs)
        ghost_sched%recv_list(nr)%src_i = domain%bc(i)%is
        ghost_sched%recv_list(nr)%src_j = domain%bc(i)%js
        ghost_sched%recv_list(nr)%src_k = domain%bc(i)%ks
        ghost_sched%recv_list(nr)%src_f = domain%bc(i)%fs
        ghost_sched%recv_list(nr)%dst_i = domain%bc(i)%i
        ghost_sched%recv_list(nr)%dst_j = domain%bc(i)%j
        ghost_sched%recv_list(nr)%dst_k = domain%bc(i)%k
        ghost_sched%recv_list(nr)%dst_f = domain%bc(i)%f
      end if

      if (is_local_block(bs) .and. (.not. is_local_block(bm))) then
        ns = ns + 1
        ghost_sched%send_list(ns)%bc_idx = i
        ghost_sched%send_list(ns)%local_block  = bs
        ghost_sched%send_list(ns)%remote_block = bm
        ghost_sched%send_list(ns)%remote_rank  = block_owner(bm)
        ghost_sched%send_list(ns)%src_i = domain%bc(i)%is
        ghost_sched%send_list(ns)%src_j = domain%bc(i)%js
        ghost_sched%send_list(ns)%src_k = domain%bc(i)%ks
        ghost_sched%send_list(ns)%src_f = domain%bc(i)%fs
        ghost_sched%send_list(ns)%dst_i = domain%bc(i)%i
        ghost_sched%send_list(ns)%dst_j = domain%bc(i)%j
        ghost_sched%send_list(ns)%dst_k = domain%bc(i)%k
        ghost_sched%send_list(ns)%dst_f = domain%bc(i)%f
      end if
    end do

    ! Sort entries by remote rank and build rank groups
    call sort_entries_by_rank(ghost_sched%send_list, ns)
    call sort_entries_by_rank(ghost_sched%recv_list, nr)
    call build_rank_groups(ghost_sched%send_list, ns, ghost_sched%send_groups, ghost_sched%n_send_ranks)
    call build_rank_groups(ghost_sched%recv_list, nr, ghost_sched%recv_groups, ghost_sched%n_recv_ranks)

    ! Build face groups for aggregated pack/unpack
    call build_face_groups_send(ghost_sched%send_list, ns, &
                                ghost_sched%send_faces, ghost_sched%n_send_faces)
    call build_face_groups_recv(ghost_sched%recv_list, nr, &
                                ghost_sched%recv_faces, ghost_sched%n_recv_faces)

    ! Pre-allocate MPI buffers (reused every RK step)
    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns > 0 .or. nr > 0) then
      block
        use mpi, only: MPI_STATUS_SIZE
        integer :: max_ranks
        max_ranks = max(ghost_sched%n_send_ranks, ghost_sched%n_recv_ranks, 1)
        allocate(ghost_sched%P_send_buf(nprim*gc * max(ns,1)))
        allocate(ghost_sched%P_recv_buf(nprim*gc * max(nr,1)))
        allocate(ghost_sched%Pg_send_buf(nprim*4 * max(ns,1)))
        allocate(ghost_sched%Pg_recv_buf(nprim*4 * max(nr,1)))
        allocate(ghost_sched%send_req(max(ghost_sched%n_send_ranks,1)))
        allocate(ghost_sched%recv_req(max(ghost_sched%n_recv_ranks,1)))
        allocate(ghost_sched%mpi_stat(MPI_STATUS_SIZE, max_ranks))
      end block
    end if

    ! Initialize persistent MPI requests for P and R exchanges
    call init_persistent_P_requests()
    call init_persistent_R_requests()

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A)') ' Ghost schedule: ', ghost_sched%n_send_ranks, &
        ' send ranks, ', ghost_sched%n_recv_ranks, ' recv ranks'
    end if
#endif

    ghost_sched%built = .true.
  end subroutine build_ghost_schedule


  !> Build pre-filtered list of BC indices where Bm is local to this rank.
  !> Avoids O(nbound) loop scanning in hot paths (Fill_Ghost_Cell, etc.).
  subroutine build_local_bc_index(domain)
    use ARES_Advanced_Types_m

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, n, ns

    ! Count local BC entries (Bm is local) and local source entries (Bs is local, type-101)
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) n = n + 1
      if (domain%bc(i)%type == 101) then
        if (is_local_block(domain%bc(i)%bs)) ns = ns + 1
      end if
    end do

    ! Build local_bc_idx (Bm is local)
    domain%n_local_bc = n
    if (allocated(domain%local_bc_idx)) deallocate(domain%local_bc_idx)
    allocate(domain%local_bc_idx(n))

    ! Build local_bs_idx (Bs is local, type-1 only — for Pg computation)
    domain%n_local_bs = ns
    if (allocated(domain%local_bs_idx)) deallocate(domain%local_bs_idx)
    allocate(domain%local_bs_idx(ns))

    ! Fill both index arrays in a single pass
    n = 0; ns = 0
    do i = 1, domain%nbound
      if (is_local_block(domain%bc(i)%b)) then
        n = n + 1
        domain%local_bc_idx(n) = i
      end if
      if (domain%bc(i)%type == 101) then
        if (is_local_block(domain%bc(i)%bs)) then
          ns = ns + 1
          domain%local_bs_idx(ns) = i
        end if
      end if
    end do
  end subroutine build_local_bc_index


  !> Free persistent MPI requests and deallocate schedule arrays.
  subroutine cleanup_ghost_schedule()
    implicit none
#ifdef USE_MPI
    call cleanup_persistent_requests()
#endif
  end subroutine cleanup_ghost_schedule


#ifdef USE_MPI
  !> Initialize persistent MPI requests for P field exchange.
  !> Called once from build_ghost_schedule. Uses MPI_SEND_INIT/RECV_INIT.
  subroutine init_persistent_P_requests()
    use ARES_Global_m, only: gc, nprim
    use mpi

    implicit none
    integer :: r, ierr, tag, buf_pos, entry_size

    entry_size = nprim * gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%P_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%P_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Initialize persistent receives
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank
      call MPI_RECV_INIT(ghost_sched%P_recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%P_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize persistent sends
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_
      call MPI_SEND_INIT(ghost_sched%P_send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%P_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_P_init = .true.
  end subroutine init_persistent_P_requests


  !> Initialize persistent MPI requests for R field exchange.
  !> Reuses P_send_buf/P_recv_buf (same size). Uses distinct tags (offset 30000).
  subroutine init_persistent_R_requests()
    use ARES_Global_m, only: gc, nprim
    use mpi

    implicit none
    integer :: r, ierr, tag, buf_pos, entry_size
    integer, parameter :: TAG_OFFSET = 30000

    entry_size = nprim * gc

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return

    allocate(ghost_sched%R_recv_req_pers(max(ghost_sched%n_recv_ranks, 1)))
    allocate(ghost_sched%R_send_req_pers(max(ghost_sched%n_send_ranks, 1)))

    ! Initialize persistent receives for R
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + TAG_OFFSET
      call MPI_RECV_INIT(ghost_sched%P_recv_buf(buf_pos), &
                         ghost_sched%recv_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%recv_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_recv_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Initialize persistent sends for R
    do r = 1, ghost_sched%n_send_ranks
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + TAG_OFFSET
      call MPI_SEND_INIT(ghost_sched%P_send_buf(buf_pos), &
                         ghost_sched%send_groups(r)%count * entry_size, &
                         MPI_DOUBLE_PRECISION, &
                         ghost_sched%send_groups(r)%rank, tag, &
                         MPI_COMM_WORLD, ghost_sched%R_send_req_pers(r), ierr)
      call check_mpi_error(ierr)
    end do

    ghost_sched%persistent_R_init = .true.
  end subroutine init_persistent_R_requests


  !> Free all persistent MPI requests.
  subroutine cleanup_persistent_requests()
    use mpi
    implicit none
    integer :: ierr, r

    if (ghost_sched%persistent_P_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%P_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%P_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_P_init = .false.
    end if

    if (ghost_sched%persistent_R_init) then
      do r = 1, ghost_sched%n_send_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_send_req_pers(r), ierr)
      end do
      do r = 1, ghost_sched%n_recv_ranks
        call MPI_REQUEST_FREE(ghost_sched%R_recv_req_pers(r), ierr)
      end do
      ghost_sched%persistent_R_init = .false.
    end if
  end subroutine cleanup_persistent_requests


  !> Sort ghost entries by remote_rank using insertion sort (lists are small).
  subroutine sort_entries_by_rank(list, n)
    implicit none
    type(ghost_entry_type), intent(inout) :: list(:)
    integer, intent(in) :: n
    ! Local
    integer :: i, j
    type(ghost_entry_type) :: tmp

    do i = 2, n
      tmp = list(i)
      j = i - 1
      do while (j >= 1 .and. list(j)%remote_rank > tmp%remote_rank)
        list(j+1) = list(j)
        j = j - 1
      end do
      list(j+1) = tmp
    end do
  end subroutine sort_entries_by_rank


  !> Build rank groups from a sorted entry list.
  subroutine build_rank_groups(list, n, groups, n_groups)
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(rank_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    ! Local
    integer :: i, ng, cur_rank

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count unique ranks
    ng = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%rank   = list(1)%remote_rank
    groups(1)%offset = 1
    groups(1)%count  = 1
    do i = 2, n
      if (list(i)%remote_rank /= list(i-1)%remote_rank) then
        ng = ng + 1
        groups(ng)%rank   = list(i)%remote_rank
        groups(ng)%offset = i
        groups(ng)%count  = 1
      else
        groups(ng)%count = groups(ng)%count + 1
      end if
    end do
  end subroutine build_rank_groups


  !> Build face groups from send entry list.
  !> Groups consecutive entries that share (local_block, src_f).
  subroutine build_face_groups_send(list, n, groups, n_groups)
    use ARES_Parameters_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%src_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%src_f, 1)
    groups(1)%dj = guide(list(1)%src_f, 2)
    groups(1)%dk = guide(list(1)%src_f, 3)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%src_f /= list(i-1)%src_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%src_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%src_f, 1)
        groups(ng)%dj = guide(list(i)%src_f, 2)
        groups(ng)%dk = guide(list(i)%src_f, 3)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_send


  !> Build face groups from recv entry list.
  !> Groups consecutive entries that share (local_block, dst_f).
  subroutine build_face_groups_recv(list, n, groups, n_groups)
    use ARES_Parameters_m, only: guide
    implicit none
    type(ghost_entry_type), intent(in) :: list(:)
    integer, intent(in) :: n
    type(face_group_type), allocatable, intent(out) :: groups(:)
    integer, intent(out) :: n_groups
    integer :: i, ng

    if (n == 0) then
      n_groups = 0
      allocate(groups(0))
      return
    end if

    ! Count face groups
    ng = 1
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) ng = ng + 1
    end do
    n_groups = ng

    allocate(groups(ng))
    ng = 1
    groups(1)%block = list(1)%local_block
    groups(1)%face  = list(1)%dst_f
    groups(1)%entry_start = 1
    groups(1)%entry_count = 1
    groups(1)%di = guide(list(1)%dst_f, 1)
    groups(1)%dj = guide(list(1)%dst_f, 2)
    groups(1)%dk = guide(list(1)%dst_f, 3)
    do i = 2, n
      if (list(i)%local_block /= list(i-1)%local_block .or. &
          list(i)%dst_f /= list(i-1)%dst_f) then
        ng = ng + 1
        groups(ng)%block = list(i)%local_block
        groups(ng)%face  = list(i)%dst_f
        groups(ng)%entry_start = i
        groups(ng)%entry_count = 1
        groups(ng)%di = guide(list(i)%dst_f, 1)
        groups(ng)%dj = guide(list(i)%dst_f, 2)
        groups(ng)%dk = guide(list(i)%dst_f, 3)
      else
        groups(ng)%entry_count = groups(ng)%entry_count + 1
      end if
    end do
  end subroutine build_face_groups_recv
#endif


  !> Post persistent P receives. Call from !$omp single.
  subroutine exchange_ghost_P_post_recv(domain)
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_post_recv()
#endif
  end subroutine exchange_ghost_P_post_recv


  !> Pack P send buffer. Call from !$omp do over face groups.
  !> n_send_faces is accessible via ghost_sched%n_send_faces.
  subroutine exchange_ghost_P_pack(domain, fg_start, fg_end)
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    integer, intent(in) :: fg_start, fg_end
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_pack_faces(domain, fg_start, fg_end)
#endif
  end subroutine exchange_ghost_P_pack


  !> Post persistent P sends. Call from !$omp single after pack is complete.
  subroutine exchange_ghost_P_post_send(domain)
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_post_send()
#endif
  end subroutine exchange_ghost_P_post_send


  !> Wait for P receives and unpack. Call wait from !$omp single, unpack from !$omp do.
  subroutine exchange_ghost_P_wait_unpack(domain)
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_wait_recv()
    call exchange_P_unpack_faces(domain, 1, ghost_sched%n_recv_faces)
#endif
  end subroutine exchange_ghost_P_wait_unpack


  !> Wait for P sends to complete. Call from !$omp single.
  subroutine exchange_ghost_P_wait_send(domain)
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_P_wait_send()
#endif
  end subroutine exchange_ghost_P_wait_send


  !> Exchange ghost cell R (residual) data between MPI ranks.
  subroutine exchange_ghost_R(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, nprim
    use ARES_Parameters_m, only: guide

    implicit none
    type(ARES_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_R_field_begin(domain)
    call exchange_R_field_end(domain)
#endif
  end subroutine exchange_ghost_R


  !> Gather all blocks' P interior data to rank 0 for I/O.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL for overlap.
  !> Temporarily allocates P on root for remote blocks if not already allocated.
  subroutine gather_P_to_root(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: nprim, gc

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, n1, n2, n3, ncells

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, r, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute per-block sizes and total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_ncells(b) = nprim * n1 * n2 * n3
        blk_offset(b) = total_buf_size
        if (mpi_is_root) then
          if (.not. is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        else
          if (is_local_block(b)) total_buf_size = total_buf_size + blk_ncells(b)
        end if
      end do
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! Recompute offsets contiguously
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do

      nreq = 0
      if (mpi_is_root) then
        ! Root: post all receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and post all sends
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%P(:, 1:n1, 1:n2, 1:n3), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data (allocate P if it was previously freed)
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%P)) &
              allocate(domain%blk(b)%P(nprim, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc))
            domain%blk(b)%P(:, 1:n1, 1:n2, 1:n3) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [nprim, n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine gather_P_to_root


  !> Gather diagnostic arrays (R, dtlocal, beta) from all ranks to root.
  !> On root, remote blocks may have had these arrays deallocated by
  !> deallocate_remote_computation_data, so they are re-allocated here
  !> to receive the data.  All ranks must call this (collective).
  subroutine gather_diagnostic_to_root(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: nprim

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    integer :: b, n1, n2, n3

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, offset, total_buf_size
      integer :: ncells_R, ncells_dt
      integer, allocatable :: reqs(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)
      ! Per-block: nprim*n1*n2*n3 for R, n1*n2*n3 for dtlocal
      integer, allocatable :: blk_offset(:), blk_size(:)

      allocate(blk_offset(domain%nb), blk_size(domain%nb))

      ! Compute buffer sizes
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_size(b) = nprim * n1 * n2 * n3 + n1 * n2 * n3  ! R + dtlocal
      end do

      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_size(b)
        else
          if (is_local_block(b)) offset = offset + blk_size(b)
        end if
      end do
      total_buf_size = offset
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      ! On root: ensure R, dtlocal are allocated for remote blocks
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            if (.not. allocated(domain%blk(b)%R)) &
              allocate(domain%blk(b)%R(nprim, n1, n2, n3))
            if (.not. allocated(domain%blk(b)%dtlocal)) &
              allocate(domain%blk(b)%dtlocal(n1, n2, n3))
          end if
        end do
      end if

      nreq = 0
      if (mpi_is_root) then
        ! Root: post receives for remote blocks
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: pack local blocks and send
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            ncells_R  = nprim * n1 * n2 * n3
            ncells_dt = n1 * n2 * n3
            offset = blk_offset(b)
            buf(offset+1 : offset+ncells_R) = &
              reshape(domain%blk(b)%R(:, 1:n1, 1:n2, 1:n3), [ncells_R])
            offset = offset + ncells_R
            buf(offset+1 : offset+ncells_dt) = &
              reshape(domain%blk(b)%dtlocal(1:n1, 1:n2, 1:n3), [ncells_dt])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_size(b), MPI_DOUBLE_PRECISION, &
                           0, 1000+b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Root: unpack received data
      if (mpi_is_root) then
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            ncells_R  = nprim * n1 * n2 * n3
            ncells_dt = n1 * n2 * n3
            offset = blk_offset(b)
            domain%blk(b)%R(:, 1:n1, 1:n2, 1:n3) = &
              reshape(buf(offset+1 : offset+ncells_R), [nprim, n1, n2, n3])
            offset = offset + ncells_R
            domain%blk(b)%dtlocal(1:n1, 1:n2, 1:n3) = &
              reshape(buf(offset+1 : offset+ncells_dt), [n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_size)
    end block
#endif
  end subroutine gather_diagnostic_to_root


  !> Scatter blocks' P data from rank 0 to owning ranks after reading.
  !> Uses non-blocking MPI_ISEND/IRECV + WAITALL for overlap.
  subroutine scatter_P_from_root(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: nprim

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, n1, n2, n3, ncells

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr, nreq, r, total_buf_size, offset
      integer, allocatable :: reqs(:), blk_offset(:), blk_ncells(:)
      integer, allocatable :: stats(:,:)
      real(R8), allocatable :: buf(:)

      ! Compute per-block sizes and total buffer size
      allocate(blk_offset(domain%nb), blk_ncells(domain%nb))
      total_buf_size = 0
      do b = 1, domain%nb
        n1 = domain%blk(b)%dim(1)
        n2 = domain%blk(b)%dim(2)
        n3 = domain%blk(b)%dim(3)
        blk_ncells(b) = nprim * n1 * n2 * n3
      end do

      ! Compute contiguous buffer offsets
      offset = 0
      do b = 1, domain%nb
        blk_offset(b) = offset
        if (mpi_is_root) then
          if (.not. is_local_block(b)) offset = offset + blk_ncells(b)
        else
          if (is_local_block(b)) offset = offset + blk_ncells(b)
        end if
      end do
      total_buf_size = offset
      allocate(buf(max(total_buf_size, 1)))
      allocate(reqs(domain%nb))
      allocate(stats(MPI_STATUS_SIZE, domain%nb))

      nreq = 0
      if (mpi_is_root) then
        ! Root: pack remote blocks and post all sends
        do b = 1, domain%nb
          if (.not. is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)) = &
              reshape(domain%blk(b)%P(:, 1:n1, 1:n2, 1:n3), [blk_ncells(b)])
            nreq = nreq + 1
            call MPI_ISEND(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           block_owner(b), b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      else
        ! Non-root: post all receives for local blocks
        do b = 1, domain%nb
          if (is_local_block(b)) then
            nreq = nreq + 1
            call MPI_IRECV(buf(blk_offset(b)+1), blk_ncells(b), MPI_DOUBLE_PRECISION, &
                           0, b, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          end if
        end do
      end if

      ! Wait for all
      if (nreq > 0) then
        call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
        call check_mpi_error(ierr)
      end if

      ! Non-root: unpack received data
      if (.not. mpi_is_root) then
        do b = 1, domain%nb
          if (is_local_block(b)) then
            n1 = domain%blk(b)%dim(1)
            n2 = domain%blk(b)%dim(2)
            n3 = domain%blk(b)%dim(3)
            domain%blk(b)%P(:, 1:n1, 1:n2, 1:n3) = &
              reshape(buf(blk_offset(b)+1 : blk_offset(b)+blk_ncells(b)), [nprim, n1, n2, n3])
          end if
        end do
      end if

      deallocate(buf, reqs, stats, blk_offset, blk_ncells)
    end block
#endif
  end subroutine scatter_P_from_root


  !> Exchange Pg(:,3:6) stencil data for inter-rank viscous connection BCs.
  !> Messages are aggregated per remote rank.
  subroutine exchange_ghost_Pg(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: nprim

    implicit none
    type(ARES_domain_type), intent(inout) :: domain

    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    call exchange_Pg_field(domain)
#endif
  end subroutine exchange_ghost_Pg


#ifdef USE_MPI
  !> Internal: exchange Pg(:,3:6) stencil using aggregated messaging.
  subroutine exchange_Pg_field(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: nprim
    use mpi

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: ns, nr, i, r, ierr, tag, bc_idx, buf_pos
    integer :: entry_size
    integer, parameter :: PG_NCOLS = 4   ! Pg columns 3..6

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns == 0 .and. nr == 0) return

    entry_size = nprim * PG_NCOLS

    ! Post one receive per remote rank (aggregated)
    do r = 1, ghost_sched%n_recv_ranks
      buf_pos = (ghost_sched%recv_groups(r)%offset - 1) * entry_size + 1
      tag = ghost_sched%recv_groups(r)%rank + 10000  ! unique tag per rank for Pg
      call MPI_IRECV(ghost_sched%Pg_recv_buf(buf_pos), &
                     ghost_sched%recv_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%recv_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%recv_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Pack all entries per rank into contiguous buffer, then one send per rank
    do r = 1, ghost_sched%n_send_ranks
      do i = ghost_sched%send_groups(r)%offset, &
             ghost_sched%send_groups(r)%offset + ghost_sched%send_groups(r)%count - 1
        bc_idx = ghost_sched%send_list(i)%bc_idx
        buf_pos = (i - 1) * entry_size
        ghost_sched%Pg_send_buf(buf_pos+1           : buf_pos+  nprim) = domain%bc(bc_idx)%Pg(:,3)
        ghost_sched%Pg_send_buf(buf_pos+  nprim+1   : buf_pos+2*nprim) = domain%bc(bc_idx)%Pg(:,4)
        ghost_sched%Pg_send_buf(buf_pos+2*nprim+1   : buf_pos+3*nprim) = domain%bc(bc_idx)%Pg(:,5)
        ghost_sched%Pg_send_buf(buf_pos+3*nprim+1   : buf_pos+4*nprim) = domain%bc(bc_idx)%Pg(:,6)
      end do
      buf_pos = (ghost_sched%send_groups(r)%offset - 1) * entry_size + 1
      tag = mpi_rank_ + 10000  ! sender's rank as tag for Pg
      call MPI_ISEND(ghost_sched%Pg_send_buf(buf_pos), &
                     ghost_sched%send_groups(r)%count * entry_size, &
                     MPI_DOUBLE_PRECISION, &
                     ghost_sched%send_groups(r)%rank, tag, &
                     MPI_COMM_WORLD, ghost_sched%send_req(r), ierr)
      call check_mpi_error(ierr)
    end do

    ! Wait for all receives
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%recv_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack into bc%Pg(:,3:6)
    do i = 1, nr
      bc_idx = ghost_sched%recv_list(i)%bc_idx
      buf_pos = (i - 1) * entry_size
      domain%bc(bc_idx)%Pg(:,3) = ghost_sched%Pg_recv_buf(buf_pos+1           : buf_pos+  nprim)
      domain%bc(bc_idx)%Pg(:,4) = ghost_sched%Pg_recv_buf(buf_pos+  nprim+1   : buf_pos+2*nprim)
      domain%bc(bc_idx)%Pg(:,5) = ghost_sched%Pg_recv_buf(buf_pos+2*nprim+1   : buf_pos+3*nprim)
      domain%bc(bc_idx)%Pg(:,6) = ghost_sched%Pg_recv_buf(buf_pos+3*nprim+1   : buf_pos+4*nprim)
    end do

    ! Wait for all sends to complete
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%send_req, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_Pg_field


  !> Post persistent receives for P field. Call from !$omp single.
  subroutine exchange_P_post_recv()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%P_recv_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_post_recv


  !> Pack P send buffer. Safe to call from !$omp do over face groups.
  subroutine exchange_P_pack_faces(domain, fg_start, fg_end)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, nprim
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    integer, intent(in) :: fg_start, fg_end
    ! Local
    integer :: fg, i, g, b, buf_pos
    integer :: Is, Js, Ks, di, dj, dk
    integer :: entry_size

    entry_size = nprim * gc

    do fg = fg_start, fg_end
      b  = ghost_sched%send_faces(fg)%block
      di = ghost_sched%send_faces(fg)%di
      dj = ghost_sched%send_faces(fg)%dj
      dk = ghost_sched%send_faces(fg)%dk
      do i = ghost_sched%send_faces(fg)%entry_start, &
             ghost_sched%send_faces(fg)%entry_start + ghost_sched%send_faces(fg)%entry_count - 1
        Is = ghost_sched%send_list(i)%src_i
        Js = ghost_sched%send_list(i)%src_j
        Ks = ghost_sched%send_list(i)%src_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          ghost_sched%P_send_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim) = &
            domain%blk(b)%P(:, Is+di*(g-1), Js+dj*(g-1), Ks+dk*(g-1))
        end do
      end do
    end do
  end subroutine exchange_P_pack_faces


  !> Post persistent sends for P field. Call from !$omp single after pack.
  subroutine exchange_P_post_send()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%P_send_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_post_send


  !> Wait for P receives and unpack into ghost cells. Safe to call unpack from !$omp do.
  subroutine exchange_P_wait_recv()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%P_recv_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_wait_recv


  !> Unpack P receive buffer. Safe to call from !$omp do over face groups.
  subroutine exchange_P_unpack_faces(domain, fg_start, fg_end)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, nprim
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    integer, intent(in) :: fg_start, fg_end
    ! Local
    integer :: fg, i, g, b, buf_pos
    integer :: Id, Jd, Kd, di, dj, dk
    integer :: entry_size

    entry_size = nprim * gc

    do fg = fg_start, fg_end
      b  = ghost_sched%recv_faces(fg)%block
      di = ghost_sched%recv_faces(fg)%di
      dj = ghost_sched%recv_faces(fg)%dj
      dk = ghost_sched%recv_faces(fg)%dk
      do i = ghost_sched%recv_faces(fg)%entry_start, &
             ghost_sched%recv_faces(fg)%entry_start + ghost_sched%recv_faces(fg)%entry_count - 1
        Id = ghost_sched%recv_list(i)%dst_i
        Jd = ghost_sched%recv_list(i)%dst_j
        Kd = ghost_sched%recv_list(i)%dst_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          domain%blk(b)%P(:, Id-di*g, Jd-dj*g, Kd-dk*g) = &
            ghost_sched%P_recv_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim)
        end do
      end do
    end do
  end subroutine exchange_P_unpack_faces


  !> Wait for P sends to complete. Call from !$omp single.
  subroutine exchange_P_wait_send()
    use mpi
    implicit none
    integer :: ierr

    if (ghost_sched%n_send == 0 .and. ghost_sched%n_recv == 0) return
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%P_send_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_P_wait_send


  !> Internal: begin non-blocking exchange of R field using persistent requests.
  subroutine exchange_R_field_begin(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, nprim
    use mpi

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: ns, nr, i, g, b, fg, ierr, buf_pos
    integer :: Is, Js, Ks, di, dj, dk
    integer :: entry_size

    ns = ghost_sched%n_send
    nr = ghost_sched%n_recv
    if (ns == 0 .and. nr == 0) return

    entry_size = nprim * gc

    ! Start persistent receives for R
    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_recv_ranks, ghost_sched%R_recv_req_pers, ierr)
      call check_mpi_error(ierr)
    end if

    ! Pack using face groups (reads R instead of P)
    do fg = 1, ghost_sched%n_send_faces
      b  = ghost_sched%send_faces(fg)%block
      di = ghost_sched%send_faces(fg)%di
      dj = ghost_sched%send_faces(fg)%dj
      dk = ghost_sched%send_faces(fg)%dk
      do i = ghost_sched%send_faces(fg)%entry_start, &
             ghost_sched%send_faces(fg)%entry_start + ghost_sched%send_faces(fg)%entry_count - 1
        Is = ghost_sched%send_list(i)%src_i
        Js = ghost_sched%send_list(i)%src_j
        Ks = ghost_sched%send_list(i)%src_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          ghost_sched%P_send_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim) = &
            domain%blk(b)%R(:, Is+di*(g-1), Js+dj*(g-1), Ks+dk*(g-1))
        end do
      end do
    end do

    ! Start persistent sends for R
    if (ghost_sched%n_send_ranks > 0) then
      call MPI_STARTALL(ghost_sched%n_send_ranks, ghost_sched%R_send_req_pers, ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_R_field_begin


  !> Internal: complete non-blocking exchange and unpack R into ghost cells.
  subroutine exchange_R_field_end(domain)
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, nprim
    use mpi

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: nr, i, g, b, fg, ierr, buf_pos
    integer :: Id, Jd, Kd, di, dj, dk
    integer :: entry_size

    nr = ghost_sched%n_recv
    if (ghost_sched%n_send == 0 .and. nr == 0) return

    entry_size = nprim * gc

    if (ghost_sched%n_recv_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_recv_ranks, ghost_sched%R_recv_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_recv_ranks), ierr)
      call check_mpi_error(ierr)
    end if

    ! Unpack using face groups
    do fg = 1, ghost_sched%n_recv_faces
      b  = ghost_sched%recv_faces(fg)%block
      di = ghost_sched%recv_faces(fg)%di
      dj = ghost_sched%recv_faces(fg)%dj
      dk = ghost_sched%recv_faces(fg)%dk
      do i = ghost_sched%recv_faces(fg)%entry_start, &
             ghost_sched%recv_faces(fg)%entry_start + ghost_sched%recv_faces(fg)%entry_count - 1
        Id = ghost_sched%recv_list(i)%dst_i
        Jd = ghost_sched%recv_list(i)%dst_j
        Kd = ghost_sched%recv_list(i)%dst_k
        buf_pos = (i - 1) * entry_size
        do g = 1, gc
          domain%blk(b)%R(:, Id-di*g, Jd-dj*g, Kd-dk*g) = &
            ghost_sched%P_recv_buf(buf_pos + (g-1)*nprim+1 : buf_pos + g*nprim)
        end do
      end do
    end do

    if (ghost_sched%n_send_ranks > 0) then
      call MPI_WAITALL(ghost_sched%n_send_ranks, ghost_sched%R_send_req_pers, &
                       ghost_sched%mpi_stat(:,1:ghost_sched%n_send_ranks), ierr)
      call check_mpi_error(ierr)
    end if
  end subroutine exchange_R_field_end
#endif


  !> MPI barrier for synchronizing ranks after I/O operations.
  !> Prevents non-root ranks from advancing to the next MPI collective
  !> while root is still writing files.
  subroutine mpi_io_barrier()
    if (mpi_size_ <= 1) return
#ifdef USE_MPI
    block
      use mpi
      integer :: ierr
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end block
#endif
  end subroutine mpi_io_barrier

  subroutine Ghost_Interrank( Im, Jm, Km, Fm, blkm, Pg )
    use ARES_Global_m, only: gc, nprim
    use ARES_Parameters_m, only: guide
    use ARES_Advanced_Types_m
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blkm
    real(R8), intent(inout)              :: Pg(nprim,6)
    integer                              :: g, Ig, Jg, Kg

    do g = 1, gc
        Ig = Im - guide(Fm,1)*g
        Jg = Jm - guide(Fm,2)*g
        Kg = Km - guide(Fm,3)*g
        Pg(:,g) = blkm % P(:,Ig,Jg,Kg)
    enddo
  end subroutine Ghost_Interrank


end module ARES_Mod_GhostExchange
