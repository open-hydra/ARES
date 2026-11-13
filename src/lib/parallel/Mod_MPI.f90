!>@brief MPI infrastructure module for ARES.
!> Provides MPI environment management, block partitioning, and
!> global reduction wrappers. When compiled without USE_MPI, all
!> routines are no-ops or identity operations (serial fallback).
module ARES_Mod_MPI
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I4 => int32, output_unit

  implicit none
  private

  ! --- Public state ---
  integer, public :: mpi_rank_ = 0        !< This process rank (0 in serial)
  integer, public :: mpi_size_ = 1        !< Total number of MPI processes
  logical, public :: mpi_is_root = .true.  !< True on rank 0

  ! --- Outcome of reading a decomposition map ---
  integer, parameter :: MAP_ABSENT = 0, MAP_READ = 1, MAP_ERROR = 2

  ! --- Block-to-rank mapping ---
  integer, allocatable, public :: block_owner(:)     !< block_owner(b) = rank owning block b
  integer, allocatable, public :: local_block_ids(:) !< Global block IDs owned by this rank
  integer, public :: n_local_blocks = 0

  ! --- Public procedures ---
  public :: mpi_init_env, mpi_finalize_env
  public :: is_local_block
  public :: partition_blocks
  public :: mpi_allreduce_sum_r8, mpi_allreduce_min_r8, mpi_allreduce_max_r8
  public :: mpi_gather_r8
  public :: mpi_allreduce_sum_r8_array
  public :: mpi_reduce_sum_r8, mpi_reduce_sum_r8_array
  public :: mpi_allreduce_norm2
  public :: mpi_bcast_logical, mpi_bcast_integer, mpi_bcast_chars
  public :: mpi_abort_all, mpi_stop_all
  public :: check_mpi_error

contains


  !> Initialize MPI environment with thread support. Call early in main program.
  !> Uses MPI_THREAD_FUNNELED: only master thread makes MPI calls.
  subroutine mpi_init_env()
#ifdef USE_MPI
    integer :: ierr, provided
    logical :: already

    call MPI_Initialized(already, ierr)
    if (.not. already) then
      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
      call check_mpi_error(ierr)
      if (provided < MPI_THREAD_FUNNELED) then
        write(*,'(A)') ' ERROR: MPI does not support the required threading level (MPI_THREAD_FUNNELED)'
        call flush_output_units()
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if
    end if
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank_, ierr)
    call check_mpi_error(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size_, ierr)
    call check_mpi_error(ierr)
    mpi_is_root = (mpi_rank_ == 0)
#else
    mpi_rank_ = 0
    mpi_size_ = 1
    mpi_is_root = .true.
#endif
  end subroutine mpi_init_env


  !> Finalize MPI environment. Call at end of main program.
  subroutine mpi_finalize_env()
#ifdef USE_MPI
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#endif
  end subroutine mpi_finalize_env


  !> Check if block b is owned by this rank.
  logical function is_local_block(b)
    integer, intent(in) :: b

    if (.not. allocated(block_owner)) then
      is_local_block = .true.  ! Serial mode or not yet partitioned
      return
    end if
    is_local_block = (block_owner(b) == mpi_rank_)
  end function is_local_block


  !> Partition nb blocks across MPI ranks.
  !> blk_ncells(b) = total number of cells in block b; blk_dims(:,b) = its cells
  !> along i, j, k.
  !> When map_file is given (the decomposition map that the key decomposition-map
  !> names), the rank of every block is read from it (read_decomposition_map); a
  !> file that cannot be read, or that does not fit the mesh or the run, stops the
  !> run. Otherwise the blocks are walked in the order of the mesh, each to the
  !> rank with the least load (the lowest such rank on a tie), and nothing is
  !> printed.
  subroutine partition_blocks(nb, blk_ncells, blk_dims, map_file)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    integer, intent(in), optional :: blk_dims(3,nb)
    character(len=*), intent(in), optional :: map_file
    ! Local
    integer :: b, r, nloc
    integer, allocatable :: rank_load(:), map_owner(:)
    logical :: from_map

    if (allocated(block_owner)) deallocate(block_owner)
    if (allocated(rank_load)) deallocate(rank_load)
    allocate(block_owner(nb))
    allocate(rank_load(0:mpi_size_-1))
    rank_load = 0

    ! Greedy: assign each block to the rank with the least load
    r = 0
    do b = 1, nb
      r = minloc(rank_load, dim=1) - 1   ! rank with minimum load (0-indexed)
      block_owner(b) = r
      rank_load(r) = rank_load(r) + blk_ncells(b)
    end do

    ! The assignment of the decomposition map, when one is named, replaces it
    from_map = .false.
    if (present(blk_dims) .and. present(map_file)) then
      allocate(map_owner(nb))
      call read_decomposition_map(map_file, nb, blk_dims, blk_ncells, map_owner, from_map)
      if (from_map) then
        block_owner = map_owner
        rank_load = 0
        do b = 1, nb
          rank_load(block_owner(b)) = rank_load(block_owner(b)) + blk_ncells(b)
        end do
        call report_partition_balance(nb, blk_ncells, rank_load, map_file)
      end if
      deallocate(map_owner)
    end if

    ! Build local block list
    nloc = count(block_owner == mpi_rank_)
    n_local_blocks = nloc
    if (allocated(local_block_ids)) deallocate(local_block_ids)
    allocate(local_block_ids(nloc))
    nloc = 0
    do b = 1, nb
      if (block_owner(b) == mpi_rank_) then
        nloc = nloc + 1
        local_block_ids(nloc) = b
      end if
    end do

    deallocate(rank_load)
  end subroutine partition_blocks


  !> Balance of a partition read from a decomposition map, printed by rank 0: the
  !> balance of the cells over the ranks, the file, and a warning about ranks the
  !> map leaves without blocks. The assignment that ARES makes itself (no map)
  !> prints nothing.
  subroutine report_partition_balance(nb, blk_ncells, rank_load, map_file)
    integer, intent(in) :: nb
    integer, intent(in) :: blk_ncells(nb)
    integer, intent(in) :: rank_load(0:)
    character(len=*), intent(in) :: map_file
    ! Local
    real(R8) :: ideal, eff

    if (mpi_rank_ /= 0) return
    if (mpi_size_ <= 1) then
      write(*,'(A,A)') '  MPI partition read from ', trim(map_file)
      return
    end if

    ideal = real(sum(blk_ncells), R8) / real(mpi_size_, R8)
    eff   = ideal / real(maxval(rank_load), R8) * 100.0_R8

    write(*,'(A,I0,A,I0,A,F5.1,A)') '  MPI partition: ', nb, ' blocks over ', &
      mpi_size_, ' ranks, balance ', eff, '% of ideal'
    write(*,'(A,A)') '  MPI partition read from ', trim(map_file)
    if (count(rank_load == 0) > 0 .and. mpi_size_ > nb) then
      write(*,'(A,I0,A,I0,A,I0,A)') '  WARNING: the decomposition map gives no block to ', count(rank_load == 0), &
        ' of ', mpi_size_, ' ranks; they idle (the mesh has only ', nb, ' blocks).'
    else if (count(rank_load == 0) > 0) then
      write(*,'(A,I0,A,I0,A)') '  WARNING: the decomposition map gives no block to ', count(rank_load == 0), &
        ' of ', mpi_size_, ' ranks; they idle.'
    else if (real(maxval(blk_ncells), R8) > ideal) then
      write(*,'(A,I0,A)') '  NOTE: balance is capped by the largest block (', &
        maxval(blk_ncells), ' cells); only a finer mesh split can improve it.'
    end if
  end subroutine report_partition_balance


  !> Block-to-rank assignment from a decomposition map, in the format of the map
  !> that the ATLAS mesh decomposer (MDB) writes: lines starting with '#' are
  !> comments; then one line per block, blocks 1 to nb in order, with the 13
  !> integers
  !>   new parent i0 i1 j0 j1 k0 k1 ni nj nk cells rank
  !> of which ARES uses new (the block), ni nj nk and cells (checked against the
  !> mesh) and rank (from 0). The first '#' line after the block lines ends them
  !> (MDB writes the face origins of the blocks after it). The ranks must lie
  !> between 0 and N-1; when the mesh has at least as many blocks as the run has
  !> ranks, the highest rank must be used too, or the map was written for another
  !> number of ranks, an error; with fewer blocks than ranks, some ranks idle
  !> whatever the map (as without one), with a warning (report_partition_balance).
  !> A directory in place of the file is an error too.
  !> Rank 0 reads the file and broadcasts what it
  !> found, so that every rank takes the same decision: the assignment in owner
  !> (found true), or an error, which every rank prints before every rank stops
  !> (mpi_stop_all, exit code 1). A file that does not exist is an error: the
  !> map is read only when the key decomposition-map names it.
  subroutine read_decomposition_map(path, nb, blk_dims, blk_ncells, owner, found)
    character(len=*), intent(in) :: path
    integer, intent(in)  :: nb
    integer, intent(in)  :: blk_dims(3,nb), blk_ncells(nb)
    integer, intent(out) :: owner(nb)
    logical, intent(out) :: found
    ! Local
    integer :: status
    character(len=256) :: message
#ifdef USE_MPI
    integer :: ierr
#endif

    owner   = -1
    status  = MAP_ABSENT
    message = ''
    if (mpi_rank_ == 0) call parse_decomposition_map(path, nb, blk_dims, blk_ncells, owner, status, message)
#ifdef USE_MPI
    if (mpi_size_ > 1) then
      call MPI_BCAST(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
      call MPI_BCAST(owner, nb, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
      call MPI_BCAST(message, len(message), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end if
#endif
    if (status == MAP_ERROR) then
      write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] ERROR: ', trim(message)
      flush(output_unit)
#ifdef USE_MPI
      if (mpi_size_ > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      ! Every rank has the error, broadcast by rank 0: the ranks stop together, so that the output of every rank reaches the
      ! log, and not with MPI_ABORT, which only makes a best attempt to stop the processes
      call mpi_stop_all('decomposition map error (see above)')
    end if
    found = (status == MAP_READ)
  end subroutine read_decomposition_map


  !> Internal: read and check the decomposition map on rank 0 (read_decomposition_map).
  subroutine parse_decomposition_map(path, nb, blk_dims, blk_ncells, owner, status, message)
    character(len=*), intent(in)    :: path
    integer, intent(in)             :: nb
    integer, intent(in)             :: blk_dims(3,nb), blk_ncells(nb)
    integer, intent(inout)          :: owner(nb)
    integer, intent(out)            :: status
    character(len=*), intent(inout) :: message
    ! Local
    integer, parameter :: UNSET = -huge(1)
    character(len=1024) :: line
    integer :: unit, ios, nline, nrow, v(13), b
    logical :: exists, is_dir, in_rows

    ! The map is read only when the key decomposition-map names it: a file that does
    ! not exist is an error (it cannot be opened, below), and so is a directory of
    ! that name. gfortran reports a directory as existing (and path/. too); ifx does
    ! not, but opens it.
    status = MAP_ERROR
    inquire(file=path, exist=exists)
    inquire(file=trim(path)//'/.', exist=is_dir)
    if (.not. exists .and. .not. is_dir) then
      open(newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios == 0) then
        close(unit)
        is_dir = .true.
      end if
    end if
    if (is_dir) then
      message = trim(path)//': is a directory, not a file'
      return
    end if
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      message = trim(path)//': cannot be opened'
      return
    end if

    nline = 0
    nrow = 0
    in_rows = .false.
    do
      read(unit, '(A)', iostat=ios) line
      if (is_iostat_end(ios)) exit
      if (ios /= 0) then
        write(message, '(A,A,I0,A,I0,A)') trim(path), ': read error after line ', nline, ' (iostat ', ios, ')'
        close(unit)
        return
      end if
      nline = nline + 1
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') then
        if (in_rows) exit          ! end of the block lines
        cycle
      end if
      in_rows = .true.
      v = UNSET
      ios = 1
      if (starts_with_13_integers(line)) read(line, *, iostat=ios) v
      if (ios /= 0 .or. any(v == UNSET)) then
        write(message, '(A,A,I0,A)') trim(path), ': line ', nline, &
          ' is not 13 integers (new parent i0 i1 j0 j1 k0 k1 ni nj nk cells rank)'
        close(unit)
        return
      end if
      nrow = nrow + 1
      if (nrow > nb) then
        write(message, '(A,A,I0,A)') trim(path), ': more block lines than the ', nb, ' blocks of the mesh'
        close(unit)
        return
      end if
      b = v(1)
      if (b /= nrow) then
        write(message, '(A,A,I0,A,I0,A,I0,A,I0,A)') trim(path), ': line ', nline, ' lists block ', b, &
          ' where block ', nrow, ' is expected (blocks 1 to ', nb, ', in order)'
        close(unit)
        return
      end if
      if (any(v(9:11) /= blk_dims(:,b))) then
        write(message, '(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') trim(path), ': block ', b, ' is ', &
          v(9), ' x ', v(10), ' x ', v(11), ' cells in the map, ', &
          blk_dims(1,b), ' x ', blk_dims(2,b), ' x ', blk_dims(3,b), ' in the mesh'
        close(unit)
        return
      end if
      if (v(12) /= blk_ncells(b)) then
        write(message, '(A,A,I0,A,I0,A,I0,A)') trim(path), ': block ', b, ' has ', v(12), &
          ' cells in the map, ', blk_ncells(b), ' in the mesh'
        close(unit)
        return
      end if
      if (v(13) < 0 .or. v(13) >= mpi_size_) then
        write(message, '(A,A,I0,A,I0,A,I0,A)') trim(path), ': block ', b, ' is on rank ', v(13), &
          ', outside the ranks of the run (0 to ', mpi_size_ - 1, ')'
        close(unit)
        return
      end if
      owner(b) = v(13)
    end do
    close(unit)

    if (nrow /= nb) then
      write(message, '(A,A,I0,A,I0,A)') trim(path), ': ', nrow, ' block lines for the ', nb, ' blocks of the mesh'
      return
    end if
    if (maxval(owner) + 1 < mpi_size_ .and. nb >= mpi_size_) then
      write(message, '(A,A,I0,A,I0,A,A,I0,A,I0,A)') trim(path), ': the largest rank is ', maxval(owner), &
        ': the map is for ', maxval(owner) + 1, trim(merge(' rank, ', ' ranks,', maxval(owner) == 0)), &
        ' but the run has ', mpi_size_, ' ranks and the mesh ', nb, &
        ' blocks (not fewer than the ranks): a map for another number of ranks'
      return
    end if
    status = MAP_READ
  end subroutine parse_decomposition_map


  !> Internal: whether line starts with 13 integers, each an optional sign and
  !> digits, separated by blanks, tabs or a carriage return. What follows the
  !> 13th is not looked at. A list-directed read alone may take 8.5 for 8.
  logical function starts_with_13_integers(line) result(ok)
    character(len=*), intent(in) :: line
    ! Local
    integer :: i, ntok, ndig
    character :: c

    ok = .false.
    ntok = 0
    ndig = -1                                 ! -1 between tokens, else the digits of the token
    do i = 1, len(line)
      c = line(i:i)
      if (c == ' ' .or. c == achar(9) .or. c == achar(13)) then
        if (ndig == 0) return                 ! a sign alone
        if (ndig > 0 .and. ntok == 13) then
          ok = .true.
          return
        end if
        ndig = -1
      else if (ndig < 0) then                 ! first character of a token
        ntok = ntok + 1
        if (c == '+' .or. c == '-') then
          ndig = 0
        else if (c >= '0' .and. c <= '9') then
          ndig = 1
        else
          return
        end if
      else if (c >= '0' .and. c <= '9') then
        ndig = ndig + 1
      else
        return
      end if
    end do
    ok = (ntok == 13 .and. ndig > 0)
  end function starts_with_13_integers


  !> MPI_ALLREDUCE with MPI_SUM for a scalar real(R8).
  subroutine mpi_allreduce_sum_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_sum_r8


  !> MPI_ALLREDUCE with MPI_MIN for a scalar real(R8).
  subroutine mpi_allreduce_min_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_min_r8


  !> MPI_ALLREDUCE with MPI_MAX for a scalar real(R8).
  subroutine mpi_allreduce_max_r8(local_val, global_val)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: global_val
#ifdef USE_MPI
    integer :: ierr
    call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    global_val = local_val
#endif
  end subroutine mpi_allreduce_max_r8


  !> Gather one real(R8) per rank onto root. `arr` must be at least mpi_size_
  !> long on rank 0; only rank 0's copy is meaningful on return.
  subroutine mpi_gather_r8(local_val, arr)
    real(R8), intent(in)  :: local_val
    real(R8), intent(out) :: arr(:)
#ifdef USE_MPI
    integer :: ierr
    call MPI_GATHER(local_val, 1, MPI_DOUBLE_PRECISION, &
                    arr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#else
    arr(1) = local_val
#endif
  end subroutine mpi_gather_r8


  !> MPI_ALLREDUCE with MPI_SUM for an array of real(R8): arr is overwritten with the global result.
  !> The send buffer is a copy of arr, not MPI_IN_PLACE.
  !> With gfortran the MPI sentinels in COMMON blocks (MPI_IN_PLACE, MPI_BOTTOM, MPI_STATUS_IGNORE,
  !> MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE, MPI_ARGV_NULL, MPI_ARGVS_NULL) are not recognised while
  !> -fno-underscoring is set (cmake/SetFortranFlags.cmake): do not use them.
  subroutine mpi_allreduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    real(R8), allocatable :: send(:)
    send = arr
    call MPI_ALLREDUCE(send, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_allreduce_sum_r8_array


  !> Distributed L2 norm: local sum of squares + MPI_ALLREDUCE SUM + sqrt.
  function mpi_allreduce_norm2(x, n) result(global_norm)
    integer, intent(in) :: n
    real(R8), intent(in) :: x(n)
    real(R8) :: global_norm
    real(R8) :: local_ss

    local_ss = dot_product(x, x)
#ifdef USE_MPI
    call mpi_allreduce_sum_r8(local_ss, global_norm)
    global_norm = sqrt(global_norm)
#else
    global_norm = sqrt(local_ss)
#endif
  end function mpi_allreduce_norm2


  !> MPI_REDUCE with MPI_SUM for a scalar real(R8) to root rank 0: val holds the sum on the root and is
  !> left as it is elsewhere, where the receive buffer is not significant. Every rank sends a copy of val:
  !> no buffer passed twice (an OUT argument cannot be aliased with another one) and no MPI_IN_PLACE.
  !> With gfortran the MPI sentinels in COMMON blocks (MPI_IN_PLACE, MPI_BOTTOM, MPI_STATUS_IGNORE,
  !> MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE, MPI_ARGV_NULL, MPI_ARGVS_NULL) are not recognised while
  !> -fno-underscoring is set (cmake/SetFortranFlags.cmake): do not use them.
  subroutine mpi_reduce_sum_r8(val)
    real(R8), intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    real(R8) :: send
    send = val
    call MPI_REDUCE(send, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8


  !> MPI_REDUCE with MPI_SUM for an array of real(R8) to root rank 0.
  !> Result is only valid on root. More efficient than ALLREDUCE when
  !> only root needs the result (e.g. for residual output). arr is left as it is elsewhere, where the
  !> receive buffer is not significant. Every rank sends a copy of arr: no buffer passed twice (an OUT
  !> argument cannot be aliased with another one) and no MPI_IN_PLACE.
  !> With gfortran the MPI sentinels in COMMON blocks (MPI_IN_PLACE, MPI_BOTTOM, MPI_STATUS_IGNORE,
  !> MPI_STATUSES_IGNORE, MPI_ERRCODES_IGNORE, MPI_ARGV_NULL, MPI_ARGVS_NULL) are not recognised while
  !> -fno-underscoring is set (cmake/SetFortranFlags.cmake): do not use them.
  subroutine mpi_reduce_sum_r8_array(arr, n)
    integer, intent(in)     :: n
    real(R8), intent(inout) :: arr(n)
#ifdef USE_MPI
    integer :: ierr
    real(R8), allocatable :: send(:)
    send = arr
    call MPI_REDUCE(send, arr, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_reduce_sum_r8_array


  !> Broadcast a logical value from root (rank 0) to all ranks.
  subroutine mpi_bcast_logical(val)
    logical, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_logical


  !> Broadcast an integer value from root (rank 0) to all ranks.
  subroutine mpi_bcast_integer(val)
    integer, intent(inout) :: val
#ifdef USE_MPI
    integer :: ierr
    call MPI_BCAST(val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine mpi_bcast_integer


  !> Broadcast a character string from root (rank 0) to all ranks: its length first, then its characters, so that every
  !> rank passes MPI_BCAST the count of root (MPI-4.1, section 6.4: the type signature must be the same on every
  !> process). The other ranks receive the string allocated to that length.
  subroutine mpi_bcast_chars(text)
    character(len=:), allocatable, intent(inout) :: text
    integer :: n
#ifdef USE_MPI
    integer :: ierr
#endif
    n = 0
    if (mpi_is_root .and. allocated(text)) n = len(text)
    call mpi_bcast_integer(n)
    if (.not. mpi_is_root) then
      if (allocated(text)) deallocate(text)
      allocate(character(len=n) :: text)
    end if
#ifdef USE_MPI
    if (n > 0) then
      call MPI_BCAST(text, n, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call check_mpi_error(ierr)
    end if
#endif
  end subroutine mpi_bcast_chars


  !> Check MPI return code. Abort all ranks on error.
  subroutine check_mpi_error(ierr)
    integer, intent(in) :: ierr
#ifdef USE_MPI
    integer :: abort_ierr
    if (ierr /= 0) then
      write(*,'(A,I4,A,I6)') ' [RANK ', mpi_rank_, '] MPI error code: ', ierr
      call flush_output_units()
      call MPI_ABORT(MPI_COMM_WORLD, 1, abort_ierr)
    end if
#endif
  end subroutine check_mpi_error


  !> MPI-aware abort: prints message, writes out the output of this process (flush_output_units), then aborts all
  !> ranks.
  subroutine mpi_abort_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] ABORT: ', trim(message)
    call flush_output_units()
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#else
    write(*,'(A,A)') ' ABORT: ', trim(message)
    error stop
#endif
  end subroutine mpi_abort_all


#ifdef USE_MPI
  !> Internal: write out what this process still holds for output_unit and error_unit, before MPI_ABORT. MPI_ABORT
  !> makes only a "best attempt" to abort the processes (MPI-4.1, section 11.4.2) and need not end this one through
  !> the end of the Fortran program, where the run-time library writes out its units: built with gfortran and run
  !> with Open MPI 4.1.4 as one process with its output on a file, a run without this call loses the message, what it
  !> wrote before and the lines on error_unit (where FLINT writes its warnings). Called before every MPI_ABORT of ARES. It
  !> writes out the units of this process only: what the other ranks still hold when the run is aborted is not
  !> written out by this call.
  subroutine flush_output_units()
    use iso_fortran_env, only: error_unit
    flush(output_unit)
    flush(error_unit)
  end subroutine flush_output_units
#endif


  !> Stop of every rank on an error that every rank has found, a decision taken
  !> on all ranks (the outcome of reading the decomposition map, of a check of
  !> rank 0 or of the reading of input.ini, which rank 0 broadcasts; a file that
  !> every rank reads, such as bc.txt):
  !> each rank prints the message and flushes its output, the ranks wait for
  !> each other, so that what any rank printed before (rank 0 the explanation of
  !> the error) is out before any rank ends, then the ranks finalize MPI and end
  !> with exit code 1. MPI_ABORT
  !> (mpi_abort_all) makes only a "best attempt" to abort the processes
  !> (MPI-4.1, section 11.4.2): with Intel MPI 2021.12 the launcher at times
  !> ends before writing out what the ranks printed (measured), and the first
  !> rank to call it can stop another one before that one prints. A program that
  !> ends normally calls MPI_FINALIZE in every process, which can then give an
  !> exit code (MPI-4.1, section 11.2.2).
  !> Precondition: every rank must reach this call on the same error, from the
  !> same data (files that every rank reads at set-up, the outcome of a
  !> broadcast). A decision that can differ between ranks, such as one on the
  !> initial condition, which each rank reads itself (Check_IC in Wrap_Setup),
  !> must be reduced first (mpi_allreduce_max_r8 of a flag); else the ranks that
  !> do not stop wait in their next collective call and the run hangs.
  subroutine mpi_stop_all(message)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer :: ierr
    write(*,'(A,I4,A,A)') ' [RANK ', mpi_rank_, '] STOP: ', trim(message)
    flush(output_unit)
    ! MPI_BARRIER returns at a rank only after every rank has entered it (MPI-4.1, section 6.3)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
#else
    write(*,'(A,A)') ' STOP: ', trim(message)
    flush(output_unit)
#endif
    error stop 1
  end subroutine mpi_stop_all


end module ARES_Mod_MPI
