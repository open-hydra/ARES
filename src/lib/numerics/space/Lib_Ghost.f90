module ARES_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Parameters_m
  use ARES_Global_m

  implicit none

  !> Ghost fills of a wall or 400 face (Ghost_Wall) in which the cut of an extrapolated quantity acted, in which the
  !> temperature of a 302 wall took the floor T_w/2, and in which the pressure or the temperature of a ghost cell was brought
  !> back within the fluid tables: reported at the end of a run in which any of them happened (Report_Wall_Ghosts)
  integer(I8) :: n_ghost_cut = 0, n_ghost_floor = 0, n_ghost_table = 0

contains

  subroutine Fill_Ghost_Cell ( domain )
    use ARES_Mod_MPI, only: is_local_block, mpi_size_
    use ARES_Mod_GhostExchange, only: exchange_ghost_P_post_recv, exchange_ghost_P_pack, &
                                       exchange_ghost_P_post_send, exchange_ghost_P_wait_unpack, &
                                       exchange_ghost_P_wait_send, &
                                       exchange_ghost_chimera_post_recv, exchange_ghost_chimera_pack, &
                                       exchange_ghost_chimera_post_send, exchange_ghost_chimera_wait_recv, &
                                       exchange_ghost_chimera_unpack, exchange_ghost_chimera_wait_send, &
                                       Ghost_Interrank, exchange_ghost_Pg, ghost_sched, &
                                       set_active_mg_level
    use ARES_Mod_Timers, only: timer_comm_begin, timer_comm_end
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    integer :: fg
    logical :: mpi_on

    ! On a single rank the exchange has nothing to do, and running it anyway
    ! costs a chain of empty master blocks, each one followed by a barrier.  The
    ! test is uniform over the team, so the worksharing constructs and barriers
    ! it guards are still encountered by every thread.
    mpi_on = mpi_size_ > 1

    if (mpi_on) then
      ! MPI: post persistent receives, pack buffer in parallel, then post sends
      ! Chimera donor cells travel in their own non-blocking exchange, started
      ! here so it overlaps with the local BC processing below.
      ! MPI is initialised with MPI_THREAD_FUNNELED: every MPI call is made by
      ! the primary thread, in a master construct (named masked since OpenMP
      ! 5.1; master is kept because older compilers lack masked), followed by
      ! an explicit barrier wherever the other threads need its result. A
      ! single construct would let any thread make the call, which FUNNELED
      ! forbids.
      !$omp master
      call set_active_mg_level(domain%mg_level)
      call exchange_ghost_P_post_recv(domain)
      call exchange_ghost_chimera_post_recv(domain)
      !$omp end master
      !$omp barrier

      ! Pack send buffers in parallel: P over face groups, chimera over donor cells
      !$omp do schedule(static) private(fg)
      do fg = 1, ghost_sched%n_send_faces
        call exchange_ghost_P_pack(domain, fg, fg)
      end do
      !$omp end do nowait
      !$omp do schedule(static) private(ii)
      do ii = 1, ghost_sched%n_chim_send
        call exchange_ghost_chimera_pack(domain, ii, ii)
      end do

      ! Post sends (must wait for all packing to complete — implicit barrier from !$omp do)
      !$omp master
      call exchange_ghost_P_post_send(domain)
      call exchange_ghost_chimera_post_send(domain)
      !$omp end master
    end if

    ! Process LOCAL BC entries while MPI communication is in flight
    ! Uses pre-filtered local_bc_idx to avoid scanning all nbound entries
    !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs)
    do ii = 1, domain % n_local_bc
      i  = domain % local_bc_idx(ii)
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type)
        case(101,201) ! block connection
          Bs = domain % bc(i) % bs
          if (.not. is_local_block(Bs)) cycle  ! inter-rank handled after MPI completes
          Is = domain % bc(i) % is
          Js = domain % bc(i) % js
          Ks = domain % bc(i) % ks
          Fs = domain % bc(i) % fs
          call Ghost_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                  Is, Js, Ks, Fs, domain % blk(Bs), &
                                  domain % bc(i) % Pg )
        case(200)
          call Ghost_Axysimmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
        case(300)
          call Ghost_Symmetry ( Im, Jm, Km, Fm, domain % blk(Bm) )
        case(0,404:406) ! inlet/outlet/extrapolation: zero-gradient
          call Ghost_ZG_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
        case(103,104,301,302,400) ! walls (103, 104: coupled to another solver) and extrapolation faces
          call Ghost_Wall ( domain % bc(i), domain % blk(Bm) )
        case(102)
          ! Chimera: processed after the donor-cell MPI exchange completes (below)
          continue
        case default
          call Ghost_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
      end select
    enddo

    ! MPI: wait for P receives to complete; unpack chimera donor cells into
    ! the kept-alive P arrays of remote donor blocks.  What is timed here is
    ! the communication the local BC work did not hide, plus the wait on
    ! slower neighbours.
    if (mpi_on) then
      !$omp master
      call timer_comm_begin()
      call exchange_ghost_P_wait_unpack(domain)
      call exchange_ghost_chimera_wait_recv(domain)
      call timer_comm_end()
      !$omp end master
      !$omp barrier

      ! Chimera donor cells are distinct, so the unpack has no write conflicts
      !$omp do schedule(static) private(ii)
      do ii = 1, ghost_sched%n_chim_recv
        call exchange_ghost_chimera_unpack(domain, ii, ii)
      end do
    end if

    ! Compute Pg(:,3:6) for type-1 connections where source block is local.
    ! The four tangential neighbours of a source cell on the edge of its face are
    ! ghost cells of Bs, and those can be filled from another rank: this loop must
    ! follow the P unpack above, or it reads them one stage old (MPI-only drift at
    ! block edges and junctions).
    !$omp do schedule (dynamic) private(ii, i, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do ii = 1, domain % n_local_bs
      i  = domain % local_bs_idx(ii)
      Bs = domain % bc(i) % bs
      Is = domain % bc(i) % is
      Js = domain % bc(i) % js
      Ks = domain % bc(i) % ks
      Fs = domain % bc(i) % fs
      d11s = domain % bc(i) % d11
      d12s = domain % bc(i) % d12
      d21s = domain % bc(i) % d21
      d22s = domain % bc(i) % d22
      call Fill_BC_Ghost_Connection ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, &
                                      domain % blk(Bs), domain % bc(i) % Pg )
    enddo

    ! Chimera ghost fill: all donor data (local and remote) is now current
    !$omp do schedule (dynamic) private(ii, i)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 102) cycle
      call Ghost_Chimera ( domain % nb, domain % blk, domain % bc(i) )
    enddo

    ! Chimera Pg(:,3:6): reads the ghost columns written just above, so it
    ! must stay a separate loop (implicit barrier in between). Each entry
    ! writes its own Pg only: beyond the edges of a block face the ghost row
    ! is extrapolated into Pg, not written in the block (Fill_BC_Ghost_Chimera),
    ! so the entries of two chimera faces of a block that meet at an edge write
    ! no common location, and one loop over all the entries has no data race
    !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 102) cycle
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      call Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Pg )
    enddo

    ! Process INTER-RANK type-1 entries (Bm local, Bs remote)
    if (mpi_on) then
      !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm, Bs)
      do ii = 1, domain % n_local_bc
        i  = domain % local_bc_idx(ii)
        if (domain % bc(i) % type /= 101 .and. domain % bc(i) % type /= 201) cycle
        Bs = domain % bc(i) % bs
        if (is_local_block(Bs)) cycle  ! already processed above
        Bm = domain % bc(i) % b
        Im = domain % bc(i) % i
        Jm = domain % bc(i) % j
        Km = domain % bc(i) % k
        Fm = domain % bc(i) % f
        call Ghost_Interrank(Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Pg)
      enddo

      ! Wait for P sends to complete before reusing buffers
      !$omp master
      call exchange_ghost_P_wait_send(domain)
      call exchange_ghost_chimera_wait_send(domain)
      call exchange_ghost_Pg(domain)
      !$omp end master
      !$omp barrier
    end if

  end subroutine Fill_Ghost_Cell


  subroutine Ghost_Connection ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, Pg )

    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm, Is, Js, Ks, Fs
    type(ARES_block_type), intent(in) :: blks
    type(ARES_block_type), intent(inout) :: blkm
    real(R8), intent(inout)              :: Pg(nprim,6)
    ! Local
    integer :: g, Ig, Jg, Kg

    do g = 1, gc
      Ig = Im - guide(Fm,1)*(g)
      Jg = Jm - guide(Fm,2)*(g)
      Kg = Km - guide(Fm,3)*(g)
      blkm % P(:,Ig,Jg,Kg) = blks % P (:, Is + guide(Fs,1)*(g-1), &
                                          Js + guide(Fs,2)*(g-1), &
                                          Ks + guide(Fs,3)*(g-1)  )
      Pg (:,g) = blkm % P (:,Ig,Jg,Kg)
    enddo

  end subroutine Ghost_Connection


  subroutine Ghost_Axysimmetry ( Im, Jm, Km, Fm, blk )
    use ARES_Lib_Metrics, only: cos_delthe, sin_delthe
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blk
    ! Local
    integer :: Ig, Jg, Kg, Ks

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    if (ndir < 3) then
      ! Rotate v and w of +/- delthe and copy other variables
      blk % P (nv,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Km) * cos_delthe + guide(Fm,3) * blk % P (nw,Im,Jm,Km) * sin_delthe
      blk % P (nw,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Km) * sin_delthe * (-guide(Fm,3)) + blk % P (nw,Im,Jm,Km) * cos_delthe
      blk % P (1:nu,Ig,Jg,Kg) = blk % P (1:nu,Im,Jm,Km)
      blk % P (nh:nprim,Ig,Jg,Kg) = blk % P (nh:nprim,Im,Jm,Km)
    else
      ! The cell next to the opposite k face, turned about x by -guide(Fm,3)*delthe: back by the
      ! sector angle behind face 5, forward behind face 6 (the same turn as the 2-D branch above)
      Ks = blk % dim(3) - Km + 1
      blk % P (nv,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Ks) * cos_delthe + guide(Fm,3) * blk % P (nw,Im,Jm,Ks) * sin_delthe
      blk % P (nw,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Ks) * sin_delthe * (-guide(Fm,3)) + blk % P (nw,Im,Jm,Ks) * cos_delthe
      blk % P (1:nu,Ig,Jg,Kg) = blk % P (1:nu,Im,Jm,Ks)
      blk % P (nh:nprim,Ig,Jg,Kg) = blk % P (nh:nprim,Im,Jm,Ks)
    endif
      
  end subroutine Ghost_Axysimmetry


  subroutine Ghost_Symmetry ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blk
    ! Local
    integer  :: Ig, Jg, Kg, Ic, Jc, Kc
    real(R8) :: normal(3), vel(3), udotn

    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    if (Fm <= 2) then
      Ic = Im - mod(Fm,2)
      normal = blk % dir(1) % f (Ic,Jm,Km) % N
    elseif (Fm <= 4) then
      Jc = Jm - mod(Fm,2)
      normal = blk % dir(2) % f (Im,Jc,Km) % N
    else
      Kc = Km - mod(Fm,2)
      normal = blk % dir(3) % f (Im,Jm,Kc) % N
    endif

    ! Formula to invert normal velocity: u = |u_n|n + |u_t|t -2*|u_n|n = -|u_n|n + |u_t|t
    vel = blk % P (nu:nw,Im,Jm,Km)   ! cell velocity
    udotn = dot_product( vel, normal )  ! normal component
    vel = vel - 2d0*udotn*normal         ! ghost cell velocity
    
    blk % P (nu:nw,Ig,Jg,Kg) = vel
    blk % P (np,Ig,Jg,Kg) = blk % P (np,Im,Jm,Km)
    blk % P (nh:nprim,Ig,Jg,Kg) = blk % P (nh:nprim,Im,Jm,Km)
      
  end subroutine Ghost_Symmetry


  subroutine Ghost_ZG_Extrapolate ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blk
    ! Local
    integer :: g

    ! Simple <zero gradient> extrapolation
    do g = 1, gc
      select case (Fm)
        case(1)
          blk % P(:,Im-g,Jm,Km) = blk % P (:,Im,Jm,Km)
        case(2)
          blk % P(:,Im+g,Jm,Km) = blk % P (:,Im,Jm,Km)
        case(3)
          blk % P(:,Im,Jm-g,Km) = blk % P (:,Im,Jm,Km)
        case(4)
          blk % P(:,Im,Jm+g,Km) = blk % P (:,Im,Jm,Km)
        case(5)
          blk % P(:,Im,Jm,Km-g) = blk % P (:,Im,Jm,Km)
        case(6)
          blk % P(:,Im,Jm,Km+g) = blk % P (:,Im,Jm,Km)
      end select
    enddo
      
  end subroutine Ghost_ZG_Extrapolate


  subroutine Ghost_Chimera ( nb, Blk, bc )
    implicit none
    integer, intent(in)                  :: nb
    type(ARES_block_type), intent(inout) :: Blk(nb)
    type(ARES_bc_type), intent(inout)    :: bc
    ! Local
    integer  :: Bm, Im, Jm, Km, Fm, Ig, Jg, Kg, Bs, Is, Js, Ks, c, c1, c2, g
    real(R8) :: primg(nprim)

    ! Preliminary assignments
    Bm = bc % b
    Im = bc % i
    Jm = bc % j
    Km = bc % k
    Fm = bc % f

    ! Ghost cells filled with the primitives { p vel h rho*r } of the donor cells,
    ! blended with the volume fractions: no conversion to conservative variables.
    ! Donors ni(1) for the first ghost row, ni(2) for the second.
    do g = 1, 2
      if (g == 1) then
        c1 = 1
        c2 = bc % ni(1)
      else
        c1 = bc % ni(1) + 1
        c2 = sum ( bc % ni )
      endif

      primg = 0.d0
      do c = c1, c2
        Bs = bc % donorID(c,1)
        Is = bc % donorID(c,2)
        Js = bc % donorID(c,3)
        Ks = bc % donorID(c,4)
        primg = primg + blk(Bs) % P (:,Is,Js,Ks) * bc % volume_fraction(c)
      enddo

      Ig = Im - guide(Fm,1)*g
      Jg = Jm - guide(Fm,2)*g
      Kg = Km - guide(Fm,3)*g
      blk(Bm) % P (:,Ig,Jg,Kg) = primg
      bc % Pg (:,g) = primg
    enddo

  end subroutine Ghost_Chimera


  subroutine Ghost_Extrapolate ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blk
    ! Local
    integer :: g

    ! Extrapolation with 2nd order accuracy. Theory: suppose we have a stencil x1,x2,x3,x4 with x4 unknown.
    ! Taylor 2nd order: x4 = x3 + f1(x3) + 1/2 f2(x3) with f1 and f2 1st and 2nd derivative. Approximate 
    ! f1(x3) = (x4 - x2)/2 central and f2(x3) = x3 - 2x2 + x1 left sided => x4 = 3*x3 - 3*x2 + x1
    do g = 1, gc
      select case (Fm)
        case(1)
          blk % P(:,Im-g,Jm,Km) = 3d0 * blk % P (:,Im-g+1,Jm,Km) - &
                                  3d0 * blk % P (:,Im-g+2,Jm,Km) + &
                                        blk % P (:,Im-g+3,Jm,Km)
        case(2)
          blk % P(:,Im+g,Jm,Km) = 3d0 * blk % P (:,Im+g-1,Jm,Km) - &
                                  3d0 * blk % P (:,Im+g-2,Jm,Km) + &
                                        blk % P (:,Im+g-3,Jm,Km)
        case(3)
          blk % P(:,Im,Jm-g,Km) = 3d0 * blk % P (:,Im,Jm-g+1,Km) - &
                                  3d0 * blk % P (:,Im,Jm-g+2,Km) + &
                                        blk % P (:,Im,Jm-g+3,Km)
        case(4)
          blk % P(:,Im,Jm+g,Km) = 3d0 * blk % P (:,Im,Jm+g-1,Km) - &
                                  3d0 * blk % P (:,Im,Jm+g-2,Km) + &
                                        blk % P (:,Im,Jm+g-3,Km)
        case(5)
          blk % P(:,Im,Jm,Km-g) = 3d0 * blk % P (:,Im,Jm,Km-g+1) - &
                                  3d0 * blk % P (:,Im,Jm,Km-g+2) + &
                                        blk % P (:,Im,Jm,Km-g+3)
        case(6)
          blk % P(:,Im,Jm,Km+g) = 3d0 * blk % P (:,Im,Jm,Km+g-1) - &
                                  3d0 * blk % P (:,Im,Jm,Km+g-2) + &
                                        blk % P (:,Im,Jm,Km+g-3)
      end select
    enddo
      
  end subroutine Ghost_Extrapolate


  subroutine Fill_BC_Ghost_Connection ( Is, Js, Ks, Fs, d11s, d12s, d21s, d22s, blks, Pg )
    implicit none
    integer, intent(in) :: Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(ARES_block_type), intent(in) :: blks
    real(R8), intent(inout)           :: Pg(nprim,6)
    ! Local
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

    select case(Fs)
      case(1:2)
        i1 = Is
        j1 = Js - d11s
        k1 = Ks - d12s
        i2 = Is
        j2 = Js + d11s
        k2 = Ks + d12s
        i3 = Is
        j3 = Js - d21s
        k3 = Ks - d22s
        i4 = Is
        j4 = Js + d21s
        k4 = Ks + d22s
      case(3:4)
        i1 = Is - d11s
        j1 = Js 
        k1 = Ks - d12s
        i2 = Is + d11s
        j2 = Js 
        k2 = Ks + d12s
        i3 = Is - d21s
        j3 = Js
        k3 = Ks - d22s
        i4 = Is + d21s
        j4 = Js 
        k4 = Ks + d22s
      case(5:6)
        i1 = Is - d11s
        j1 = Js - d12s
        k1 = Ks 
        i2 = Is + d11s
        j2 = Js + d12s
        k2 = Ks
        i3 = Is - d21s
        j3 = Js - d22s
        k3 = Ks
        i4 = Is + d21s
        j4 = Js + d22s
        k4 = Ks
    end select

    Pg (:,3) = blks % P (:,i1,j1,k1)
    Pg (:,4) = blks % P (:,i2,j2,k2)
    Pg (:,5) = blks % P (:,i3,j3,k3)
    Pg (:,6) = blks % P (:,i4,j4,k4)

  end subroutine Fill_BC_Ghost_Connection


  !> Pg(:,3:6) of a chimera face: the tangential neighbours of its ghost cell, for the viscous stencil of the face
  !> (BC_Connection_Visc), which takes on the side of the ghost cells the central differences along the face of the row of
  !> the ghost cells, as on the side of the cell those of the row of the cell. Pg(:,3:4) are the lower and the upper
  !> neighbour in the first tangential direction of the face, Pg(:,5:6) in the second (j and k for faces 1-2, i and k for
  !> 3-4, i and j for 5-6). Within the block face a neighbour is the ghost cell of the row, filled from the donors
  !> (Ghost_Chimera). Beyond an edge of the block face the row has no ghost cell: the neighbour is extrapolated linearly
  !> along the row, in the variables whose differences the stencil takes, so that the central difference there is the
  !> one-sided one, exact for a linear field as that of the row of the cell is (the ghost cell copied there halves it, and
  !> the face takes 3/4 of the derivative of a linear field along that direction). A block one cell wide in a direction has
  !> no row to extrapolate along: in the k direction of a 2-D mesh the neighbour is the ghost cell turned by the wedge
  !> angle, as the k neighbours of a cell are (Ghost_Axysimmetry; no turn on a plane mesh), elsewhere the ghost cell
  !> itself. Nothing is written in the block: two chimera faces that meet at an edge of a block share no written location.
  subroutine Fill_BC_Ghost_Chimera ( Im, Jm, Km, Fm, blk, Pg )
    use ARES_Lib_Metrics, only: delthe
    implicit none
    integer, intent(in)               :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(in) :: blk
    real(R8), intent(inout)           :: Pg(nprim,6)
    ! Local
    integer :: c(3), g(3), e(3), t(2), n, s, m

    c = [ Im, Jm, Km ]
    g = c - guide(Fm,:)
    select case (Fm)
      case (1:2)
        t = [ 2, 3 ]
      case (3:4)
        t = [ 1, 3 ]
      case default
        t = [ 1, 2 ]
    end select

    do n = 1, 2
      e = 0
      e(t(n)) = 1
      do s = -1, 1, 2
        m = 2*n + 1 + (s+1)/2   ! Pg(:,3) and Pg(:,4) for n = 1, Pg(:,5) and Pg(:,6) for n = 2
        if ( c(t(n)) + s >= 1 .and. c(t(n)) + s <= blk % dim(t(n)) ) then
          ! Within the block face: the ghost cell of the row
          Pg(:,m) = blk % P (:,g(1)+s*e(1),g(2)+s*e(2),g(3)+s*e(3))
        elseif ( blk % dim(t(n)) > 1 ) then
          ! Beyond an edge of the block face: extrapolated along the row
          Pg(:,m) = Visc_Extrapolated ( blk % P (:,g(1),g(2),g(3)), blk % P (:,g(1)-s*e(1),g(2)-s*e(2),g(3)-s*e(3)) )
        else
          ! One cell in this direction: the ghost cell, in the k direction of a 2-D mesh turned by the wedge angle as by
          ! Ghost_Axysimmetry beyond face 5 (s = -1) or face 6 (s = 1)
          Pg(:,m) = blk % P (:,g(1),g(2),g(3))
          if ( t(n) == 3 .and. ndir < 3 ) then
            Pg(nv,m) = blk % P (nv,g(1),g(2),g(3)) * cos(delthe) - s * blk % P (nw,g(1),g(2),g(3)) * sin(delthe)
            Pg(nw,m) = blk % P (nv,g(1),g(2),g(3)) * sin(delthe) * s + blk % P (nw,g(1),g(2),g(3)) * cos(delthe)
          end if
        end if
      enddo
    enddo

  contains

    !> The state whose viscous variables (velocity, temperature, and the turbulence variables per unit mass: those of
    !> Visc_Variables in BC_Connection_Visc) are twice those of the state a minus those of the state b. The state is
    !> taken at the pressure of a, which is not a viscous variable: the enthalpy of the extrapolated temperature at that
    !> pressure, and the turbulence variables times the density of that state
    function Visc_Extrapolated ( a, b ) result ( q )
      use FLINT_Lib_Thermodynamic
      implicit none
      real(R8), intent(in) :: a(nprim), b(nprim)
      real(R8)             :: q(nprim)
      ! Local
      real(R8) :: rho_a, rho_b, rho_q, temperature

      rho_a = ph2vars( a(np), a(nh), rho_tab )
      rho_b = ph2vars( b(np), b(nh), rho_tab )
      temperature = 2d0 * ph2vars( a(np), a(nh), T_tab ) - ph2vars( b(np), b(nh), T_tab )
      q(np)    = a(np)
      q(nu:nw) = 2d0 * a(nu:nw) - b(nu:nw)
      q(nh)    = pT2h( q(np), temperature )
      if ( nprim > nh ) then
        rho_q = ph2vars( q(np), q(nh), rho_tab )
        q(nt:nprim) = ( 2d0 * a(nt:nprim) / rho_a - b(nt:nprim) / rho_b ) * rho_q
      end if

    end function Visc_Extrapolated

  end subroutine Fill_BC_Ghost_Chimera


  !> Ghost cells of a wall (103, 104, 301, 302) or of an extrapolation face (400), both layers
  !> (docs/theory/numerics.md, "Ghost cells of walls and 400 faces"). The quantities whose differences the viscous fluxes
  !> and the turbulence sources take, T, u and the turbulence variables per unit mass phi = P/rho, are extrapolated, and
  !> the pressure; the ghost state follows from them, h = h(p, T) and rho = rho(p, h) of the fluid tables. A quantity with
  !> a wall value q_w takes 2 q_w - q_1: the velocity of a wall (zero) and the turbulence variables of the walls (their
  !> wall values, as the wall routines take them, so that the ghost values written here are the only ones: the fluxes and
  !> the source terms read the same); the temperature of a 302 wall, max(2 T_w - T_1, T_w/2). The others
  !> q_1 + psi (q_1 - q_2), psi <= s = l0/l1 (rec%sx), cut where needed to keep the pressure, the temperature and each
  !> turbulence variable within [q_1/2, 2 q_1] (Bounded_Extrapolation), the velocity not cut. Exact for the linear fields
  !> whose wall values are their values at the centre of the face, at the positions of ARES_Lib_Stencil_Metrics
  !> (2 x_f - x_1 for the quantities with a wall value, x_1 + s (x_1 - x_2) for the extrapolated ones). In an inviscid run
  !> the walls take no wall values: every quantity is extrapolated, as behind a 400 face. With one cell along the direction
  !> of the face the extrapolated quantities are copied. The pressure and the temperature of a ghost cell are then kept
  !> within the fluid tables (Table_Bounds), where h(p, T) and rho(p, h) read them.
  subroutine Ghost_Wall ( rec, blk )
    use ARES_Config_Types_m, only: obj_rans
    use ARES_Lib_RANS, only: RANS_Set_Wall_Values
    use ARES_Lib_BC_Fluxes, only: Face_Index
    use FLINT_Lib_Thermodynamic, only: ph2vars, pT2h, rho_tab, T_tab, mi_tab2D
    implicit none
    type(ARES_bc_type), intent(in)       :: rec
    type(ARES_block_type), intent(inout) :: blk
    ! Local
    integer  :: Im, Jm, Km, Fm, Ig, Jg, Kg, c(3,2), l, s, Dir, Face_i, Face_j, Face_k
    real(R8) :: p(2), T(2), rho(2), u(3,2), phi(nrans,2)
    real(R8) :: pG(2), TG(2), uG(3,2), phiG(nrans,2), phiw(nrans), rw(nrans)
    real(R8) :: sx, h_w, rho_w, mil, hs, metric, hG
    logical  :: wall, cut, cutq, inside, ok, table

    Im = rec % i ; Jm = rec % j ; Km = rec % k ; Fm = rec % f
    sx = rec % sx
    wall = model > 0 .and. Is_Wall(rec % type)

    ! The boundary cell (1) and the next one along the direction of the face (2; the cell itself in a block of one cell)
    c(:,1) = [Im, Jm, Km]
    c(:,2) = c(:,1)
    if ( blk % dim((Fm+1)/2) > 1 ) c(:,2) = c(:,1) + guide(Fm,:)
    do l = 1, 2
      p(l)     = blk % P(np,c(1,l),c(2,l),c(3,l))
      rho(l)   = ph2vars ( p(l), blk % P(nh,c(1,l),c(2,l),c(3,l)), rho_tab )
      T(l)     = ph2vars ( p(l), blk % P(nh,c(1,l),c(2,l),c(3,l)), T_tab )
      u(:,l)   = blk % P(nu:nw,c(1,l),c(2,l),c(3,l))
      phi(:,l) = blk % P(nt:nprim,c(1,l),c(2,l),c(3,l)) / rho(l)
    enddo

    cut = .false.

    ! Pressure
    call Bounded_Extrapolation ( p(1), p(2), sx, pG(1), pG(2), cutq )
    cut = cut .or. cutq

    ! Temperature
    if ( wall .and. rec % type == 302 ) then
      TG(1) = max ( 2d0 * rec % Tw - T(1), 0.5d0 * rec % Tw )
      TG(2) = max ( 2d0 * rec % Tw - T(2), 0.5d0 * rec % Tw )
      if ( 2d0 * rec % Tw - T(1) < 0.5d0 * rec % Tw ) then
        !$omp atomic
        n_ghost_floor = n_ghost_floor + 1
      end if
    else
      call Bounded_Extrapolation ( T(1), T(2), sx, TG(1), TG(2), cutq )
      cut = cut .or. cutq
    end if

    ! Velocity: zero at the walls
    if ( wall ) then
      uG(:,1) = - u(:,1)
      uG(:,2) = - u(:,2)
    else
      uG(:,1) = u(:,1) + sx * ( u(:,1) - u(:,2) )
      uG(:,2) = u(:,1) + 2d0 * sx * ( u(:,1) - u(:,2) )
    end if

    ! Turbulence variables: the wall values of the wall routines (302: BC_Wall_Temperature, the wall state at T_w; 301:
    ! BC_Wall_Heat, the state of the cell; the coupled walls 103 and 104, whose viscous flux the other solver gives: the
    ! state of the cell, as an adiabatic wall). Every wall, the coupled ones too, takes the roughness of the nearest wall
    ! (blk%ks, with SA-rough only), as the wall flux and the source terms of the cell take it: ks = 0 gives the smooth value
    if ( wall .and. model == 2 ) then
      call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
      metric = dot_product ( blk % M(Im,Jm,Km) % c(Dir,:), blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n )
      hs = 0d0
      if ( obj_rans % rough ) hs = blk % ks(Im,Jm,Km)
      if ( rec % type == 302 ) then
        h_w   = pT2h ( p(1), rec % Tw )
        rho_w = ph2vars ( p(1), h_w, rho_tab )
        mil   = ph2vars ( p(1), h_w, mi_tab2D )
      else
        rho_w = rho(1)
        mil   = ph2vars ( p(1), blk % P(nh,Im,Jm,Km), mi_tab2D )
      end if
      call RANS_Set_Wall_Values ( rans_cell=blk % P(nt:nprim,Im,Jm,Km), mil_wall=mil, rans_wall=rw, metric=metric, &
                                  dist=blk % yn(Im,Jm,Km), ks=hs )
      phiw = rw / rho_w
      phiG(:,1) = 2d0 * phiw - phi(:,1)
      phiG(:,2) = 2d0 * phiw - phi(:,2)
    else
      do s = 1, nrans
        call Bounded_Extrapolation ( phi(s,1), phi(s,2), sx, phiG(s,1), phiG(s,2), cutq )
        cut = cut .or. cutq
      enddo
    end if

    if ( cut ) then
      !$omp atomic
      n_ghost_cut = n_ghost_cut + 1
    end if

    ! Ghost state of the two layers
    table = .false.
    do l = 1, 2
      call Table_Bounds ( pG(l), TG(l), inside, ok )
      table = table .or. .not. inside
      if ( ok ) then
        hG = pT2h ( pG(l), TG(l) )
      else
        pG(l) = p(1)                                   ! no temperature is in both rows of the tables: the cell state
        hG = blk % P(nh,Im,Jm,Km)
      end if
      Ig = Im - l*guide(Fm,1)
      Jg = Jm - l*guide(Fm,2)
      Kg = Km - l*guide(Fm,3)
      blk % P(np,Ig,Jg,Kg)    = pG(l)
      blk % P(nu:nw,Ig,Jg,Kg) = uG(:,l)
      blk % P(nh,Ig,Jg,Kg)    = hG
      if ( nrans > 0 ) blk % P(nt:nprim,Ig,Jg,Kg) = phiG(:,l) * ph2vars ( pG(l), hG, rho_tab )
    enddo

    if ( table ) then
      !$omp atomic
      n_ghost_table = n_ghost_table + 1
    end if

  end subroutine Ghost_Wall


  !> The pressure and the temperature of a ghost state brought within the fluid tables, where h(p, T) (pT2h) and rho(p, h)
  !> (ph2vars) read them without leaving their arrays: p within the pressures of the tables, T within the temperatures of
  !> both rows of the tables that bracket p, from h_min to h_max. inside: they were already within. ok: the two rows have
  !> temperatures in common (false only on tables whose T range jumps from one pressure to the next, where the caller
  !> takes the state of the cell).
  subroutine Table_Bounds ( p, T, inside, ok )
    use FLINT_Lib_Thermodynamic, only: pmin, deltap, Tmin2, deltaT, h_tab2D
    implicit none
    real(R8), intent(inout) :: p, T
    logical,  intent(out)   :: inside, ok
    ! Local
    integer  :: i, ip, ih                ! ip, ih: the last pressure and enthalpy indices of the tables
    real(R8) :: q, lo, hi

    ip = ubound ( h_tab2D, 1 )
    ih = ubound ( h_tab2D, 2 )
    q  = min ( max ( p, pmin ), pmin + ( ip - 1d-6 ) * deltap )
    i  = int ( ( q - pmin ) / deltap )
    lo = max ( Tmin2(i), Tmin2(i+1) )
    hi = min ( Tmin2(i) + ( ih - 1d-6 ) * deltaT(i), Tmin2(i+1) + ( ih - 1d-6 ) * deltaT(i+1) )
    ok = lo <= hi
    inside = ok .and. q == p .and. T >= lo .and. T <= hi
    p = q
    T = min ( max ( T, lo ), hi )

  end subroutine Table_Bounds


  !> q_1 + psi (q_1 - q_2) on the first ghost layer and q_1 + 2 psi_2 (q_1 - q_2) on the second, with psi <= sx and
  !> psi_2 <= psi the largest values that keep the ghost value within [q_1/2, 2 q_1] (with the sign of q_1; a zero q_1 is
  !> copied). cut: the bound acted on the first layer.
  pure subroutine Bounded_Extrapolation ( q1, q2, sx, g1, g2, cut )
    implicit none
    real(R8), intent(in)  :: q1, q2, sx
    real(R8), intent(out) :: g1, g2
    logical,  intent(out) :: cut
    ! Local
    real(R8) :: d, b

    d = q1 - q2
    b = huge(1d0)
    if ( q1 * d > 0d0 ) then
      b = q1 / d                    ! away from zero: |g| <= 2 |q1|
    elseif ( q1 * d < 0d0 ) then
      b = -0.5d0 * q1 / d           ! towards zero: |g| >= |q1|/2
    elseif ( abs(d) > 0d0 ) then
      b = 0d0                       ! q1 = 0
    end if
    g1 = q1 + min ( sx, b ) * d
    g2 = q1 + 2d0 * min ( sx, 0.5d0 * b ) * d
    cut = b < sx

  end subroutine Bounded_Extrapolation


  !> One line at the end of a run in which the cut of an extrapolated quantity acted in a ghost fill of a wall or 400 face,
  !> the temperature of a 302 wall took the floor T_w/2, or the pressure or the temperature of a ghost cell was brought
  !> back within the fluid tables: the fills of each kind, over the whole run, its grid levels and ranks. Each says that the
  !> mesh does not resolve the profile next to the face (the last one, also that the tables do not reach the ghost state);
  !> a run in which none acted prints no line. Collective: every rank calls it, and the three reductions, in every run, so
  !> that no rank takes another MPI path; rank 0 writes the line when one of the totals is not zero.
  subroutine Report_Wall_Ghosts ()
    use ARES_Mod_MPI, only: mpi_is_root, mpi_allreduce_sum_r8
    implicit none
    real(R8) :: ncut, nfloor, ntable

    call mpi_allreduce_sum_r8 ( real(n_ghost_cut, R8), ncut )
    call mpi_allreduce_sum_r8 ( real(n_ghost_floor, R8), nfloor )
    call mpi_allreduce_sum_r8 ( real(n_ghost_table, R8), ntable )
    if ( mpi_is_root .and. ( ncut > 0d0 .or. nfloor > 0d0 .or. ntable > 0d0 ) ) &
      write(*,'(A,I0,A,I0,A,I0)') ' Ghost cells of walls and 400 faces, fills with the extrapolation cut: ', &
                                  nint(ncut, I8), ', with the 302 temperature at T_w/2: ', nint(nfloor, I8), &
                                  ', brought within the fluid tables: ', nint(ntable, I8)

  end subroutine Report_Wall_Ghosts


  function Is_Wall(type) result (ans)
      implicit none
      integer, intent(in) :: type
      logical :: ans

      select case (type)
        case (103, 104, 301, 302)  ! wall BCs (heat flux, coupled)
          ans = .true.
        case default
          ans = .false.
      end select
    end function Is_Wall

end module ARES_Lib_Ghost