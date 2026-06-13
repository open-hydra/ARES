module ARES_Lib_Ghost
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Parameters_m
  use ARES_Global_m

  implicit none

contains

  subroutine Fill_Ghost_Cell ( domain )
    use ARES_Mod_MPI, only: is_local_block, mpi_size_
    use ARES_Mod_GhostExchange, only: exchange_ghost_P_post_recv, exchange_ghost_P_pack, &
                                       exchange_ghost_P_post_send, exchange_ghost_P_wait_unpack, &
                                       exchange_ghost_P_wait_send, &
                                       Ghost_Interrank, exchange_ghost_Pg, ghost_sched
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: ii, i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    integer :: fg


    ! MPI: post persistent receives, pack buffer in parallel, then post sends
    !$omp single
    call exchange_ghost_P_post_recv(domain)
    !$omp end single

    ! Pack send buffer in parallel over face groups
    !$omp do schedule(static) private(fg)
    do fg = 1, ghost_sched%n_send_faces
      call exchange_ghost_P_pack(domain, fg, fg)
    end do

    ! Post sends (must wait for all packing to complete — implicit barrier from !$omp do)
    !$omp single
    call exchange_ghost_P_post_send(domain)
    !$omp end single nowait

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
        case default
          call Ghost_Extrapolate ( Im, Jm, Km, Fm, domain % blk(Bm) )
      end select
    enddo

    ! Compute Pg(:,3:6) for type-1 connections where source block is local.
    ! This only reads local blk%P data, so it can run before the MPI exchange.
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

    ! Compute Pg
    !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm)
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
    enddo

    ! MPI: wait for P receives to complete
    !$omp single
    call exchange_ghost_P_wait_unpack(domain)
    !$omp end single

    ! Process INTER-RANK type-1 entries (Bm local, Bs remote)
    !$omp do schedule (dynamic) private(ii, i, Bm, Im, Jm, Km, Fm, Bs)
    do ii = 1, domain % n_local_bc
      i  = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 101) cycle
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
    !$omp single
    call exchange_ghost_P_wait_send(domain)
    call exchange_ghost_Pg(domain)
    !$omp end single

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
    use ARES_Lib_Metrics, only: delthe
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
      blk % P (nv,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Km) * cos(delthe) + &
                guide(Fm,3) * blk % P (nw,Im,Jm,Km) * sin(delthe)
      blk % P (nw,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Km) * sin(delthe) * (-guide(Fm,3)) &
                            + blk % P (nw,Im,Jm,Km) * cos(delthe)
      blk % P (1:nu,Ig,Jg,Kg) = blk % P (1:nu,Im,Jm,Km)
      blk % P (nh:nprim,Ig,Jg,Kg) = blk % P (nh:nprim,Im,Jm,Km)
    else
      ! cell is connected with opposite K-face, with v/w rotated of +-delthe
      Ks = blk % dim(3) - Km + 1
      blk % P (nv,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Ks) * cos(delthe) + &
                    guide(Fm,3) * blk % P (nw,Im,Jm,Ks) * sin(delthe)
      blk % P (nw,Ig,Jg,Kg) = blk % P (nv,Im,Jm,Ks) * sin(delthe) + &
                 (-guide(Fm,3)) * blk % P (nw,Im,Jm,Ks) * cos(delthe)
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


  subroutine Ghost_Wall_Extrapolation ( domain )
    implicit none
    type(ARES_domain_type) :: domain
    integer :: ii, i, ig, jg, kg
    integer :: bm, im, jm, km, fm

    !$omp do private ( ii, i, bm, im, jm, km, fm, ig, jg, kg )
    ! BC array processing — iterate only over local BC entries
    do ii = 1, domain%n_local_bc
      i = domain%local_bc_idx(ii)

      ! Preliminary assignments
      bm = domain % bc(i) % b
      im = domain % bc(i) % i
      jm = domain % bc(i) % j
      km = domain % bc(i) % k
      fm = domain % bc(i) % f

      ! ghost cell coordinates
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)

      if ( Is_Wall( domain%bc(i)%type ) ) then
        ! Invert the velocity components
        domain%blk(bm)%P (nu:nw,ig,jg,kg) = - domain%blk(bm)%P (nu:nw,im,jm,km)
        ! NOTE: the extrapolation of the turbulence variables is done in BC routines, due
        ! to the fact that omega_wall needs to be set to extrapolate its value in ghost cell

        ! Other variables copied
        domain%blk(bm)%P (np,ig,jg,kg) = domain%blk(bm)%P (np,im,jm,km)
        domain%blk(bm)%P (nh:nprim,ig,jg,kg) = domain%blk(bm)%P (nh:nprim,im,jm,km)
      endif
    enddo

  end subroutine Ghost_Wall_Extrapolation


  function Is_Wall(type) result (ans)
      implicit none
      integer, intent(in) :: type
      logical :: ans

      select case (type)
        case (103, 301, 302)  ! wall BCs (heat flux, coupled)
          ans = .true.
        case default
          ans = .false.
      end select
    end function Is_Wall

end module ARES_Lib_Ghost