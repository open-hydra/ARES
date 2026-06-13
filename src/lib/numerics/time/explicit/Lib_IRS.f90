!>@brief: Implicit Residual Smoothing related subroutines.
module ARES_Lib_IRS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Residual_Smoothing

  integer, parameter, private :: njacobi = 2

contains

  !% Compute the Implicit Residual Smooting operator: R* = R + beta*LAPLACE(R*) with Jacobi iterations.
  subroutine Residual_Smoothing ( domain, beta )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Lib_Ghost
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: beta
    ! Local
    integer :: b, d

    do d = 1, ndir

      ! Residual extrapolation in the ghost cells.
      call Ghost_Residual_Extrapolation ( domain )

      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle

        call Residual_Smoothing_Blk ( domain % blk(b) % r, domain % blk(b) % rs1, &
                                      domain % blk(b) % rs2, domain % blk(b) % dim, d, beta )

      enddo ! blocks
    
    enddo ! direction
  
  end subroutine Residual_Smoothing


  subroutine Residual_Smoothing_Blk ( Res, Res_Star, Res_Star_New, n, Dir, beta )
    use ARES_Global_m
    implicit none
    integer, intent(in) :: n(3), Dir
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res, Res_Star, Res_Star_New
    real(R8), intent(in) :: beta
    ! Local
    integer :: i, j, k, jacobi
    integer :: i1, j1, k1, i2, j2, k2
    real(R8), dimension(nprim) :: r, rm, rp

    ! Res* initialized as Res in the block, including ghost cells.
    !$omp do collapse (3)
    do k = 0, n(3) + 1
    do j = 0, n(2) + 1
    do i = 0, n(1) + 1

      Res_Star(:,i,j,k) = Res(:,i,j,k)
      
    enddo; enddo; enddo

    ! Smoothing with jacobi iterations.
    do jacobi = 1, njacobi

      !$omp do collapse (3) private (r, rm, rp, i1, i2, j1, j2, k1, k2)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)

        call Jacobi_Stencil ( i1, j1, k1, i2, j2, k2, i, j, k, Dir )
        r  = Res (:,i,j,k)
        rm = Res_Star (:,i1,j1,k1)
        rp = Res_Star (:,i2,j2,k2)
            
        ! Jacobi: r*_New(i) = [ r(i) + e*( r*(i-1) + r*(i+1) ) ] / [ 1 + 2e ]
        Res_Star_New(:,i,j,k) = ( r + beta * ( rm + rp ) ) / ( 1d0 + 2d0*beta )

      enddo; enddo; enddo

      !$omp do collapse (3)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)

        Res_Star(:,i,j,k) = Res_Star_New(:,i,j,k)

      enddo; enddo; enddo

    enddo ! jacobi

    ! R* after i-direction smoothing becomes the new R.
    !$omp do collapse(3)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      Res(:,i,j,k) = Res_Star(:,i,j,k)

    enddo; enddo; enddo

    contains
      
      subroutine Jacobi_Stencil ( i1, j1, k1, i2, j2, k2, i, j, k, d )
        implicit none
        integer, intent(in)   :: i, j, k, d
        integer, intent(out)  :: i1, j1, k1, i2, j2, k2

        select case (d)
          case (1)
            i1 = i-1
            j1 = j
            k1 = k
            i2 = i+1
            j2 = j
            k2 = k
          case (2)
            i1 = i
            j1 = j-1
            k1 = k
            i2 = i
            j2 = j+1
            k2 = k
          case (3)
            i1 = i
            j1 = j
            k1 = k-1
            i2 = i
            j2 = j
            k2 = k+1
        end select
      
      end subroutine Jacobi_Stencil

  end subroutine Residual_Smoothing_Blk

  
  subroutine Ghost_Residual_Extrapolation ( domain )
    use ARES_Advanced_Types_m
    use ARES_Parameters_m, only: guide
    use ARES_Mod_MPI, only: is_local_block
    use ARES_Mod_GhostExchange, only: exchange_ghost_R
    implicit none
    type(ARES_domain_type) :: domain
    ! Local
    integer :: ii, i, ig, jg, kg
    integer :: bm, im, jm, km, fm, bs, is, js, ks

    ! MPI: exchange residual data for inter-rank type-1 connections
    !$omp single
    call exchange_ghost_R(domain)
    !$omp end single

    ! BC array processing — iterate only over local BC entries
    !$omp do private ( ii, i, bm, im, jm, km, fm, bs, is, js, ks, ig, jg, kg )
    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)

      ! Preliminary assignments.
      bm = domain % bc(i) % b
      im = domain % bc(i) % i
      jm = domain % bc(i) % j
      km = domain % bc(i) % k
      fm = domain % bc(i) % f

      ! Ghost cell coordinates.
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)

      select case( domain % bc(i) % type)

        case (101) ! connection
          bs = domain % bc(i) % bs
          if (.not. is_local_block(bs)) cycle  ! inter-rank already handled by exchange_ghost_R
          is = domain % bc(i) % is
          js = domain % bc(i) % js
          ks = domain % bc(i) % ks

          domain % blk(bm) % r (:,ig,jg,kg) = domain % blk(bs) % r (:,is,js,ks)

        case default  ! wall, inlet, outlet — zero-gradient ghost
          domain % blk(bm) % r (:,ig,jg,kg) = domain % blk(bm) % r (:,im,jm,km)

      end select
    enddo

  end subroutine Ghost_Residual_Extrapolation

end module ARES_Lib_IRS