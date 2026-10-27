!>@brief: Module to evaluate correction terms (depending on DSij/Dt) in -RC turbulence corrections:
!> the factor 2 r*/(1 + r*) of f_r1 (r* = S/Omega) and the numerator of r~, 2 W_ik S_jk DS_ij/Dt, shared by
!> SA-RC (Spalart & Shur 1997; Shur et al. 2000) and SST-RC (Smirnov & Menter 2009). DS_ij/Dt is the
!> steady convective part u.grad(S_ij).
module ARES_Lib_SpalartShur
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains
  
  subroutine compute_velocity_gradient ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Mod_MPI, only: is_local_block

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, n
    real(kind=8), dimension(3,3) :: vel_diff

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3) private(vel_diff)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        ! Cell-centered velocity gradient
        vel_diff(:,1) = ( domain%blk(b)%P(nu:nw,i+1,j,k) - domain%blk(b)%P(nu:nw,i-1,j,k) )/2d0
        vel_diff(:,2) = ( domain%blk(b)%P(nu:nw,i,j+1,k) - domain%blk(b)%P(nu:nw,i,j-1,k) )/2d0
        vel_diff(:,3) = ( domain%blk(b)%P(nu:nw,i,j,k+1) - domain%blk(b)%P(nu:nw,i,j,k-1) )/2d0
        domain%blk(b)%vel_gradient(i,j,k)%c = matmul ( vel_diff + 1d-40, domain%blk(b)%m(i,j,k)%c )
      end do ; end do ; end do

      ! Extrapolation in ghost cells -as second derivatives must be computed
      ! i faces
      !$omp do collapse(2) private(n)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
        n = domain % blk(b) % dim(1) + 1
        domain % blk(b) % vel_gradient(0,j,k) % c = domain % blk(b) % vel_gradient(1,j,k) % c
        domain % blk(b) % vel_gradient(n,j,k) % c = domain % blk(b) % vel_gradient(n-1,j,k) % c
      end do ; end do

      ! j faces
      !$omp do collapse(2) private(n)
      do k = 1, domain % blk(b) % dim(3)
      do i = 1, domain % blk(b) % dim(1)   
        n = domain % blk(b) % dim(2) + 1
        domain % blk(b) % vel_gradient(i,0,k) % c = domain % blk(b) % vel_gradient(i,1,k) % c
        domain % blk(b) % vel_gradient(i,n,k) % c = domain % blk(b) % vel_gradient(i,n-1,k) % c
      end do ; end do

      ! k faces
      !$omp do collapse(2) private(n)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)       
        n = domain % blk(b) % dim(3) + 1
        domain % blk(b) % vel_gradient(i,j,0) % c = domain % blk(b) % vel_gradient(i,j,1) % c
        domain % blk(b) % vel_gradient(i,j,n) % c = domain % blk(b) % vel_gradient(i,j,n-1) % c
      end do ; end do

    end do

  end subroutine compute_velocity_gradient


  ! Computation of 2 r*/(1 + r*) and rtilde (missing denominator) for the Spalart-Shur correction.
  ! fr1 (the overall correction term) is not computed here as the model constants and rtilde
  ! denominator definition are different between SA and SST.
  subroutine compute_rc_terms ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Lib_Fluid, only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    !use ARES_Lib_RotatingFrame, only: obj_rot
    use ARES_Mod_MPI, only: is_local_block

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, ii, jj, kk
    real(kind=8), dimension(3,3,3) :: strain_gradient
    real(kind=8), dimension(3,3) :: Sij, Wij, strain_transport
    real(kind=8) :: omega(3), S, O, rtilde

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3) private (i, j, k, ii, jj, kk, strain_gradient), &
      !$omp private (Sij, Wij, strain_transport, omega, S, O, rtilde)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        ! Finite difference strain rate tensor
        strain_gradient(:,:,1) = ( strain_tensor(domain%blk(b)%vel_gradient(i+1,j,k)%c) & 
                                 - strain_tensor(domain%blk(b)%vel_gradient(i-1,j,k)%c) )*0.5d0
        strain_gradient(:,:,2) = ( strain_tensor(domain%blk(b)%vel_gradient(i,j+1,k)%c) & 
                                 - strain_tensor(domain%blk(b)%vel_gradient(i,j-1,k)%c) )*0.5d0
        strain_gradient(:,:,3) = ( strain_tensor(domain%blk(b)%vel_gradient(i,j,k+1)%c) & 
                                 - strain_tensor(domain%blk(b)%vel_gradient(i,j,k-1)%c) )*0.5d0

        ! Transformation component by component
        do ii = 1, 3
        do jj = 1, 3
          strain_gradient(ii,jj,:) = matmul ( strain_gradient(ii,jj,:), domain%blk(b)%m(i,j,k)%c )
        end do
        end do

        Sij = strain_tensor( domain%blk(b)%vel_gradient(i,j,k)%c )
        omega = vorticity_vector(domain%blk(b)%vel_gradient(i,j,k)%c)
        Wij = vorticity_tensor( omega )
        S = sqrt ( 2d0*sum (Sij**2) )
        O = sqrt ( sum(omega**2) )
        ! The factor 2 r*/(1 + r*) of fr1, r* = S/O, as 2 S/(S + O): where the vorticity vanishes it takes its limit 2.
        ! Where strain and vorticity are both below 1e-30 1/s, in a fluid at rest or a uniform flow, r* is not defined
        ! (the velocity gradient there is the 1d-40 added to the differences in compute_velocity_gradient, times the
        ! metric of the cell): the factor is 1, the value of simple shear, and rtilde is 0 below, so the correction
        ! does not act and the production is that of the model without it
        if ( S + O > 1d-30 ) then
          domain%blk(b)%rc_term1(i,j,k) = 2d0*S/(S + O)
        else
          domain%blk(b)%rc_term1(i,j,k) = 1d0
        end if
        ! Material derivative of strain rate tensor, neglecting the time term
        do ii = 1, 3
        do jj = 1, 3
          strain_transport(ii,jj) = dot_product ( domain%blk(b)%p(nu:nw,i,j,k), &
                                                  strain_gradient(ii,jj,:) )
        end do
        end do
        rtilde = 0d0
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
          rtilde = rtilde + Wij(ii,kk)*Sij(jj,kk)*strain_transport(ii,jj)
        enddo
        enddo
        enddo
        ! Note ! Denominator is missing as its different between models
        ! Below the threshold of the factor, rtilde is a ratio of the 1d-40 of the velocity gradient, not zero where the
        ! metric of the cells changes along the flow: it is 0 there, as in simple shear
        if ( S + O <= 1d-30 ) rtilde = 0d0
        domain%blk(b)%rc_term2(i,j,k) = 2d0*rtilde 
        ! Note !
      end do
      end do
      end do
    end do

  end subroutine compute_rc_terms

end module ARES_Lib_SpalartShur