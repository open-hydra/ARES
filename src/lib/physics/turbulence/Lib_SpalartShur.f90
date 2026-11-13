!>@brief: Module to evaluate correction terms (depending on DSij/Dt) in -RC turbulence corrections:
!> the factor 2 r*/(1 + r*) of f_r1 (r* = S/Omega) and the numerator of r~, 2 W_ik S_jk DS_ij/Dt, shared by
!> SA-RC (Spalart & Shur 1997; Shur et al. 2000) and SST-RC (Smirnov & Menter 2009). DS_ij/Dt is the
!> steady convective part u.grad(S_ij).
module ARES_Lib_SpalartShur
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains
  
  !> Velocity gradient of every cell, from the central differences of the velocity along the index directions and the
  !> cell metric of the source terms (Source_Metrics), and in the first layer of ghost cells the gradient that the
  !> strain-rate gradient of compute_rc_terms takes there: behind the faces between blocks and of rotational periodicity
  !> the gradient of the cell across (Ghost_Velocity_Gradient), behind the other faces a copy of the cell.
  subroutine compute_velocity_gradient ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Mod_MPI, only: is_local_block
    use ARES_Base_Types_m, only: ARES_tensor_3D_type
    use ARES_Lib_Stencil_Metrics, only: Source_Metrics

    implicit none
    type(ARES_domain_type), intent(inout), target :: domain
    ! Local
    integer :: b, i, j, k, n
    real(kind=8), dimension(3,3) :: vel_diff
    type(ARES_tensor_3D_type), pointer :: Mi(:,:,:)   ! cell metric of the sources (Source_Metrics)

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Source_Metrics ( domain % blk(b), Mi )
      !$omp do collapse(3) private(vel_diff)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        ! Cell-centered velocity gradient
        vel_diff(:,1) = ( domain%blk(b)%P(nu:nw,i+1,j,k) - domain%blk(b)%P(nu:nw,i-1,j,k) )/2d0
        vel_diff(:,2) = ( domain%blk(b)%P(nu:nw,i,j+1,k) - domain%blk(b)%P(nu:nw,i,j-1,k) )/2d0
        vel_diff(:,3) = ( domain%blk(b)%P(nu:nw,i,j,k+1) - domain%blk(b)%P(nu:nw,i,j,k-1) )/2d0
        domain%blk(b)%vel_gradient(i,j,k)%c = matmul ( vel_diff + 1d-40, Mi(i,j,k)%c )
      end do ; end do ; end do

      ! Ghost cells, as second derivatives must be computed: a copy of the cell (replaced below behind the faces between
      ! blocks and of rotational periodicity)
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

    call Ghost_Velocity_Gradient ( domain )

  end subroutine compute_velocity_gradient


  !> First layer of ghost cells of the velocity gradient behind the faces between blocks (101), of a periodic pair (201)
  !> and of rotational periodicity (200): the gradient of the cell across, as Ghost_Connection and Ghost_Axysimmetry fill
  !> the velocity, so that next to these faces the strain-rate gradient of compute_rc_terms is that of the undivided mesh.
  !> The gradient is a tensor in Cartesian components: it is copied whatever the orientation of the two blocks, and from
  !> another rank by exchange_velocity_gradient. Behind a 200 face it is that of the image cell turned about x as the
  !> velocity is, R G R^T; with a zero sector angle (a planar mesh, a translated sector) the copy of the cell stays.
  !> Cell_Source_Metrics puts each of these values at the position of its cell.
  subroutine Ghost_Velocity_Gradient ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: ndir
    use ARES_Parameters_m, only: guide
    use ARES_Lib_Metrics, only: delthe, cos_delthe, sin_delthe
    use ARES_Mod_MPI, only: is_local_block, mpi_size_
    use ARES_Mod_GhostExchange, only: exchange_velocity_gradient

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer  :: ii, r, b, bs, f, g(3), s(3)
    real(R8) :: Rot(3,3), cs, sn
    logical  :: turned

    turned = .not. ( delthe >= 0d0 .and. delthe <= 0d0 )   ! a wedge or a sector turned about x
    cs = cos_delthe                                         ! computed once with delthe (Check_Mesh_Type)
    sn = sin_delthe

    !$omp do schedule (dynamic) private(ii, r, b, bs, f, g, s, Rot)
    do ii = 1, domain % n_local_bc
      r = domain % local_bc_idx(ii)
      b = domain % bc(r) % b
      f = domain % bc(r) % f
      g = [ domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k ] - guide(f,:)
      select case ( domain % bc(r) % type )
        case (101, 201)
          bs = domain % bc(r) % bs
          if ( .not. is_local_block(bs) ) cycle              ! from another rank, below
          domain % blk(b) % vel_gradient(g(1),g(2),g(3)) % c = &
            domain % blk(bs) % vel_gradient(domain % bc(r) % is, domain % bc(r) % js, domain % bc(r) % ks) % c
        case (200)
          if ( .not. turned ) cycle
          ! the image of Ghost_Axysimmetry: the cell itself on a wedge, the cell next to the opposite k face of a sector,
          ! turned about x by -guide(f,3) delthe (v' = cos v + guide(f,3) sin w, w' = -guide(f,3) sin v + cos w)
          s = [ domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k ]
          if ( ndir == 3 ) s(3) = domain % blk(b) % dim(3) - s(3) + 1
          Rot = reshape ( [ 1d0, 0d0, 0d0,  0d0, cs, -guide(f,3)*sn,  0d0, guide(f,3)*sn, cs ], [3,3] )
          domain % blk(b) % vel_gradient(g(1),g(2),g(3)) % c = &
            matmul ( Rot, matmul ( domain % blk(b) % vel_gradient(s(1),s(2),s(3)) % c, transpose ( Rot ) ) )
      end select
    enddo
    !$omp end do

    ! Faces between ranks: the gradient of the source cell, by message (MPI_THREAD_FUNNELED: the primary thread)
    if ( mpi_size_ > 1 ) then
      !$omp master
      call exchange_velocity_gradient ( domain )
      !$omp end master
      !$omp barrier
    end if

  end subroutine Ghost_Velocity_Gradient


  ! Computation of 2 r*/(1 + r*) and rtilde (missing denominator) for the Spalart-Shur correction.
  ! fr1 (the overall correction term) is not computed here as the model constants and rtilde
  ! denominator definition are different between SA and SST.
  subroutine compute_rc_terms ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Lib_Fluid, only: Strain_Tensor, Vorticity_Vector, Vorticity_Tensor
    !use ARES_Lib_RotatingFrame, only: obj_rot
    use ARES_Mod_MPI, only: is_local_block
    use ARES_Base_Types_m, only: ARES_tensor_3D_type
    use ARES_Lib_Stencil_Metrics, only: Source_Metrics

    implicit none
    type(ARES_domain_type), intent(inout), target :: domain
    ! Local
    integer :: b, i, j, k, ii, jj, kk
    real(kind=8), dimension(3,3,3) :: strain_gradient
    real(kind=8), dimension(3,3) :: Sij, Wij, strain_transport
    real(kind=8) :: omega(3), S, O, rtilde
    type(ARES_tensor_3D_type), pointer :: Mi(:,:,:), Ma(:,:,:)   ! Ma: the positions of the ghost velocity gradients

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Source_Metrics ( domain % blk(b), Mi, Ma=Ma )
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
          strain_gradient(ii,jj,:) = matmul ( strain_gradient(ii,jj,:), Ma(i,j,k)%c )
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