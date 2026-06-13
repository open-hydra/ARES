!> @brief Module for SSGLRR-RSM-w2019 model. see https://turbmodels.larc.nasa.gov/rsm-ssglrr.html.
module ARES_Lib_SSGLRR
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  real(kind=8), dimension(3,3), parameter :: &
  deltaij = reshape((/ 1d0, 0d0, 0d0, &
                       0d0, 1d0, 0d0, &
                       0d0, 0d0, 1d0  /), shape(deltaij), order=(/2,1/) )
  real(kind=8), parameter, public :: &
  C_mu = 0.09d0, &
  ! omega coefficients
  alphaw_w = 0.5556d0, &
  betaw_w = 0.075d0, &
  sigmaw_w = 0.5d0, &
  sigmad_w = 0d0, &
  C1_w = 1.8d0, &
  C1star_w = 0d0, &
  C2_w = 0d0, &
  C3_w = 0.8d0, &
  C3star_w = 0d0, &
  C2_LRR = 0.52d0, &
  C4_w = 0.5d0*(18d0*C2_LRR + 12d0)/11d0, &
  C5_w = 0.5d0*(-14d0*C2_LRR+20d0)/11d0, &
  D_w = 0.75d0*C_mu, &
  ! epsilon coefficients
  alphaw_e = 0.44d0, &
  betaw_e = 0.0828d0, &
  sigmaw_e = 0.856d0, &
  sigmad_e = 1.712d0, &
  C1_e = 1.7d0, &
  C1star_e = 0.9d0, &
  C2_e = 1.05d0, &
  C3_e = 0.8d0, &
  C3star_e = 0.65d0, &
  C4_e = 0.625d0, &
  C5_e = 0.2d0, &
  D_e = 0.22d0

contains


  subroutine SSGLRR_Source_Terms ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Mod_MPI, only: is_local_block

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call SSGLRR_Blk ( domain % blk(b) % P,    &
                        domain % blk(b) % r,    &
                        domain % blk(b) % M,    &
                        domain % blk(b) % vol,  &
                        domain % blk(b) % yn,   &
                        domain % blk(b) % dim   )
    
    end do

  end subroutine SSGLRR_Source_Terms


  subroutine SSGLRR_Blk ( Prim, Res, M, Volume, walldist, n )
    use ARES_Base_Types_m
    use ARES_Global_m
    use ARES_Lib_Fluid
    use FLINT_Lib_Thermodynamic
    !use ARES_Lib_RotatingFrame, only: obj_rot

    implicit none
    integer, intent(in) :: n(3)
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: walldist
    type(ARES_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    ! Local
    integer :: i, j, k, ii, jj, kk
    real(kind=8) :: rho, mil, rhoRij(3,3), omega, kappa, rhoe, dist
    real(kind=8) :: Gradvel(3,3), Grad(8,3), Divel, Sij(3,3), Sij_star(3,3), Wij(3,3)
    real(kind=8) :: F1, F2, alphaw, betaw, sigmaw, sigmad, C1, C2, C3, C4, C5, C1star, C3star
    real(kind=8), dimension(3,3) :: rhoP, rhoeij, aij, C1_term, C2_term, C3_term, C4_term, C5_term, aikakj, rhoPI 
    real(kind=8) :: rhoPkk, aklakl, aklSkl, Diff, Source(7) 

    !$omp do collapse (3) private (rho, mil, rhoRij, omega, kappa, rhoe, dist) &
    !$omp private (Gradvel, Grad, Divel, Sij, Sij_star, Wij) &
    !$omp private (F1, F2, alphaw, betaw, sigmaw, sigmad, C1, C2, C3, C4, C5, C1star, C3star) &
    !$omp private (rhoP, rhoeij, aij, C1_term, C2_term, C3_term, C4_term, C5_term, aikakj, rhoPI) &
    !$omp private (rhoPkk, aklakl, aklSkl, Diff, Source) &
    !$omp private (i,j,k, ii, jj, kk)

    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      ! Local variables
      rho = ph2vars( Prim(np,i,j,k), Prim(nh,i,j,k),  rho_tab )
      mil = ph2vars( Prim(np,i,j,k), Prim(nh,i,j,k),  mi_tab2D )

      rhoRij(1,1) = Prim(nt,i,j,k)
      rhoRij(2,2) = Prim(nt+1,i,j,k)
      rhoRij(3,3) = Prim(nt+2,i,j,k)
      rhoRij(1,2) = Prim(nt+3,i,j,k)
      rhoRij(1,3) = Prim(nt+4,i,j,k)
      rhoRij(2,3) = Prim(nt+5,i,j,k)
      rhoRij(2,1) = rhoRij(1,2)
      rhoRij(3,1) = rhoRij(1,3)
      rhoRij(3,2) = rhoRij(2,3)
      omega = Prim(nt+6,i,j,k) / rho
      kappa = 0.5d0 * ( rhoRij(1,1) + rhoRij(2,2) + rhoRij(3,3) ) / rho
      rhoe = rho * C_mu * kappa * omega
      dist = walldist(i,j,k)

      ! Velocity gradient computation
      Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
      Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j,k) % c )

      ! Gradient of Rij, omega (1-7) and kappa (8)
      Grad(1:7,1) = ( Prim(nt:nt+6,i+1,j,k) / ph2vars( Prim(np,i+1,j,k), Prim(nh,i+1,j,k),  rho_tab ) - &
                      Prim(nt:nt+6,i-1,j,k) / ph2vars( Prim(np,i-1,j,k), Prim(nh,i-1,j,k),  rho_tab ) ) * 0.5d0
      Grad(1:7,2) = ( Prim(nt:nt+6,i,j+1,k) / ph2vars( Prim(np,i,j+1,k), Prim(nh,i,j+1,k),  rho_tab ) - &
                      Prim(nt:nt+6,i,j-1,k) / ph2vars( Prim(np,i,j-1,k), Prim(nh,i,j-1,k),  rho_tab ) ) * 0.5d0
      Grad(1:7,3) = ( Prim(nt:nt+6,i,j,k+1) / ph2vars( Prim(np,i,j,k+1), Prim(nh,i,j,k+1),  rho_tab ) - &
                      Prim(nt:nt+6,i,j,k-1) / ph2vars( Prim(np,i,j,k-1), Prim(nh,i,j,k-1),  rho_tab ) ) * 0.5d0
      Grad(8,1) = 0.5d0 * Sum ( Grad(1:3,1) )
      Grad(8,2) = 0.5d0 * Sum ( Grad(1:3,2) )
      Grad(8,3) = 0.5d0 * Sum ( Grad(1:3,3) )
      Grad = Matmul ( Grad, M(i,j,k) % c )

      ! Strain tensor Sij, Sij*, and vorticity tensor Wij 
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3) 
      Sij = Strain_Tensor ( Gradvel )
      Sij_star = Sij - 1d0/3d0 * Divel * deltaij
      Wij = Vorticity_Tensor ( Vorticity_Vector ( Gradvel ) )

      ! Blending of coefficients
      F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mil, gradw=Grad(7,:), gradk=Grad(8,:), d=dist )
      F2 = 1d0 - F1
      alphaw = F1*alphaw_w + F2*alphaw_e 
      betaw = F1*betaw_w + F2*betaw_e 
      sigmaw = F1*sigmaw_w + F2*sigmaw_e
      sigmad = F1*sigmad_w + F2*sigmad_e
      C1 = F1*C1_w + F2*C1_e 
      C2 = F1*C2_w + F2*C2_e 
      C3 = F1*C3_w + F2*C3_e 
      C4 = F1*C4_w + F2*C4_e 
      C5 = F1*C5_w + F2*C5_e 
      C1star = F1*C1star_w + F2*C1star_e
      C3star = F1*C3star_w + F2*C3star_e

      ! Production term
      rhoP = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        rhoP(ii,jj) = rhoP(ii,jj) + rhoRij(ii,kk)*Gradvel(jj,kk)
      enddo
      enddo 
      enddo

      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        rhoP(ii,jj) = rhoP(ii,jj) + rhoRij(jj,kk)*Gradvel(ii,kk)
      enddo
      enddo
      enddo
      rhoP = -rhoP

      ! Dissipation term
      rhoeij = 2d0/3d0 * rhoe * deltaij

      ! Pressure-strain correlation
      aij = rhoRij/(rho*kappa) - 2d0/3d0 * deltaij
      aklakl = Sum (aij*aij) 

      ! C1 term
      rhoPkk = rhoP(1,1) + rhoP(2,2) + rhoP(3,3)
      C1_term = -( C1*rhoe + 0.5d0*C1star*rhoPkk )*aij

      ! C2 term
      aikakj = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        aikakj(ii,jj) = aikakj(ii,jj) + aij(ii,kk)*aij(kk,jj)
      enddo
      enddo
      enddo
      C2_term = C2*rhoe * ( aikakj - 1d0/3d0*aklakl*deltaij )

      ! C3 term
      C3_term = ( C3 - C3star * Sqrt( aklakl ) ) * rho*kappa*Sij_star

      ! C4 term
      C4_term = 0d0
      aklSkl = Sum ( aij*Sij )
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C4_term(ii,jj) = C4_term(ii,jj) + aij(ii,kk)*Sij(jj,kk)
      enddo
      enddo
      enddo

      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C4_term(ii,jj) = C4_term(ii,jj) + aij(jj,kk)*Sij(ii,kk) 
      enddo
      enddo
      enddo

      C4_term = C4_term - 2d0/3d0*aklSkl*deltaij
      C4_term = C4*rho*kappa * C4_term

      ! C5 term
      C5_term = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C5_term(ii,jj) = C5_term(ii,jj) + aij(ii,kk)*Wij(jj,kk)
      enddo
      enddo
      enddo

      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        C5_term(ii,jj) = C5_term(ii,jj) + aij(jj,kk)*Wij(ii,kk)
      enddo
      enddo
      enddo
      C5_term = C5*rho*kappa * C5_term

      rhoPI = C1_term + C2_term + C3_term + C4_term + C5_term

      ! Cross Diffusion term
      Diff = sigmad * rho/omega * Max( Dot_Product ( Grad(8,:), Grad(7,:) ), 0d0 )

      ! Source terms
      Source(1) = rhoP(1,1) + rhoPI(1,1) - rhoeij(1,1) 
      Source(2) = rhoP(2,2) + rhoPI(2,2) - rhoeij(2,2) 
      Source(3) = rhoP(3,3) + rhoPI(3,3) - rhoeij(3,3)
      Source(4) = rhoP(1,2) + rhoPI(1,2) - rhoeij(1,2)
      Source(5) = rhoP(1,3) + rhoPI(1,3) - rhoeij(1,3)
      Source(6) = rhoP(2,3) + rhoPI(2,3) - rhoeij(2,3)
      Source(7) = alphaw*omega/kappa*0.5d0*rhoPkk - betaw*rho*omega**2 + Diff
      ! print*, i,j,k
      ! print*, source
      
      Res(nt:nt+6,i,j,k) = Res(nt:nt+6,i,j,k) - Source * Volume(i,j,k)

    enddo ; enddo ; enddo ! (i, j, k) loop

  end subroutine SSGLRR_Blk

  
  subroutine SSGLRR_Eddy_Viscosity ( mi_t, var, mi_l, rho, Gradvel, dist, ks )
    use ARES_Global_m
    use ARES_Lib_Fluid
    implicit none
    real(kind=8), intent(in), dimension(nRANS) :: var        ! : [rho*Rij rho*w]
    real(kind=8), intent(in)                   :: mi_l       ! : Molecular viscosity
    real(kind=8), intent(in)                   :: rho        ! : Density
    real(kind=8), intent(in)                   :: dist       ! : Wall distance (not used)
    real(kind=8), intent(in), dimension(3,3)   :: Gradvel    ! : Velocity gradient
    real(kind=8), intent(in)                   :: ks         ! : sand grain roughenss of nearest wall ( not used)
    real(kind=8), intent(out)                  :: mi_t       ! : Eddy viscosity
    ! Local
    real(kind=8) :: rhok, w

    w = var(7) / rho
    rhok = 0.5d0 * Sum( var(1:3) )
    mi_t = rhok / w
  end subroutine SSGLRR_Eddy_Viscosity


  subroutine SSGLRR_RANS_Diffusive_Flux ( Flux, var, mi_l, vGrad, Grad, rho, Area, Normal, Dist )
    use ARES_Global_m
    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: var                    ! : [rho*Rij rho*w]
    real(kind=8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(kind=8), intent(in), dimension(3,3)     :: vGrad                  ! : Velocity gradient (not used)
    real(kind=8), intent(in), dimension(nRANS,3) :: Grad                   ! : Gradient of [ Rij w ]
    real(kind=8), intent(in)                     :: rho                    ! : Density
    real(kind=8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(kind=8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    integer :: i
    real(kind=8) :: omega, kappa, Gradk(3), F1, F2, D_, sigmaw, Diff(6,3)

    omega = var(7)/rho
    kappa = 0.5d0 * Sum ( var(1:3) ) / rho
    Gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) ) 

    ! Blending of coefficients
    F1 = compute_F1 ( rho=rho, w=omega, k=kappa, mil=mi_l, gradw=Grad(7,:), gradk=Gradk, d=Dist )
    F2 = 1d0 - F1
    D_ = ( D_w * F1 + D_e * F2 ) / ( C_mu*omega )
    sigmaw = F1 * sigmaw_w + F2 * sigmaw_e

    Diff(:,1) = mi_l*Grad(1:6,1) + D_*( var(1)*Grad(1:6,1) + var(4)*Grad(1:6,2) + var(5)*Grad(1:6,3) )
    Diff(:,2) = mi_l*Grad(1:6,2) + D_*( var(4)*Grad(1:6,1) + var(2)*Grad(1:6,2) + var(6)*Grad(1:6,3) )
    Diff(:,3) = mi_l*Grad(1:6,3) + D_*( var(5)*Grad(1:6,1) + var(6)*Grad(1:6,2) + var(3)*Grad(1:6,3) )

    do i = 1, 6
      Flux(i) = Dot_Product ( [ Diff(i,1), Diff(i,2), Diff(i,3) ], Normal )
    end do

    Flux(7) = ( mi_l + sigmaw*rho*kappa/omega ) * Dot_Product ( Grad(7,:), Normal )
    Flux = Flux * Area

  end subroutine SSGLRR_RANS_Diffusive_Flux


  subroutine SSGLRR_SD_RANS_Diffusive_Flux ( Flux, var, mi_l, vGrad, Grad, rho, Area, Normal, Dist )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: var                    ! : [rho*Rij rho*w]
    real(kind=8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(kind=8), intent(in), dimension(3,3)     :: vGrad                  ! : Velocity gradient (not used)
    real(kind=8), intent(in), dimension(nRANS,3) :: Grad                   ! : Gradient of [ Rij w ]
    real(kind=8), intent(in)                     :: rho                    ! : Density
    real(kind=8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(kind=8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    integer :: i
    real(kind=8) :: omega, kappa, Gradk(3), F1, F2, D_, sigmaw, mi_t, Diff(6,3)

    omega = var(7)/rho
    kappa = 0.5d0 * Sum ( var(1:3) ) / rho
    Gradk = 0.5d0 * ( Grad(1,:) + Grad(2,:) + Grad(3,:) ) 

    ! Blending of coefficients
    F1 = compute_F1 (rho=rho, w=omega, k=kappa, mil=mi_l, gradw=Grad(7,:), gradk=Gradk, d=Dist)
    F2 = 1d0 - F1
    D_ = 0.5d0*C_mu*F1 + 2d0/3d0*0.22d0*F2
    sigmaw = F1 * sigmaw_w + F2 * sigmaw_e
    mi_t = rho*kappa/omega

    Diff(1:6,:) = (mi_l + D_/C_mu*mi_t) * Grad(1:6,:)
    do i = 1, 6
      Flux(i) = Dot_Product ( [ Diff(i,1), Diff(i,2), Diff(i,3) ], Normal )
    end do
    Flux(7) = ( mi_l + sigmaw*mi_t ) * Dot_Product ( Grad(7,:), Normal )
    Flux = Flux * Area

  end subroutine SSGLRR_SD_RANS_Diffusive_Flux

  subroutine SSGLRR_Set_Wall_Values ( var_cell, mil_wall, var_wall, metric, dist, ks )
    use ARES_Global_m
    use ARES_Config_Types_m, only: obj_rans

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: var_cell              ! : rans variable at wall
    real(kind=8), intent(in)                     :: mil_wall              ! : laminar viscosity at wall
    real(kind=8), intent(out), dimension(nRANS)  :: var_wall              ! : rans variable at wall
    real(kind=8), intent(in)                         :: metric            ! : MetricTensor(Dir) dot_product FaceNormal
    real(kind=8), intent(in)                         :: dist              ! : cell center distance wall
    real(kind=8), intent(in)                         :: ks                ! : roughness of nearest wall

    var_wall(1:6) = 0d0 ! Solid surface condition on rhoRij
    var_wall(7) = 8d2 * mil_wall / ( dist**2 ) ! approximate BC for smooth surface (Menter kw-SST)

  end subroutine SSGLRR_Set_Wall_Values


  subroutine SSGLRR_Extrapolate_Wall ( var, wall_var, rho, wall_rho, ghost_var )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(in),  dimension(nRANS) :: var, wall_var  ! : boundary cell and wall variables
    real(kind=8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(kind=8), intent(out), dimension(nRANS) :: ghost_var ! : ghost cell RANS variables

    ghost_var(1:6) = - var(1:6) ! enforcing rhoRij_wall == 0
    ghost_var(7) = rho * 2d0 * wall_var(7)/wall_rho - var(7) ! extrapolating from omega_wall

  end subroutine SSGLRR_Extrapolate_Wall


  subroutine SSGLRR_Enforce_Realizability ( rhoRij )

    implicit none
    real(kind=8), intent(inout),  dimension(:) :: rhoRij

    ! 1st condition: u'u', v'v' and w'w' > 0
    if ( rhoRij(1) < 1d-16 ) rhoRij(1) = 1d-16
    if ( rhoRij(2) < 1d-16 ) rhoRij(2) = 1d-16
    if ( rhoRij(3) < 1d-16 ) rhoRij(3) = 1d-16

    ! 2nd condition: abs( u'v' ) < sqrt ( u'u' * v'v' ), etc
    if ( Abs( rhoRij(4) ) > Sqrt( rhoRij(1)*rhoRij(2) ) ) then
      rhoRij(4) = Sign( 1d0, rhoRij(4) ) * Sqrt( rhoRij(1)*rhoRij(2) )
    end if

    if ( Abs( rhoRij(5) ) > Sqrt( rhoRij(1)*rhoRij(3) ) ) then
      rhoRij(5) = Sign( 1d0, rhoRij(5) ) * Sqrt( rhoRij(1)*rhoRij(3) )
    end if

    if ( Abs( rhoRij(6) ) > Sqrt( rhoRij(2)*rhoRij(3) ) ) then
      rhoRij(6) = Sign( 1d0, rhoRij(6) ) * Sqrt( rhoRij(2)*rhoRij(3) )
    end if

  end subroutine SSGLRR_Enforce_Realizability


  function compute_F1 ( rho, w, k, mil, Gradw, Gradk, d ) result ( F1 )

    implicit none
    real(kind=8), intent(in) :: rho, w, k, mil, Gradw(3), Gradk(3), d
    real(kind=8) :: CD, zeta, F1

    CD = sigmad_e * rho / w * Max ( Dot_Product( Gradw, Gradk ), 1d-20 )
    zeta = Min ( Max( Sqrt(k)/(C_mu*w*d), 5d2*mil/(rho*w*d**2) ), &
                 4d0*sigmaw_e*rho*k/(CD*d**2) )
    F1 = Tanh ( zeta**4 )

  end function compute_F1
  
  
end module ARES_Lib_SSGLRR