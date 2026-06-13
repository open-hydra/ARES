!> @brief Module for Kappa-Omega SST (2003) model. see https://turbmodels.larc.nasa.gov/sst.html.
module ARES_Lib_SST
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  real(kind=8), parameter, public :: &
  sigma_k1 = 0.85d0, &
  sigma_k2 = 1.0d0, &
  sigma_w1 = 0.5d0, &
  sigma_w2 = 0.856d0, &
  beta_1 = 0.075d0, &
  beta_2 = 0.0828d0, &
  beta_star = 0.09d0, &
  karman = 0.41d0, &
  gamma_1 = 5d0/9d0, &
  gamma_2 = 0.44d0, &
  cr1 = 1.d0, &
  cr2 = 2.d0, &
  cr3 = 1.d0
  real(kind=8), public :: alpha_1 = 0.31d0

contains


  subroutine SST_Source_Terms ( domain )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_rans
    use ARES_Lib_SpalartShur, only : compute_velocity_gradient, compute_rc_terms
    use ARES_Mod_MPI, only: is_local_block

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    logical :: SpalartShur, k_energy_coupling
      
    SpalartShur       = obj_rans%SpalartShur
    k_energy_coupling = obj_rans%k_energy_coupling

    if ( obj_rans%SpalartShur ) then
      call compute_velocity_gradient ( domain )
      call compute_rc_terms ( domain )
    end if
    
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call SST_Blk ( domain % blk(b) % P,            &
                     domain % blk(b) % r,            &
                     domain % blk(b) % M,            &
                     domain % blk(b) % vol,          &
                     domain % blk(b) % yn,           &
                     domain % blk(b) % vel_gradient, &
                     domain % blk(b) % rc_term1,     &
                     domain % blk(b) % rc_term2,     &
                     domain % blk(b) % dim,          &
                     SpalartShur, k_energy_coupling )
    
    end do

  end subroutine SST_Source_Terms


  subroutine SST_Blk ( Prim, Res, M, Volume, WDist, gradv, rc1, rc2, n, SpalartShur, k_energy_coupling )
    use ARES_Base_Types_m
    use ARES_Global_m
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Fluid
    use ARES_Config_Types_m, only: obj_rans

    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: SpalartShur, k_energy_coupling
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: WDist
    type(ARES_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    type(ARES_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: gradv
    real(kind=8), dimension(n(1), n(2), n(3)), intent(in) :: rc1, rc2
    ! Local
    integer :: i, j, k, ii, jj
    real(kind=8) :: rho, mil, kap, ome, dist, Gradvel(3,3), Divel, Sij(3,3), S, Wij(3,3), diag, mi_t, F(2)
    real(kind=8) :: vort(3), O, rstar, D, rtilde, frot, fr1
    real(kind=8) :: Tij(3,3), Prod(2), Grad(2,3), dkDotdw, Diff, beta, gamma, Diss(2), Source(2)

    !$omp do collapse (3) private ( rho, mil, kap, ome, dist, Gradvel, Divel, Sij, S, Wij, diag, mi_t, F ), &
    !$omp private ( Tij, Prod, Grad, dkDotdw, Diff, beta, gamma, Diss, Source, i, j, k, ii, jj )
    
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      ! Local variables
      rho = ph2vars( Prim(np,i,j,k), Prim(nh,i,j,k),  rho_tab )
      mil = ph2vars( Prim(np,i,j,k), Prim(nh,i,j,k),  mi_tab2D )
      kap = Prim(nt,i,j,k) / rho
      ome = Prim(nt+1,i,j,k) / rho
      dist = WDist(i,j,k)

      ! Velocity gradient computation
      if ( obj_rans%SpalartShur ) then
        Gradvel = gradv(i,j,k) % c
      else
        Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
        Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
        Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
        Gradvel = Matmul ( Gradvel, M(i,j,k) % c )
      end if
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)

      ! Gradient of kappa and omega
      Grad(:,1) = ( Prim(nt:nt+1,i+1,j,k) / ph2vars( Prim(np,i+1,j,k), Prim(nh,i+1,j,k),  rho_tab ) - &
                    Prim(nt:nt+1,i-1,j,k) / ph2vars( Prim(np,i-1,j,k), Prim(nh,i-1,j,k),  rho_tab ) ) * 0.5d0
      Grad(:,2) = ( Prim(nt:nt+1,i,j+1,k) / ph2vars( Prim(np,i,j+1,k), Prim(nh,i,j+1,k),  rho_tab ) - &
                    Prim(nt:nt+1,i,j-1,k) / ph2vars( Prim(np,i,j-1,k), Prim(nh,i,j-1,k),  rho_tab ) ) * 0.5d0
      Grad(:,3) = ( Prim(nt:nt+1,i,j,k+1) / ph2vars( Prim(np,i,j,k+1), Prim(nh,i,j,k+1),  rho_tab ) - &
                    Prim(nt:nt+1,i,j,k-1) / ph2vars( Prim(np,i,j,k-1), Prim(nh,i,j,k-1),  rho_tab ) ) * 0.5d0
      Grad = Matmul ( Grad, M(i,j,k) % c )

      ! Blending of coefficients
      F = compute_F ( rho=rho, w=ome, k=kap, mil=mil, gradk=Grad(1,:), gradw=Grad(2,:), d=dist )
      gamma = gamma_1*F(1) + gamma_2*(1d0-F(1))
      beta = beta_1*F(1) + beta_2*(1d0-F(1))

      ! Strain and vorticity
      Sij = Strain_Tensor ( Gradvel )
      S = Sqrt ( 2d0 * sum ( Sij**2 ) )

      ! Eddy viscosity
      vort = vorticity_vector ( gradvel )
      Wij = Vorticity_Tensor ( vort )
      O = Sqrt ( 2d0 * sum ( Wij**2 ) )
      mi_t = rho*kap*alpha_1 / Max ( alpha_1*ome, S*F(2) )

      ! Turbulent stress tensor
      ! Tij = 2d0 * mi_t * Sij
      ! diag = 2d0/3d0*mi_t*Divel + 2d0/3d0*rho*kap
      ! Tij(1,1) = Tij(1,1) - diag
      ! Tij(2,2) = Tij(2,2) - diag
      ! Tij(3,3) = Tij(3,3) - diag

      ! Production term
      ! Prod(1) = 0d0
      ! do ii = 1, 3
      ! do jj = 1, 3
      !   Prod(1) = Prod(1) + Tij(ii,jj)*Gradvel(ii,jj)
      ! end do
      ! end do
      Prod(1) = mi_t * S**2 ! SST-2003m

      if ( obj_rans%SpalartShur ) then
        rstar = rc1(i,j,k)
        D = sqrt( max ( S**2, 0.09d0*ome**2) )
        rtilde = rc2(i,j,k) / ( O*D**3 )
        frot = (1.d0 + cr1)*(2.d0*rstar/(1.d0 + rstar))*(1.d0 - cr3*atan(cr2*rtilde)) - cr1
        fr1 = max ( min(frot, 1.25d0), 0d0 )
        Prod(1) = fr1*Prod(1)
      end if

      Prod(1) = Min( Prod(1), 10d0*beta_star*rho*ome*kap ) ! k-eqn production limiter
      Prod(2) = ( gamma*rho / mi_t ) * Prod(1)

      ! Cross Diffusion term
      dkDotdw = Dot_Product ( Grad(1,:), Grad(2,:) )
      Diff = 2d0*(1d0-F(1))*rho*sigma_w2/ome*dkDotdw

      ! Dissipation term
      Diss(1) = beta_star*rho*kap*ome
      Diss(2) = beta*rho*ome**2

      ! Source terms
      Source(1) = Prod(1) - Diss(1) 
      Source(2) = Prod(2) - Diss(2) + Diff
      
      Res(nt:nt+1,i,j,k) = Res(nt:nt+1,i,j,k) - Source * Volume(i,j,k)
      if (obj_rans%k_energy_coupling) Res(nh,i,j,k) = Res(nh,i,j,k) + Source(1) * Volume(i,j,k)

    enddo ; enddo ; enddo ! (i, j, k) loop

  end subroutine SST_Blk

  
  subroutine SST_Eddy_Viscosity ( mi_t, rkw, mi_l, rho, Gradvel, dist, ks )
    use ARES_Global_m
    use ARES_Lib_Fluid

    implicit none
    real(kind=8), intent(in), dimension(nRANS) :: rkw        ! : rkw = [rho*k rho*w]
    real(kind=8), intent(in)                   :: mi_l       ! : Molecular viscosity
    real(kind=8), intent(in)                   :: rho        ! : Density
    real(kind=8), intent(in)                   :: dist       ! : Wall distance
    real(kind=8), intent(in)                   :: ks         ! : roughness ( not used)
    real(kind=8), intent(in), dimension(3,3)   :: Gradvel    ! : Velocity gradient
    real(kind=8), intent(out)                  :: mi_t       ! : Eddy viscosity
    ! Local
    real(kind=8) :: kap, ome, Sij(3,3), S, F2, vort(3), Wij(3,3), O

    kap = rkw(1) / rho
    ome = rkw(2) / rho
    F2 = compute_F2 ( rho=rho, w=ome, k=kap, mil=mi_l, d=dist )
    Sij = Strain_Tensor ( Gradvel )
    S = Sqrt ( 2d0 * sum ( Sij**2 ) )
    vort = vorticity_vector ( gradvel )
    Wij = Vorticity_Tensor ( vort )
    O = Sqrt ( 2d0 * sum ( Wij**2 ) )
    mi_t = rho*kap*alpha_1 / Max ( alpha_1*ome, S*F2 )

  end subroutine SST_Eddy_Viscosity


  subroutine SST_RANS_Diffusive_Flux ( Flux, rkw, mi_l, Gradvel, Gradkw, rho, Area, Normal, Dist )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: rkw                    ! : [ rho*k rho*w ]
    real(kind=8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(kind=8), intent(in), dimension(3,3)     :: Gradvel                ! : Velocity gradient
    real(kind=8), intent(in), dimension(nRANS,3) :: Gradkw                 ! : Gradient of [ k w ]
    real(kind=8), intent(in)                     :: rho                    ! : Density
    real(kind=8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(kind=8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    real(kind=8) :: kap, ome, F1, sigma_k, sigma_w, mi_t

    kap = rkw(1)/rho
    ome = rkw(2)/rho
    F1 = compute_F1 (rho=rho, w=ome, k=kap, mil=mi_l, &
                     gradk=Gradkw(1,:), gradw=Gradkw(2,:), d=dist)
    sigma_k = sigma_k1*F1 + sigma_k2*(1d0-F1)
    sigma_w = sigma_w1*F1 + sigma_w2*(1d0-F1)
    call SST_Eddy_Viscosity( mi_t=mi_t, rkw=rkw, mi_l=mi_l, rho=rho, &
                             Gradvel=Gradvel, dist=dist, ks=0d0 )
    Flux(1) = ( mi_l + sigma_k*mi_t ) * Dot_Product ( Gradkw(1,:), Normal )
    Flux(2) = ( mi_l + sigma_w*mi_t ) * Dot_Product ( Gradkw(2,:), Normal )
    Flux = Flux * Area

  end subroutine SST_RANS_Diffusive_Flux


  subroutine SST_Set_Wall_Values ( rkw_cell, mil_wall, rkw_wall, metric, dist, ks )
    use ARES_Global_m
    use ARES_Config_Types_m, only: obj_rans

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: rkw_cell              ! : rans variable at wall
    real(kind=8), intent(in)                     :: mil_wall              ! : laminar viscosity at wall
    real(kind=8), intent(out), dimension(nRANS)  :: rkw_wall              ! : rans variable at wall
    real(kind=8), intent(in)                     :: metric                ! : MetricTensor(Dir) dot_product FaceNormal
    real(kind=8), intent(in)                     :: dist                  ! : cell center distance wall
    real(kind=8), intent(in)                     :: ks                    ! : roughness of nearest wall
    
    rkw_wall(1) = 0d0 ! Solid surface condition on k
    rkw_wall(2) = 8d2 * mil_wall / ( dist**2 ) ! approximate BC for smooth surface (Menter kw-SST)

  end subroutine SST_Set_Wall_Values


  subroutine SST_Blowing_Correction ( rho, mil, tau, rkw, mdot, dist )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(out), dimension(nRANS)  :: rkw
    real(kind=8), intent(in)                     :: rho, mil, mdot, dist
    real(kind=8), dimension(3), intent(in)       :: tau
    ! Local
    real(kind=8) :: tauw, uT, v_w, v_plus, SB
    
    rkw(1) = 0d0 ! Solid surface condition on k
    tauw = sqrt( sum(tau**2) )
    uT = sqrt ( tauw / rho )
    v_w = abs(mdot)/rho
    v_plus = v_w / uT
    SB = 25d0 / ( v_plus*( 1d0 + 5d0 * v_plus ) )
    rkw(2) = (rho*uT)**2 * SB / mil ! Wilcox BC for mass injection

  end subroutine SST_Blowing_Correction


  subroutine SST_Blowing_noCorrection ( rho, mil, tau, rkw, mdot, dist )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(out), dimension(nRANS)  :: rkw
    real(kind=8), intent(in)                     :: rho, mil, mdot, dist
    real(kind=8), dimension(3), intent(in)       :: tau
    
    rkw(1) = 0d0
    rkw(2) = 8d2 * mil / ( dist**2 )

  end subroutine SST_Blowing_noCorrection


  subroutine SST_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use ARES_Global_m
    implicit none
    real(kind=8), intent(in),  dimension(nRANS) :: prim, wall_prim  ! : boundary cell and wall variables
    real(kind=8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(kind=8), intent(out), dimension(nRANS) :: ghost_prim ! : ghost cell RANS variables

    ghost_prim(1) = - prim(1) ! enforcing k_wall == 0
    ghost_prim(2) = rho * 2d0 * wall_prim(2)/wall_rho - prim(2) ! extrapolating from omega_wall

  end subroutine SST_Extrapolate_Wall


  subroutine SST_Enforce_Realizability ( rhokrhow )

    implicit none
    real(kind=8), intent(inout),  dimension(:) :: rhokrhow
    if ( rhokrhow(1) < 1d-10 ) rhokrhow(1) = 1d-10
    if ( rhokrhow(2) < 1d-10 ) rhokrhow(2) = 1d-10

  end subroutine SST_Enforce_Realizability


  function compute_F ( rho, w, k, mil, Gradw, Gradk, d ) result ( F )
    implicit none
    real(kind=8), intent(in) :: rho, w, k, mil, Gradw(3), Gradk(3), d
    real(kind=8) :: CD, arg1, arg2, F(2)

    CD = Max ( 2d0*rho*sigma_w2/w*Dot_Product( Gradw, Gradk ), 1d-10 )
    arg1 = Min ( Max ( Sqrt(k)/(beta_star*w*d), &
                       5d2*mil/(rho*w*d**2) ), &
                       4d0*sigma_w2*rho*k/(CD*d**2) )
    arg2 = Max ( 2d0*Sqrt(k)/(beta_star*w*d), 5d2*mil/(rho*w*d**2) )
    F(1) = Tanh ( arg1**4 )
    F(2) = Tanh ( arg2**2 )
  end function compute_F

  function compute_F1 ( rho, w, k, mil, Gradw, Gradk, d ) result ( F1 )
    implicit none
    real(kind=8), intent(in) :: rho, w, k, mil, Gradw(3), Gradk(3), d
    real(kind=8) :: CD, arg1, F1

    CD = Max ( 2d0*rho*sigma_w2/w*Dot_Product( Gradw, Gradk ), 1d-10 )
    arg1 = Min ( Max ( Sqrt(k)/(beta_star*w*d), &
                       5d2*mil/(rho*w*d**2) ), &
                 4d0*sigma_w2*rho*k/(CD*d**2) )
    F1 = Tanh ( arg1**4 )
  end function compute_F1

  function compute_F2 ( rho, w, k, mil, d ) result ( F2 )
    implicit none
    real(kind=8), intent(in) :: rho, w, k, mil, d
    real(kind=8) :: arg2, F2

    arg2 = Max ( 2d0*Sqrt(k)/(beta_star*w*d), 5d2*mil/(rho*w*d**2) )
    F2 = Tanh ( arg2**2 )

  end function compute_F2
  
end module ARES_Lib_SST