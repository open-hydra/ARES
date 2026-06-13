!> @brief Module for Wilcox 2006 Kappa-Omega model. see https://turbmodels.larc.nasa.gov/wilcox.html.
module ARES_Lib_Wilcox2006
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  real(kind=8), parameter, public :: &
  sigma_k = 0.6d0, &
  sigma_w = 0.5d0, &
  beta_star = 0.09d0, &
  beta_0 = 0.0708d0, &
  gamma = 13d0/25d0, &
  Clim = 7d0/8d0

contains


  subroutine Wilcox2006_Source_Terms ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Mod_MPI, only: is_local_block
    use ARES_Config_Types_m, only: obj_rans

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b
    logical :: k_energy_coupling
    
    k_energy_coupling = obj_rans%k_energy_coupling

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle

      call Wilcox_Blk ( domain % blk(b) % P,    &
                        domain % blk(b) % r,    &
                        domain % blk(b) % M,    &
                        domain % blk(b) % vol,  &
                        domain % blk(b) % dim,  &
                        k_energy_coupling       )
    
    end do

  end subroutine Wilcox2006_Source_Terms


  subroutine Wilcox_Blk ( Prim, Res, M, Volume, n, k_energy_coupling )
    use ARES_Base_Types_m
    use ARES_Global_m
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Fluid
    !use ARES_Lib_RotatingFrame, only: obj_rot

    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: k_energy_coupling
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(kind=8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    real(kind=8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Volume
    type(ARES_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
    ! Local
    integer :: i, j, k, ii, jj, kk
    real(kind=8) :: rho, kap, ome, Gradvel(3,3), Divel, Sij(3,3), Sij_bar(3,3), diag, ome_hat, mi_t, Tij(3,3), Prod(2)
    real(kind=8) :: Gradkw(2,3), dkDotdw, sigma_d, Diff, Sij_hat(3,3), Wij(3,3), Num, X_w, f_beta, beta, Diss(2), Source(2)

    !$omp do collapse (3) private (rho, kap, ome, Gradvel, Divel, Sij, Sij_bar, diag, ome_hat, mi_t, Tij, Prod), &
    !$omp private ( Gradkw, dkDotdw, sigma_d, Diff, Sij_hat, Wij, Num, X_w, f_beta, beta, Diss, Source ), &
    !$omp private ( i, j, k, ii, jj, kk )

    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      ! Local variables
      !rho = Sum ( Prim(1:np,i,j,k) )
      rho = ph2vars( Prim(np,i,j,k), Prim(nh,i,j,k), rho_tab )
      kap = Prim(nt,i,j,k) / rho
      ome = Prim(nt+1,i,j,k) / rho

      ! Velocity gradient computation
      Gradvel(:,1) = ( Prim(nu:nw,i+1,j,k) - Prim(nu:nw,i-1,j,k) ) * 0.5d0
      Gradvel(:,2) = ( Prim(nu:nw,i,j+1,k) - Prim(nu:nw,i,j-1,k) ) * 0.5d0
      Gradvel(:,3) = ( Prim(nu:nw,i,j,k+1) - Prim(nu:nw,i,j,k-1) ) * 0.5d0
      Gradvel = Matmul ( Gradvel, M(i,j,k) % c )

      ! Strain tensor and similar stuff
      Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3) 
      Sij = Strain_Tensor ( Gradvel )
      Sij_bar = Sij
      diag = 1d0/3d0 * Divel
      Sij_bar(1,1) = Sij_bar(1,1) - diag
      Sij_bar(2,2) = Sij_bar(2,2) - diag
      Sij_bar(3,3) = Sij_bar(3,3) - diag

      ! Eddy viscosity
      ome_hat = Max ( ome, Clim * Sqrt( 2d0 * Sum ( Sij_bar**2 ) / beta_star ) )
      mi_t = rho * kap / ome_hat

      ! Turbulent stress tensor
      Tij = 2d0 * mi_t * Sij
      diag = 2d0/3d0 * mi_t * Divel + 2d0/3d0 * rho * kap
      Tij(1,1) = Tij(1,1) - diag
      Tij(2,2) = Tij(2,2) - diag
      Tij(3,3) = Tij(3,3) - diag

      ! Production term
      Prod(1) = 0d0
      do ii = 1, 3
      do jj = 1, 3
        Prod(1) = Prod(1) + Tij(ii,jj) * Gradvel(ii,jj)
      end do
      end do
      Prod(2) = ( gamma * ome / kap ) * Prod(1)
      Prod(1) = Min( Prod(1), 20d0*beta_star*rho*ome*kap ) ! k-eqn production limiter

      ! Gradient of kap and ome
      Gradkw(:,1) = ( Prim(nt:nt+1,i+1,j,k) / ph2vars( Prim(np,i+1,j,k), Prim(nh,i+1,j,k), rho_tab ) - &
                      Prim(nt:nt+1,i-1,j,k) / ph2vars( Prim(np,i-1,j,k), Prim(nh,i-1,j,k), rho_tab ) ) * 0.5d0
      Gradkw(:,2) = ( Prim(nt:nt+1,i,j+1,k) / ph2vars( Prim(np,i,j+1,k), Prim(nh,i,j+1,k), rho_tab ) - &
                      Prim(nt:nt+1,i,j-1,k) / ph2vars( Prim(np,i,j-1,k), Prim(nh,i,j-1,k), rho_tab ) ) * 0.5d0
      Gradkw(:,3) = ( Prim(nt:nt+1,i,j,k+1) / ph2vars( Prim(np,i,j,k+1), Prim(nh,i,j,k+1), rho_tab ) - &
                      Prim(nt:nt+1,i,j,k-1) / ph2vars( Prim(np,i,j,k-1), Prim(nh,i,j,k-1), rho_tab ) ) * 0.5d0
      Gradkw = Matmul ( Gradkw, M(i,j,k) % c )

      ! Cross Diffusion term
      dkDotdw = Dot_Product ( Gradkw(1,:), Gradkw(2,:) )
      sigma_d = 1d0/8d0
      if ( dkDotdw <= 0d0 ) sigma_d = 0d0
      Diff = sigma_d * rho / ome * dkDotdw

      ! Dissipation term
      Sij_hat(1,1) = Sij(1,1) - 0.5d0 * Divel 
      Sij_hat(2,2) = Sij(2,2) - 0.5d0 * Divel 
      Sij_hat(3,3) = Sij(3,3) - 0.5d0 * Divel 

      Wij = Vorticity_Tensor ( Vorticity_Vector ( Gradvel ) )

      Num = 0d0
      do ii = 1, 3
      do jj = 1, 3
      do kk = 1, 3
        Num = Num + Wij(ii,jj) * Wij(jj,kk) * Sij_hat(kk,ii)
      end do
      end do
      end do

      X_w = Abs ( Num / ( beta_star * ome )**3 )
      f_beta = ( 1d0 + 85d0 * X_w ) / ( 1d0 + 100d0 * X_w )
      beta = beta_0 * f_beta

      Diss(1) = beta_star * rho * kap * ome
      Diss(2) = beta * rho * ome**2

      ! Source terms
      Source(1) = Prod(1) - Diss(1) 
      Source(2) = Prod(2) - Diss(2) + Diff
      
      Res(nt:nt+1,i,j,k) = Res(nt:nt+1,i,j,k) - Source * Volume(i,j,k)
      if (k_energy_coupling) Res(nh,i,j,k) = Res(nh,i,j,k) + Source(1) * Volume(i,j,k)

    enddo ; enddo ; enddo ! (i, j, k) loop

  end subroutine Wilcox_Blk

  
  subroutine Wilcox2006_Eddy_Viscosity ( mi_t, rkw, mi_l, rho, Gradvel, dist, ks )
    use ARES_Global_m
    use ARES_Lib_Fluid

    implicit none
    real(kind=8), intent(in), dimension(nRANS) :: rkw        ! : rkw = [rho*k rho*w]
    real(kind=8), intent(in)                   :: mi_l       ! : Molecular viscosity
    real(kind=8), intent(in)                   :: rho        ! : Density
    real(kind=8), intent(in)                   :: dist       ! : Wall distance (not used)
    real(kind=8), intent(in)                   :: ks       ! : roughness (not used)
    real(kind=8), intent(in), dimension(3,3)   :: Gradvel    ! : Velocity gradient
    real(kind=8), intent(out)                  :: mi_t       ! : Eddy viscosity
    ! Local
    real(kind=8) :: w, Divel, Sij(3,3), diag, w_hat

    w = rkw(2) / rho
    Divel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3) 
    Sij = Strain_Tensor ( Gradvel )
    diag = 1d0/3d0 * Divel
    Sij(1,1) = Sij(1,1) - diag
    Sij(2,2) = Sij(2,2) - diag
    Sij(3,3) = Sij(3,3) - diag
    w_hat = Max ( w, Clim * Sqrt( 2d0 * Sum( Sij**2 ) / beta_star ) )
    mi_t = rkw(1) / w_hat

  end subroutine Wilcox2006_Eddy_Viscosity


  subroutine Wilcox2006_RANS_Diffusive_Flux ( Flux, rkw, mi_l, vGrad, Gradkw, rho, Area, Normal, Dist )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: rkw                    ! : [ rho*k rho*w ]
    real(kind=8), intent(in)                     :: mi_l                   ! : laminar viscosity
    real(kind=8), intent(in), dimension(3,3)     :: vGrad                  ! : Velocity gradient (not used)clear
    real(kind=8), intent(in), dimension(nRANS,3) :: Gradkw                 ! : Gradient of [ k w ]
    real(kind=8), intent(in)                     :: rho                    ! : Density
    real(kind=8), intent(in)                     :: Area, Normal(3), Dist  ! : Metrics
    real(kind=8), intent(out), dimension(nRANS)  :: Flux                   ! : Diffusive flux
    ! Local
    real(kind=8) :: rhok, w

    rhok = rkw(1)
    w = rkw(2)/rho
    Flux(1) = ( mi_l + sigma_k * rhok / w ) * Dot_Product ( Gradkw(1,:), Normal )
    Flux(2) = ( mi_l + sigma_w * rhok / w ) * Dot_Product ( Gradkw(2,:), Normal )
    Flux = Flux * Area

  end subroutine Wilcox2006_RANS_Diffusive_Flux


  subroutine Wilcox2006_Set_Wall_Values ( rkw_cell, mil_wall, rkw_wall, metric, dist, ks )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(in), dimension(nRANS)   :: rkw_cell              ! : rans variable at wall
    real(kind=8), intent(in)                     :: mil_wall              ! : laminar viscosity at wall
    real(kind=8), intent(out), dimension(nRANS)  :: rkw_wall              ! : rans variable at wall
    real(kind=8), intent(in)                     :: metric                ! : MetricTensor(Dir) dot_product FaceNormal
    real(kind=8), intent(in)                     :: dist                  ! : cell center distance wall
    real(kind=8), intent(in)                     :: ks                    ! : roughness of nearest wall
    
    rkw_wall(1) = 0d0 ! Solid surface condition on k
    rkw_wall(2) = 8d2 * mil_wall / ( dist**2 ) ! approximate BC for smooth surface (Menter kw-SST)

  end subroutine Wilcox2006_Set_Wall_Values

  subroutine Wilcox2006_Extrapolate_Wall ( prim, wall_prim, rho, wall_rho, ghost_prim )
    use ARES_Global_m

    implicit none
    real(kind=8), intent(in),  dimension(nRANS) :: prim, wall_prim  ! : boundary cell and wall variables
    real(kind=8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
    real(kind=8), intent(out), dimension(nRANS) :: ghost_prim ! : ghost cell RANS variables

    ghost_prim(1) = - prim(1) ! enforcing k_wall == 0
    ghost_prim(2) = rho * 2d0 * wall_prim(2)/wall_rho - prim(2) ! extrapolating from omega_wall

  end subroutine Wilcox2006_Extrapolate_Wall


  subroutine Wilcox2006_Enforce_Realizability ( rhokrhow )

    implicit none
    real(kind=8), intent(inout),  dimension(:) :: rhokrhow
    if ( rhokrhow(1) < 1d-10 ) rhokrhow(1) = 1d-10
    if ( rhokrhow(2) < 1d-10 ) rhokrhow(2) = 1d-10

  end subroutine Wilcox2006_Enforce_Realizability
  
  
end module ARES_Lib_Wilcox2006