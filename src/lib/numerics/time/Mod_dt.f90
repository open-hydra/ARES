module ARES_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_dt, Set_Global_dt

contains

  !> Local time step of every cell (blk%dtlocal) and its minimum over the domain (domain%dtglobal, on input an upper
  !> bound). In time-accurate RANS runs the turbulence sources are explicit and the step includes a bound of their
  !> stiffness from the current state of the cell (RANS_Stiffness_Bound). Call it outside OpenMP parallel regions.
  subroutine Compute_dt ( domain, cfl, vnn, rampa_iter )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Lib_RANS
    use ARES_Config_Types_m, only: obj_time_scheme, obj_sim_param, obj_rans
    use FLINT_Lib_Thermodynamic, only: ph2vars, rho_tab
    use ARES_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8
    use ARES_Lib_Metrics, only: delthe
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: cfl, vnn
    integer, intent(in)  :: rampa_iter
    ! Local
    integer  :: i, j, k, b
    real(R8) :: dtcell, dtglobal, dtglobal_mpi, dstiff, ks
    real(R8) :: rc, er(3), cen(3)
    logical  :: explicit_sources, rough, wedge

    dtglobal = domain % dtglobal

    ! Time-accurate RANS runs (time-accurate, or a time-accurate HYDRA coupling): the turbulence sources are explicit
    ! (the point-implicit factor serves steady runs only), and a bound of their stiffness on the current state of each
    ! cell enters the step. The roughness of the nearest wall counts with SA-rough only, as in the source terms
    explicit_sources = model == 2 .and. ( obj_time_scheme%time_accurate .or. obj_sim_param%HYDRA_time_accurate )
    rough = obj_rans%rough

    ! 2-D axisymmetric wedge (every 2-D mesh with a wedge angle, whatever the type of its k faces): the viscous stress on
    ! the wedge faces enters the rows of u_r and u_theta with terms that grow as 1/r^2 next to the axis, and the step
    ! must see them (compute, s_h). Faces of type 200 (Lib_BC_Fluxes_Rotational) and symmetry faces 300
    ! (Lib_BC_Fluxes_Symmetry: without swirl the mirror image of the cell is the cell turned by the wedge angle) carry
    ! that stress; null k faces carry none, and there s_h only makes the step smaller than needed. Same guard as the
    ! stress of the faces of type 200: delthe is negative on a right-handed mesh with the axis on the j = jm face.
    ! delthe exactly zero written as two ordered comparisons (no real equality test).
    wedge = ndir == 2 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) .and. model > 0

    do b = 1, domain % nb ! loop over blocks
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3) private ( dtcell, dstiff, ks, rc, er, cen ), reduction ( min : dtglobal )
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        dstiff = 0d0
        if ( explicit_sources ) then
          ks = 0d0
          if ( rough ) ks = domain%blk(b)%ks(i,j,k)
          dstiff = RANS_Stiffness_Bound ( domain%blk(b)%P(nt:nt+nRANS-1,i,j,k), &
                                          ph2vars( domain%blk(b)%P(np,i,j,k), domain%blk(b)%P(nh,i,j,k), rho_tab ), &
                                          domain%blk(b)%yn(i,j,k), ks )
        end if

        ! Wedge: distance of the cell centre (mean of the nodes) from the x axis and radial unit vector
        rc = 0d0
        er = 0d0
        if ( wedge ) then
          cen = 0.125d0 * ( domain%blk(b)%node(i  ,j  ,k  )%c + domain%blk(b)%node(i-1,j  ,k  )%c + &
                            domain%blk(b)%node(i  ,j-1,k  )%c + domain%blk(b)%node(i-1,j-1,k  )%c + &
                            domain%blk(b)%node(i  ,j  ,k-1)%c + domain%blk(b)%node(i-1,j  ,k-1)%c + &
                            domain%blk(b)%node(i  ,j-1,k-1)%c + domain%blk(b)%node(i-1,j-1,k-1)%c )
          rc = sqrt ( cen(2)**2 + cen(3)**2 )
          if ( rc > 0d0 ) er = [ 0d0, cen(2), cen(3) ] / rc
        end if

        ! Compute local cell dt according to CFL and VNN numbers
        call compute ( pressure = domain%blk(b)%p(np,i,j,k), &
                       vel = domain%blk(b)%p(nu:nw,i,j,k), &
                       entalpy = domain%blk(b)%p(nh,i,j,k), &
                       rans_ = domain%blk(b)%P(nt:nprim,i,j,k), &
                       met = domain%blk(b)%m(i,j,k)%c, &
                       rc = rc, er = er, &
                       stiff = dstiff, &
                       dtmin = dtcell, &
                       cfl = cfl, &
                       vnn = vnn )

        ! Apply CFL reduction if required
        if ( domain%iter < rampa_iter ) dtcell = dtcell * domain%iter / rampa_iter

        ! Update local cell dt and global minimum dt
        domain % blk(b) % dtlocal(i,j,k) = dtcell
        dtglobal = min ( dtcell, dtglobal )

      enddo; enddo; enddo
      !$omp end parallel
    enddo ! end of loop over blocks

    ! MPI: global minimum across all ranks
    call mpi_allreduce_min_r8(dtglobal, dtglobal_mpi)
    domain % dtglobal = dtglobal_mpi

    contains

      ! Local time step of a cell, one additive rule for convection and diffusion (Hundsdorfer 2000, notes on the
      ! numerical solution of advection-diffusion-reaction equations, sec. 9.3 and 6; Arnone, Liou, Povinelli 1993,
      ! NASA TM-106373, eq. 26-28):
      !   1/dt = sum_d lambda_d / cfl  +  D_max sum_d |grad xi_d|^2 / vnn,   lambda_d = |u . grad xi_d| + c |grad xi_d|,
      ! with grad xi_d the row d of the cell metric (the normal of the faces in direction d over the cell width) and D_max
      ! the largest diffusivity of the equations. Without dimensional splitting the convective limits of the directions
      ! add up (an explicit scheme on square 2-D cells needs the sum, not the largest of the two); the convective and the
      ! diffusive rates add up as well. For the linearised MUSCL-upwind (Fromm's scheme: the van Leer slope, and the van
      ! Albada and MC ones, in smooth flow; not minmod and superbee) and central-diffusion operators the rule is stable for
      ! cfl <= 1 and vnn <= 1/2 with Heun (RK2), cfl <= 1.17 and vnn <= 0.62 with RK3 (von Neumann limits, proven; at the
      ! worst cell Peclet number the additive law is exact). In 3-D, on a cell skewed in all three pairs of directions, the
      ! viscous rate takes the factor (2 + m + 1/m)/4 below, which bounds the excess of the spectral radius of the viscous
      ! operator over 4 D_max sum_d |grad xi_d|^2 there.
      ! On the 2-D axisymmetric wedge (rc > 0: distance of the cell centre from the x axis, er the radial unit vector) the
      ! stress on the wedge faces adds to the rows of u_r and u_theta diagonal terms and couplings to the neighbours,
      ! per unit mass and with nu_e = (mu + mu_t)/rho, lambda = -2/3 mu, M(d,.) = grad xi_d:
      !   u_r      -(4/3) nu_e/rc^2;  u_r of the neighbours d+-1: -+(2/3) nu_e M(d,.).er/(2 rc);
      !            u_x of the neighbours d+-1: +-(2/3) nu_e M(d,x)/(2 rc)   (from lambda div u in tau_thetatheta);
      !   u_theta  -nu_e/rc^2;  u_theta of the neighbours d+-1: +-nu_e M(d,.).er/(2 rc);   u_x: none.
      ! They raise the Gershgorin row sum of the viscous operator by at most s_h = max(s_r, s_t),
      !   s_r = nu_e [ (4/3)/rc^2 + (2/3) sum_d ( |M(d,x)| + |M(d,.).er| ) / rc ],
      !   s_t = nu_e [ 1/rc^2 + sum_d |M(d,.).er| / rc ],
      ! which enters the diffusive rate as s_h/4, so that dt (4 D_max sum_d |grad xi_d|^2 + s_h) = 4 vnn. Support of the
      ! limits: away from the axis, on a uniform mesh of equal cells, the viscous operator has a real symmetric Fourier
      ! symbol, so it is normal; its spectral radius is at most 4 D_max sum_d |grad xi_d|^2 in 2-D for any skew and aspect
      ! ratio, and in 3-D where the cosines a, b, c between the grad xi_d give a^2 + b^2 + c^2 + 2abc <= 1, elsewhere within
      ! the factor (2 + m + 1/m)/4 < 9/8 that the rate takes in 3-D; so vnn = 1/2 (Heun) and 0.62 (RK3) from their real
      ! stability intervals (Hundsdorfer 2000, sec. 2.4, eq. 2.10, and sec. 3.3, pp. 25-28). Next to the axis, on equal
      ! rectangular cells, the scalars are stable on every row (self-adjoint in the volume-weighted inner product, row sums
      ! within the rate); the momentum operator is not self-adjoint, and its support is evidence (the eigenvalues of the
      ! operator for dx/dr 0.01 to 128, and the largest stable vnn of a diffusion-dominated 1-degree wedge, dx/dr 1 to 8: at
      ! least 0.50 with Heun and 0.63 with RK3).
      ! The ghost cell of the axis counts with its modulus, so the bound holds whatever the condition on the axis; next to
      ! the axis of square cells s_h = 8 nu/dr^2 (the folded exact sum is 20/3).
      ! The real fluid: cv = cp/gamma, gamma = c^2 (drho/dp)_T, with (drho/dp)_T = (drho/dp)_h + (drho/dh)_p (dh/dp)_T and
      ! (dh/dp)_T = 1/rho + T (drho/dT)_p / rho^2 from the tables (rp_tab, rh_tab, dT_tab); for a perfect gas cv = cp - R.
      ! With explicit turbulence sources (time-accurate RANS runs) a bound d of the stiffness of their destruction, at
      ! least its derivative on the state of the cell (RANS_Stiffness_Bound), adds to the rates with the constant 1:
      ! 1/dt += d, so that dt d <= 1. A forward-Euler stage then keeps nu~, k and omega positive, and so do Heun and
      ! RK3, convex combinations of forward-Euler stages (Hundsdorfer 2000, sec. 6.2, theorem 6.4 and the explicit
      ! trapezoidal rule; C.-W. Shu, S. Osher, J. Comput. Phys. 77, 1988).
      subroutine compute ( pressure, vel, entalpy, rans_, met, rc, er, stiff, dtmin, cfl, vnn )
        use ARES_Global_m
        use FLINT_Lib_Thermodynamic
        use ARES_Lib_Preconditioning
        use ARES_Config_Types_m, only: obj_prec

        implicit none
        real(R8), intent(in)  :: pressure, vel(3), entalpy, rans_(:)
        real(R8), intent(in)  :: met(3,3), cfl, vnn
        real(R8), intent(in)  :: rc, er(3)  ! wedge: distance of the cell centre from the x axis (else 0), radial direction
        real(R8), intent(in)  :: stiff    ! bound of the stiffness of the explicit turbulence sources, 0 without [1/s]
        real(R8), intent(out) :: dtmin
        ! Local
        real(R8), parameter :: sigma_SA = 2d0/3d0   ! diffusion constant of Spalart-Allmaras (Lib_Spalart.f90)
        integer  :: d
        real(R8) :: rho, Sound, Alpha, Beta, U_r, c_star, gn, un, lambda
        real(R8) :: rate_c, rate_v, temperature, cp, cv, drdp_T, mil, kl, mie, Dmax
        real(R8) :: nu_e, sum_x, sum_r, s_h
        real(R8) :: dummy(3,3)
        real(R8) :: cos12, cos13, cos23, ssum, qprod, mroot

        Sound = ph2vars( pressure, entalpy, sound_tab )
        rho   = ph2vars( pressure, entalpy, rho_tab   )
        if ( obj_prec%enabled ) then
          U_r = comp_Ur ( vel=Norm2( Vel ), sound=Sound )
          Beta = 1.d0/Sound**2
          if ( Uref > 0.d0 ) Beta = 1.d0/Uref**2 ! in case of reference sound velocity
          Alpha = 5.d-1 * ( 1.d0 - Beta * U_r**2 )
        endif

        ! Convective spectral radii and squared metrics, summed over the directions
        rate_c = 0d0
        rate_v = 0d0
        do d = 1, ndir
          gn = norm2 ( met(d,:) )                          ! |grad xi_d|: 1 / width of the cell across its d faces
          un = abs ( dot_product ( vel, met(d,:) ) ) / gn  ! velocity normal to the d faces

          if ( obj_prec%enabled ) then
            c_star = Sqrt ( Alpha**2 * un**2 + U_r**2 )
            lambda = un * ( 1d0 - Alpha ) + c_star
          else
            lambda = un + Sound
          endif

          rate_c = rate_c + lambda * gn
          rate_v = rate_v + gn**2
        enddo

        if ( model==0 ) then
          dtmin = cfl / rate_c
          return
        endif

        ! 3-D: on a cell whose three directions are skewed against each other, so that the cosines cos12, cos13, cos23
        ! between the rows grad xi_d of the metric give cos12^2 + cos13^2 + cos23^2 + 2 cos12 cos13 cos23 > 1, the
        ! spectral radius of the viscous operator exceeds 4 D_max sum_d |grad xi_d|^2 by a factor below 9/8. The rate
        ! takes the bound (2 + m + 1/m)/4 of that factor, m + 1 the largest eigenvalue of the matrix of the cosines (m the
        ! largest root of m^3 - s m - 2 cos12 cos13 cos23 = 0, s the sum of the squares of the cosines): exact on
        ! rhombohedra, and not taken where the rate already holds (orthogonal cells, and cells with one direction
        ! orthogonal to the other two)
        if ( ndir == 3 ) then
          cos12 = dot_product ( met(1,:), met(2,:) ) / ( norm2 ( met(1,:) ) * norm2 ( met(2,:) ) )
          cos13 = dot_product ( met(1,:), met(3,:) ) / ( norm2 ( met(1,:) ) * norm2 ( met(3,:) ) )
          cos23 = dot_product ( met(2,:), met(3,:) ) / ( norm2 ( met(2,:) ) * norm2 ( met(3,:) ) )
          ssum  = cos12**2 + cos13**2 + cos23**2
          qprod = cos12 * cos13 * cos23
          if ( ssum + 2d0 * qprod > 1d0 ) then
            mroot  = 2d0 * sqrt ( ssum / 3d0 ) * &
                     cos ( acos ( max ( -1d0, min ( 1d0, 3d0 * sqrt ( 3d0 ) * qprod / ssum**1.5d0 ) ) ) / 3d0 )
            rate_v = rate_v * 0.25d0 * ( 2d0 + mroot + 1d0 / mroot )
          end if
        end if

        ! Largest diffusivity: momentum 4/3 (mu + mu_t)/rho; energy, with the conductivity acting on rho cv T,
        ! (k + mu_t cp/Pr_t)/(rho cv) (gamma nu/Pr for a perfect gas); the turbulence variables. The k and omega
        ! diffusivities of SST and Wilcox 2006, (mu + sigma mu_t)/rho with sigma <= 1 and mu_t = rho k/omega at zero
        ! strain, stay below the momentum one.
        temperature = ph2vars( pressure, entalpy, T_tab  )
        cp          = ph2vars( pressure, entalpy, hT_tab )
        drdp_T      = ph2vars( pressure, entalpy, rp_tab ) + ph2vars( pressure, entalpy, rh_tab ) * &
                      ( 1d0 / rho + temperature * ph2vars( pressure, entalpy, dT_tab ) / rho**2 )
        cv          = cp / ( Sound**2 * drdp_T )
        mil = ph2vars( pressure, entalpy, mi_tab2D )
        kl  = ph2vars( pressure, entalpy, k_tab2D  )

        mie = 0.d0
        if (model==2) then
          ! Note: approximate mit for 2 equation models.
          ! velocity gradient assumed 0 and small distance from wall
          dummy = 0d0
          call Eddy_Viscosity ( mut=mie, rans_variables=rans_, mul=mil, rho=rho, &
                                vel_gradient=dummy, walldist=1d-6, ks=0d0 )
        endif

        Dmax = max ( 4d0/3d0 * ( mil + mie ) / rho, ( kl + mie * cp / obj_rans%Prt ) / ( rho * cv ) )

        if ( model==2 .and. nRANS==1 ) Dmax = max ( Dmax, ( mil + rans_(1) ) / ( sigma_SA * rho ) )

        ! Wedge: increase of the row sum from the stress on the wedge faces (see above)
        s_h = 0d0
        if ( rc > 0d0 ) then
          nu_e  = ( mil + mie ) / rho
          sum_x = 0d0
          sum_r = 0d0
          do d = 1, ndir
            sum_x = sum_x + abs ( met(d,1) )
            sum_r = sum_r + abs ( dot_product ( met(d,:), er ) )
          enddo
          s_h = nu_e * max ( 4d0/3d0 / rc**2 + 2d0/3d0 * ( sum_x + sum_r ) / rc, 1d0 / rc**2 + sum_r / rc )
        end if

        dtmin = 1d0 / ( rate_c / cfl + ( Dmax * rate_v + 0.25d0 * s_h ) / vnn )

        if ( stiff > 0d0 ) dtmin = 1d0 / ( 1d0 / dtmin + stiff )

      end subroutine compute

  end subroutine Compute_dt


  subroutine Set_Global_dt ( domain )
    use ARES_Advanced_Types_m
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        domain % blk(b) % dtlocal(i,j,k) = domain % dtglobal
      enddo; enddo; enddo
      !$omp end parallel
    enddo

  end subroutine Set_Global_dt

end module ARES_Mod_dt
