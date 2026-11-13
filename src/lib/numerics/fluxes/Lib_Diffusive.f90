module ARES_Lib_Diffusive
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Diffusive_Flux, Compute_Diffusive_Flux, Face_Metric, Inverse3, Visc_State

contains

    !> Viscous flux of the face between cells 1 and 2, along the normal from 1 to 2: the difference across the face (index
    !> direction a) and the mean central differences of the two cells along the face (b: cells 3..6, c: cells 7..10) of
    !> their viscous variables Visc1..Visc10 (Visc_State, which Viscous_Variables computes once per stage for every cell),
    !> mapped to Cartesian components with the face metric Mf (Fluxes_blk: with diffusive-metric = stencil the metric
    !> stored on the face at setup, blk%dir(d)%Mf of Setup_Stencil_Metrics; with mean (M_L + M_R)/2); the state
    !> of the face is the mean of the primitive variables Prim1 and Prim2 of the two cells. Res1 takes -flux, Res2 +flux.
    !> With a ghost cell of a wall in the stencil, Mg: one face metric per group of quantities (1 species, unused with a
    !> single species, 2 temperature, 3 velocity, 4 turbulence variables), whose ghost values sit at different positions
    !> (Ghost_Wall).
    subroutine Diffusive_Flux ( normal, area, waldis1, waldis2, rough1, rough2, Prim1, Prim2, Visc1, Visc2, Visc3, &
                                Visc4, Visc5, Visc6, Visc7, Visc8, Visc9, Visc10, Mf, Res1, Res2, a, b, c, &
                                Sc, Sct, Prt, Prt_corr, Mg )
      
      use ARES_Global_m

      implicit none
      integer, intent(in) :: a, b, c
      real(R8), intent(in) :: normal(3), area, waldis1, waldis2, rough1, rough2
      real(R8), intent(in), dimension(nprim) :: Prim1, Prim2
      real(R8), intent(in), dimension(nprim) :: Visc1, Visc2, Visc3, Visc4, Visc5, Visc6
      real(R8), intent(in), dimension(nprim) :: Visc7, Visc8, Visc9, Visc10
      real(R8), intent(in), dimension(3,3) :: Mf       ! face metric
      real(R8), intent(in), optional :: Mg(3,3,4)      ! face metrics per group of quantities
      real(R8), intent(inout), dimension(nprim) :: Res1, Res2
      real(R8), intent(in) :: Sc, Sct, Prt
      logical , intent(in) :: Prt_corr
      ! Local
      real(R8) :: Gradient(nprim,3), Prim(nprim), waldis, Flux(nprim), rough
      integer  :: v, g
      real(R8) :: g1, g2, g3

      waldis = 0.0_R8
      rough  = 0.0_R8

      ! Gradient in the same direction of the face: 1 and 2 (velocity, temperature, turbulence variables per unit mass;
      ! no concentration gradient with a single species)
      Gradient ( :, a ) = Visc2 - Visc1
      if (model==2) then
        waldis = 0.5d0 * ( waldis1 + waldis2 ) ! distance to nearest wall
        rough = 0.5d0 * (rough1 + rough2 ) ! roughness nearest wall
      end if

      ! Gradient in tangential directions: 3-10
      call Tangential_Gradient ( Visc3, Visc4, Visc5, Visc6,  Gradient(:,b) )
      call Tangential_Gradient ( Visc7, Visc8, Visc9, Visc10, Gradient(:,c) )

      !> Transform the gradients from computational to physical space.
      !> Explicit loop, not `Gradient = matmul(Gradient,M)`: Gradient aliases
      !> itself there, so the compiler must build the whole (nprim,3) result in a
      !> temporary on every face.
      if ( present(Mg) ) then
        do v = 1, nprim
          g = 1                                        ! the pressure slot: no concentration gradient
          if ( v == nh ) g = 2                         ! temperature
          if ( v >= nu .and. v <= nw ) g = 3           ! velocity
          if ( nrans > 0 .and. v >= nt ) g = 4         ! turbulence variables
          g1 = Gradient(v,1) ; g2 = Gradient(v,2) ; g3 = Gradient(v,3)
          Gradient(v,1) = g1*Mg(1,1,g) + g2*Mg(2,1,g) + g3*Mg(3,1,g)
          Gradient(v,2) = g1*Mg(1,2,g) + g2*Mg(2,2,g) + g3*Mg(3,2,g)
          Gradient(v,3) = g1*Mg(1,3,g) + g2*Mg(2,3,g) + g3*Mg(3,3,g)
        end do
      else
        do v = 1, nprim
          g1 = Gradient(v,1) ; g2 = Gradient(v,2) ; g3 = Gradient(v,3)
          Gradient(v,1) = g1*Mf(1,1) + g2*Mf(2,1) + g3*Mf(3,1)
          Gradient(v,2) = g1*Mf(1,2) + g2*Mf(2,2) + g3*Mf(3,2)
          Gradient(v,3) = g1*Mf(1,3) + g2*Mf(2,3) + g3*Mf(3,3)
        end do
      end if

      Prim = 0.5d0 * ( Prim1 + Prim2 )
      call Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, rough, Flux, Sc, Sct, Prt, Prt_corr )
      
      Res1 = Res1 - Flux
      Res2 = Res2 + Flux

    end subroutine Diffusive_Flux


    !> Metric tensor at the face between two cells of metrics M1 and M2:
    !> inverse of the mean of the inverse metrics (edge vectors dx/dxi).
    pure function Face_Metric ( M1, M2 ) result ( M )
      implicit none
      real(R8), intent(in) :: M1(3,3), M2(3,3)
      real(R8)             :: M(3,3)

      M = Inverse3 ( 0.5d0 * ( Inverse3 ( M1 ) + Inverse3 ( M2 ) ) )

    end function Face_Metric


    !> Inverse of a 3x3 matrix by its adjugate (transposed cofactors). Works for M
    !> as stored by Compute_Metric_Tensor (cofactors/det = transposed inverse of the
    !> edge vectors): transposition commutes with averaging and inversion.
    pure function Inverse3 ( A ) result ( B )
      implicit none
      real(R8), intent(in) :: A(3,3)
      real(R8)             :: B(3,3)
      real(R8) :: det

      B(1,1) =  ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
      B(1,2) = -( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
      B(1,3) =  ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
      B(2,1) = -( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
      B(2,2) =  ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
      B(2,3) = -( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
      B(3,1) =  ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
      B(3,2) = -( A(1,1)*A(3,2) - A(1,2)*A(3,1) )
      B(3,3) =  ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )
      det = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
      B = B / det

    end function Inverse3


    !> Mean of the central differences of two cells along a tangential direction, from viscous variables: Visc1 and Visc2
    !> the neighbours of the first cell (at -1 and +1), Visc3 and Visc4 those of the second
    subroutine Tangential_Gradient ( Visc1, Visc2, Visc3, Visc4, Gradient )
      
      use ARES_Global_m

      implicit none
      real(R8), intent(in), dimension(nprim) :: Visc1, Visc2, Visc3, Visc4
      real(R8), intent(out), dimension(nprim) :: Gradient

      Gradient = ( Visc2 - Visc1 + Visc4 - Visc3 ) * 0.25d0

    end subroutine Tangential_Gradient


    !> Viscous variables of a cell of primitive variables P: velocity, temperature and turbulence variables per unit mass
    !> (no concentration with a single species: zero in the pressure slot), the quantities whose differences make the
    !> viscous gradients. Viscous_Variables computes them once per stage for every cell (blk%Visc); the faces between
    !> blocks for the cells of the other block (bc%Pg).
    subroutine Visc_State ( P, Visc )
      use ARES_Global_m
      use FLINT_Lib_Thermodynamic, only: ph2vars, rho_tab, T_tab
      implicit none
      real(R8), intent(in)  :: P(nprim)
      real(R8), intent(out) :: Visc(nprim)
      ! Local
      real(R8) :: rho

      Visc(np) = 0.0_R8
      Visc(nu:nw) = P(nu:nw)
      Visc(nh) = ph2vars( P(np), P(nh), T_tab )
      if (nprim>nh) then
        rho = ph2vars( P(np), P(nh), rho_tab )
        Visc(nt:nprim) = P(nt:nprim) / rho
      endif

    end subroutine Visc_State


    subroutine Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, rough, Flux, Sc, Sct, Prt, Prt_corr )
    use ARES_Global_m
    use ARES_Lib_Fluid
    use ARES_Lib_RANS
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Prt_Correction, only: delta_Prt

      implicit none
      real(R8), intent(in) :: Prim(nprim), Gradient(nprim,3), area, normal(3), waldis, rough, Sc, Sct, Prt
      real(R8), intent(out) :: Flux(nprim)
      logical, intent(in) :: Prt_corr
      ! Local
      !integer :: s
      real(R8) :: rho, T, cp, mil, kl, mie, kappa
      real(R8) :: stress(3), DiffHFlux !, DmGradYi ,Dm(np) 
      real(R8) :: dPrt

      ! Thermodynamic and transport properties at the interface
      T   = ph2vars( Prim(np), Prim(nh), T_tab )
      rho = ph2vars( Prim(np), Prim(nh), rho_tab )
      cp  = ph2vars( Prim(np), Prim(nh), hT_tab )

      mil = ph2vars( Prim(np), Prim(nh), mi_tab2D )
      kl  = ph2vars( Prim(np), Prim(nh), k_tab2D )

      ! Eddy viscosity
      mie = 0d0
      if (model==2) then
        call Eddy_Viscosity ( mut=mie, rans_variables=Prim(nt:nprim), &
                              mul=mil, rho=rho, vel_gradient=Gradient(nu:nw,:), &
                              walldist=waldis, ks=rough )
      end if

      if ( Prt_Corr ) then
        dPrt = delta_Prt ( mil, cp, kl, Prim(nt), rough, waldis )
      else
        dPrt = 0.00d0
      end if
      kappa = kl + mie*cp/(Prt+dPrt) ! Laminar + turbulent conductivity
      Stress = stress_vector ( Gradient(nu:nw,:), normal, mil, mie, prim(nt:) ) ! Stress tensor in cartesian components

      ! Fluxes computation
      DiffHFlux = 0.d0
      Flux(np) = 0.0_R8 ! Non ho gradiente di concentrazione di specie -> una specie soltanto 

      ! Momentum flux
      Flux(nu:nw) = area * Stress

      ! Diffusive enthalpy flux
      Flux(nh) = area * ( dot_product ( Stress, Prim(nu:nw) ) + &
                 kappa * dot_product ( Gradient(nh,:), normal ) + DiffHFlux * rho )

      ! Turbulence variables diffusive flux
      if (model==2) then
        call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                   rans_variables=Prim(nt:nprim), &
                                   vel_gradient=Gradient(nu:nw,:), &
                                   rans_gradient=Gradient(nt:nprim,:), &
                                   mul=mil, rho=rho, &
                                   area=area, normal=normal, dist=waldis)
      end if

    end subroutine Compute_Diffusive_Flux

end module ARES_Lib_Diffusive