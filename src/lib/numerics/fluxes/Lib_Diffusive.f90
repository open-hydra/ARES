module ARES_Lib_Diffusive
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Diffusive_Flux, Compute_Diffusive_Flux

contains

    subroutine Diffusive_Flux ( normal, area, waldis1, waldis2, rough1, rough2, Prim1, Prim2, Prim3, Prim4, Prim5, &
                                Prim6, Prim7, Prim8, Prim9, Prim10, M1, M2, Res1, Res2, a, b, c, Sc, Sct, Prt, Prt_corr )
      
      use ARES_Global_m
      use FLINT_Lib_Thermodynamic

      implicit none
      integer, intent(in) :: a, b, c
      real(R8), intent(in) :: normal(3), area, waldis1, waldis2, rough1, rough2
      real(R8), intent(in), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4, Prim5, Prim6
      real(R8), intent(in), dimension(nprim) :: Prim7, Prim8, Prim9, Prim10
      real(R8), intent(in), dimension(3,3) :: M1, M2
      real(R8), intent(inout), dimension(nprim) :: Res1, Res2
      real(R8), intent(in) :: Sc, Sct, Prt
      logical , intent(in) :: Prt_corr
      ! Local
      real(R8) :: rho1, T1, rho2, T2, Gradient(nprim,3), Prim(nprim), M(3,3), waldis, Flux(nprim), rough

      waldis = 0.0_R8
      rough  = 0.0_R8

      ! Gradient in the same direction of the face: 1 and 2
      Gradient ( np, a ) = 0.0_R8 ! Solo una specie

      T1 = ph2vars( Prim1(np), Prim1(nh), T_tab )
      T2 = ph2vars( Prim2(np), Prim2(nh), T_tab )
      Gradient ( nh, a ) = T2 - T1 ! temperature gradient

      Gradient ( nu:nw, a ) = Prim2 ( nu:nw ) - Prim1 ( nu:nw ) ! velocity gradient
      if (model==2) then
        rho1 = ph2vars( Prim1(np), Prim1(nh), rho_tab )
        rho2 = ph2vars( Prim2(np), Prim2(nh), rho_tab )
        Gradient ( nt:nprim, a ) = Prim2 ( nt:nprim ) / rho2 - Prim1 ( nt:nprim ) / rho1 ! RANS variable gradient
        waldis = 0.5d0 * ( waldis1 + waldis2 ) ! distance to nearest wall
        rough = 0.5d0 * (rough1 + rough2 ) ! roughness nearest wall
      end if

      ! Gradient in tangential directions: 3-10

      call Tangential_Gradient ( Prim3, Prim4, Prim5, Prim6,  Gradient(:,b) )
      call Tangential_Gradient ( Prim7, Prim8, Prim9, Prim10, Gradient(:,c) )
      
      M = 0.5d0 * ( M1 + M2 )
      Gradient = matmul ( Gradient, M )

      Prim = 0.5d0 * ( Prim1 + Prim2 )
      call Compute_Diffusive_Flux ( Prim, Gradient, area, normal, waldis, rough, Flux, Sc, Sct, Prt, Prt_corr )
      
      Res1 = Res1 - Flux
      Res2 = Res2 + Flux

    end subroutine Diffusive_Flux


    subroutine Tangential_Gradient ( Prim1, Prim2, Prim3, Prim4, Gradient )
      
      use ARES_Global_m
      use FLINT_Lib_Thermodynamic

      implicit none
      real(R8), intent(in), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4
      real(R8), intent(out), dimension(nprim) :: Gradient
      ! Local
      real(R8) :: rho1, rho2, rho3, rho4, T1, T2, T3, T4

      Gradient(np) = 0.0_R8 ! Solo una specie -> no species concentration gradient

      T1 = ph2vars( Prim1(np), Prim1(nh), T_tab )
      T2 = ph2vars( Prim2(np), Prim2(nh), T_tab )
      T3 = ph2vars( Prim3(np), Prim3(nh), T_tab )
      T4 = ph2vars( Prim4(np), Prim4(nh), T_tab )

      Gradient ( nh ) = ( T2 - T1 + T4 - T3 ) * 0.25d0

      Gradient ( nu:nw ) = ( Prim2 ( nu:nw ) - Prim1 ( nu:nw ) + Prim4 ( nu:nw ) - Prim3 ( nu:nw ) ) * 0.25d0

      if (model==2) then
        rho1 = ph2vars( Prim1(np), Prim1(nh), rho_tab )
        rho2 = ph2vars( Prim2(np), Prim2(nh), rho_tab )
        rho3 = ph2vars( Prim3(np), Prim3(nh), rho_tab )
        rho4 = ph2vars( Prim4(np), Prim4(nh), rho_tab )
        Gradient ( nt:nprim ) = ( Prim2 ( nt:nprim ) / rho2 - Prim1 ( nt:nprim ) / rho1 + &
                                  Prim4 ( nt:nprim ) / rho4 - Prim3 ( nt:nprim ) / rho3 ) * 0.25d0
      end if

    end subroutine Tangential_Gradient


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