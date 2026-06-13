module ARES_Lib_Convective
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Convective_Flux

contains

    subroutine Convective_Flux ( dl, normal, area, Prim, Res )
      use ARES_Global_m
      use FLINT_Lib_Thermodynamic
      use ARES_Mod_Riemann
      use ARES_Lib_Reconstruction, only: Reconstruction

      implicit none
      real(R8), intent(in), dimension(-1:2) :: dl
      real(R8), intent(in) :: normal(3), area
      real(R8), intent(in), dimension(nprim,-1:2) :: Prim
      real(R8), intent(inout), dimension(nprim,0:1) :: Res
      ! Local
      real(R8) :: l0, l1, l2, lm, lp
      real(R8), dimension(nprim) :: Prim1, Prim4
      real(R8) :: rho1, rho4
      real(R8) :: F_r, F_u, F_v, F_w, F_E, su, sp, sm, Flux(nprim)


      l0 = 0.5d0 * ( dl(-1) + dl(0) )
      l1 = 0.5d0 * ( dl(0) + dl(1) )
      l2 = 0.5d0 * ( dl(1) + dl(2) )
      lm = 0.5d0 * dl(0)
      lp = 0.5d0 * dl(1)
      
      ! Reconstruction phase. Stencil around interface /i (cells /i and /i+1): (i-1),(i),(i+1),(i+2)
      call Reconstruction ( Prim(:,-1), Prim(:,0), Prim(:,1), Prim(:,2), l0, l1, l2, lm, lp, Prim1, Prim4 )

      ! Compute density 1&4
      rho1 = ph2vars( Prim1(np), Prim1(nh), rho_tab ) 
      rho4 = ph2vars( Prim4(np), Prim4(nh), rho_tab ) 

      ! Riemann solve
      call Riemann (prim1(np), prim1(nu), prim1(nv), prim1(nw), prim1(nh), &
                    prim4(np), prim4(nu), prim4(nv), prim4(nw), prim4(nh), &
                    normal(1),normal(2),normal(3), F_r, F_u, F_v, F_w, F_E )

      ! Sign of velocity at the interface
      su = sign ( 0.5d0, F_r )
      sp = 0.5d0 + su
      sm = su - 0.5d0

      ! Riemann fluxes are multiplied by the interface area
      Flux(1:np) = F_r * Area * ( sp - sm )
      Flux(nu)    = F_u * Area
      Flux(nv)    = F_v * Area
      Flux(nw)    = F_w * Area
      Flux(nh)    = F_E * Area

      if (nt>nh) then
        Flux(nt:nprim) = F_r * Area * ( sp * Prim1(nt:nprim)/rho1 - sm * Prim4(nt:nprim)/rho4 )
      end if

      ! Conservation-form residuals
      Res (:,0) = Res (:,0) + Flux
      Res (:,1) = Res (:,1) - Flux

    end subroutine Convective_Flux

end module ARES_Lib_Convective