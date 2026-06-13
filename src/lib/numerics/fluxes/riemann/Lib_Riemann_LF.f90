module ARES_Lib_Riemann_LF
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_PLLF, riemann_LLF

contains

  !> @brief Local Lax-Friedrichs (alias Rusanov)
subroutine Riemann_LLF(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
  use FLINT_Lib_Thermodynamic

  implicit none
  real(kind=8), intent(in)  :: hl,ul,vl,wl,pl ! : density(s), velocity, pressure and sound velocity of left state
  real(kind=8), intent(in)  :: hr,ur,vr,wr,pr ! : density(s), velocity, pressure and sound velocity of right state
  !real(kind=8), intent(in)  :: dltot,drtot
  real(kind=8), intent(out) :: F_r, F_u, F_v, F_w, F_e
  real(kind=8), intent(in)  :: nx, ny, nz
  ! Local
  real(8) :: uln, urn
  real(8) :: Frr, Fur, Fvr, Fwr, Fer, Frl, Ful, Fvl, Fwl, Fel
  real(8) :: A, al,ar, dl, dr

  uln=ul*nx+vl*ny+wl*nz
  urn=ur*nx+vr*ny+wr*nz

  dl = ph2vars ( pl, hl, rho_tab )
  dr = ph2vars ( pr, hr, rho_tab )
  ! Set spectral radius
  al = ph2vars ( pl, hl, sound_tab )
  ar = ph2vars ( pr, hr, sound_tab )
  A = MAX(ABS(uln-al),ABS(urn-ar),ABS(uln+al),ABS(urn+ar))

  call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,Frr,Fur,Fvr,Fwr,Fer,hr)
  call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,Frl,Ful,Fvl,Fwl,Fel,hl)
 
  F_r = 0.5*(Frl + Frr - A*(dr - dl))
  F_u = 0.5*(Fur + Ful - A*(dr*ur-dl*ul))
  F_v = 0.5*(Fvr + Fvl - A*(dr*vr-dl*vl))
  F_w = 0.5*(Fwr + Fwl - A*(dr*wr-dl*wl))
  F_E = 0.5*(Fel + Fer - A*(dr*(hr - pr/dr + 0.5d0*(ur**2+vr**2+wr**2)) - dl*(hl - pl/dl + 0.5d0*(ul**2+vl**2+wl**2)) ) )
     
end subroutine Riemann_LLF

! Preconditioned Local Lax-Friedrichs (alias Rusanov)
subroutine Riemann_PLLF (pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E)
  use FLINT_Lib_Thermodynamic

  implicit none
  real(kind=8), intent(in)  :: hl, ul, vl, wl, pl ! : density(s), velocity, pressure and sound velocity of left state
  real(kind=8), intent(in)  :: hr, ur, vr, wr, pr ! : density(s), velocity, pressure and sound velocity of right state
  real(kind=8), intent(out) :: F_r, F_u, F_v, F_w, F_e
  real(kind=8), intent(in)  :: nx, ny, nz
  ! common
  real(kind=8) :: dr,dl, ar, al
  ! specific
  real(kind=8) :: uln, urn, rho_Tl, rho_Tr, cpl, cpr, U_rl, U_rr
  real(kind=8) :: thetal, thetar, alphal, alphar, betal, betar, ul_mod, ur_mod
  real(kind=8) :: al_mod, ar_mod, A, Frr, Fur, Fvr, Fwr, Fer, Frl, Ful, Fvl, Fwl, Fel

  uln = ul*nx + vl*ny + wl*nz
  urn = ur*nx + vr*ny + wr*nz

  ! Thermodynamic stuff
  rho_Tl = ph2vars ( pl, hl, dT_tab )
  rho_Tr = ph2vars ( pr, hr, dT_tab )
  
  cpl = ph2vars ( pl, hl, hT_tab )
  cpr = ph2vars ( pr, hr, hT_tab )

  al = ph2vars ( pl, hl, sound_tab )
  ar = ph2vars ( pr, hr, sound_tab )

  dl = ph2vars ( pl, hl, rho_tab )
  dr = ph2vars ( pr, hr, rho_tab )

  ! Preconditioning stuff
  U_rl = Max ( 1d-5 * al, Min ( Abs( uln ) , al ) )
  U_rr = Max ( 1d-5 * ar, Min ( Abs( urn ) , ar ) )
  thetal = 1d0 / U_rl**2 - rho_Tl / ( dl * cpl )
  thetar = 1d0 / U_rr**2 - rho_Tr / ( dr * cpr )
  betal = 1d0 / al**2
  betar = 1d0 / ar**2
  alphal = 0.5d0 * ( 1d0 - betal * U_rl**2 )
  alphar = 0.5d0 * ( 1d0 - betar * U_rr**2 )
  ul_mod = Abs ( uln * ( 1d0 - alphal ) )
  ur_mod = Abs ( urn * ( 1d0 - alphar ) )
  al_mod = Sqrt ( alphal**2 * uln**2 + U_rl**2 )
  ar_mod = Sqrt ( alphar**2 * urn**2 + U_rr**2 )

  ! Set spectral radius
  A = Max ( ul_mod + al_mod, ur_mod + ar_mod )

  call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,Frr,Fur,Fvr,Fwr,Fer,hr)
  call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,Frl,Ful,Fvl,Fwl,Fel,hl)
 
  F_r = 0.5*(Frl + Frr - A*(dr - dl))
  F_u = 0.5*(Fur + Ful - A*(dr*ur-dl*ul))
  F_v = 0.5*(Fvr + Fvl - A*(dr*vr-dl*vl))
  F_w = 0.5*(Fwr + Fwl - A*(dr*wr-dl*wl))
  F_E = 0.5*(Fel + Fer -  A*(dr*(hr - pr/dr + 0.5d0*(ur**2+vr**2+wr**2)) - dl*(hl - pl/dl + 0.5d0*(ul**2+vl**2+wl**2)) ) )
     
end subroutine Riemann_PLLF

  ! Subroutine for computing the conservative fluxes from primitive variables.
  pure subroutine fluxes(p,r,u,v,w,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,h)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing the conservative fluxes from primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(kind=8), intent(IN)::  p           ! Pressure.
  real(kind=8), intent(IN)::  r           ! Density.
  real(kind=8), intent(IN)::  u,v,w       ! Velocity.
  real(kind=8), intent(IN)::  nx,ny,nz    ! Normals.
  real(kind=8), intent(OUT):: F_r         ! Flux of mass conservation.
  real(kind=8), intent(OUT):: F_u,F_v,F_w ! Flux of momentum conservation.
  real(kind=8), intent(OUT):: F_E         ! Flux of energy conservation.
  real(kind=8), intent(IN):: h            ! Entalpy.
  !---------------------------------------------------------------------------------------------------------------------------------
  F_r = r*(u*nx+v*ny+w*nz)
  F_u = F_r*u + p*nx
  F_v = F_r*v + p*ny
  F_w = F_r*w + p*nz
  F_E = F_r*(h+0.5d0*(u*u+v*v+w*w))
  !---------------------------------------------------------------------------------------------------------------------------------
end subroutine fluxes

end module ARES_Lib_Riemann_LF