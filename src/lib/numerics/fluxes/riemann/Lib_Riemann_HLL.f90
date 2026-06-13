! OVERFLOW2.1 Riemann solvers
! Reference:
!   Tramel, R., Nichols, R., & Buning, P. (2009). Addition of improved shock-capturing schemes to OVERFLOW 2.1. In 19th AIAA Computational Fluid Dynamics (p. 3988).
! Note:
! - HLLC+  | The minimum allowed beta is set to 0.0 instead of 0.4 to allow pure HLLE flows
! - HLLE++ | The Roe lambda_1 is set to un_roe instead of max(un_roe,a_roe). 
!            The adopted forumlation preserves pure Eulerian shear flow, contrary to the reference one.

module ARES_Lib_Riemann_HLL
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: riemann_HLLE, riemann_HLLC
  public :: riemann_HLLCprec
  public :: riemann_HLLCHLLE

contains

  ! HLLE
  subroutine riemann_HLLE(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E )
    use FLINT_Lib_Thermodynamic

    implicit none
    integer s
    real(kind=8), intent(in)  :: pl,ul,vl,wl,hl ! : density(s), velocity, pressure and sound velocity of left state
    real(kind=8), intent(in)  :: pr,ur,vr,wr,hr ! : density(s), velocity, pressure and sound velocity of right state
    real(kind=8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    real(kind=8), intent(in)  :: nx, ny, nz
    ! common
    real(kind=8) :: rhol, rhor, al, ar, dl, dr, E0r, E0l
    ! specific
    real(kind=8) :: S1,S4,Frl,Frr,ful,fur,fvl,fvr,fwl,fwr,fel,fer
    real(kind=8) :: unl,unr
    real(kind=8) :: drtot_ROE,dltot_ROE, u_ROE,v_ROE,w_ROE,dtot_ROEsum, a_ROE, un_ROE
    
    ! stati L e R
    rhol = ph2vars ( pl, hl, rho_tab )
    rhor = ph2vars ( pr, hr, rho_tab )

    al = ph2vars ( pl, hl, sound_tab )
    ar = ph2vars ( pr, hr, sound_tab )

    dl = rhol
    dr = rhor

    unl=ul*nx+vl*ny+wl*nz
    unr=ur*nx+vr*ny+wr*nz

    !-------------------------------------------
    ! medie di Roe

    drtot_ROE=dsqrt(rhor)
    dltot_ROE=dsqrt(rhol)
    dtot_ROEsum=drtot_ROE+dltot_ROE

    u_ROE=(drtot_ROE*ur+dltot_ROE*ul)/dtot_ROEsum
    v_ROE=(drtot_ROE*vr+dltot_ROE*vl)/dtot_ROEsum
    w_ROE=(drtot_ROE*wr+dltot_ROE*wl)/dtot_ROEsum
    a_ROE=(drtot_ROE*ar+dltot_ROE*al)/dtot_ROEsum
    un_ROE=u_ROE*nx+v_ROE*ny+w_ROE*nz

    S1 = min(0.d0,unl-al,un_ROE-a_ROE)
    S4 = max(0.d0,unr+ar,un_ROE+a_ROE)

    ! computing fluxes
    select case(minloc([-S1,S1*S4,S4],dim=1))
    case(1)
      call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hl)
    case(2)
      call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,Frl,Ful,Fvl,Fwl,Fel,hl)
      call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,Frr,Fur,Fvr,Fwr,Fer,hr)
      E0l = hl - pl/rhol + 0.5d0*(ul**2+vl**2+wl**2)
      E0r = hr - pr/rhor + 0.5d0*(ur**2+vr**2+wr**2)
      F_r = (S4*Frl-S1*Frr+S1*S4*(rhor-rhol))/(S4-S1)
      F_u = (S4*Ful-S1*Fur+S1*S4*(rhor*ur-rhol*ul))/(S4-S1)
      F_v = (S4*Fvl-S1*Fvr+S1*S4*(rhor*vr-rhol*vl))/(S4-S1)
      F_w = (S4*Fwl-S1*Fwr+S1*S4*(rhor*wr-rhol*wl))/(S4-S1)
      F_E = (S4*Fel-S1*Fer+S1*S4*(rhor*E0r-rhol*E0l))/ &
            (S4-S1)
    case(3)
      call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hr)
    endselect

  end subroutine riemann_HLLE


  ! HLLC Batten
  subroutine riemann_HLLC(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E )
    use FLINT_Lib_Thermodynamic

    implicit none
    integer :: s
    real(kind=8), intent(in)  :: hl,ul,vl,wl,pl ! : density(s), velocity, pressure and sound velocity of left state
    real(kind=8), intent(in)  :: hr,ur,vr,wr,pr ! : density(s), velocity, pressure and sound velocity of right state
    real(kind=8), intent(out) :: F_r, F_u, F_v, F_w, F_e
    real(kind=8), intent(in)  :: nx, ny, nz
    ! common
    real(kind=8) :: al, ar
    ! specific
    real(kind=8) :: S1, S4, Sstar, pstar
    real(kind=8) :: E1, E4, U1S, U2S, U3S, U4S, U5S
    real(kind=8) :: unl,unr
    real(kind=8) :: drtot,dltot
    real(kind=8) :: drtot_ROE,dltot_ROE, a_ROE, un_ROE, &
                        u_ROE,v_ROE,w_ROE,dtot_ROEsum

    ! stati L e R
    dltot = ph2vars ( pl, hl, rho_tab )
    drtot = ph2vars ( pr, hr, rho_tab )

    al = ph2vars ( pl, hl, sound_tab )
    ar = ph2vars ( pr, hr, sound_tab )

    unl = ul*nx+vl*ny+wl*nz
    unr = ur*nx+vr*ny+wr*nz

    !-------------------------------------------
    ! medie di Roe
    drtot_ROE=dsqrt(drtot)
    dltot_ROE=dsqrt(dltot)
    dtot_ROEsum=drtot_ROE+dltot_ROE

    u_ROE=(drtot_ROE*ur+dltot_ROE*ul)/dtot_ROEsum
    v_ROE=(drtot_ROE*vr+dltot_ROE*vl)/dtot_ROEsum
    w_ROE=(drtot_ROE*wr+dltot_ROE*wl)/dtot_ROEsum
    a_ROE=(drtot_ROE*ar+dltot_ROE*al)/dtot_ROEsum
    un_ROE=u_ROE*nx+v_ROE*ny+w_ROE*nz

    S1 = min(0.d0,unl-al,un_ROE-a_ROE)
    S4 = max(0.d0,unr+ar,un_ROE+a_ROE)

    Sstar = (pr-pl+dltot*unl*(s1-unl)-drtot*unr*(s4-unr))/(dltot*(s1-unl)-drtot*(S4-unr))
    pstar = dltot*(unl-S1)*(unl-Sstar)+pl

    ! computing fluxes
    select case(minloc([-S1,S1*Sstar,Sstar*S4,S4],dim=1))
    case(1)
      call fluxes(pl,dltot,ul,vl,wl,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hl)
    case(2)
      call fluxes(pl,dltot,ul,vl,wl,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hl)
      !E1  = E0(pl,dl,sqrt(ul**2+vl**2+wl**2))
      E1 = hl - pl/dltot + 0.5d0*(ul**2+vl**2+wl**2)
      U1S = dltot*(S1-unl)/(S1-Sstar)
      U2S = ((S1-unl)*dltot*ul+(pstar-pl)*nx)/(S1-Sstar)
      U3S = ((S1-unl)*dltot*vl+(pstar-pl)*ny)/(S1-Sstar)
      U4S = ((S1-unl)*dltot*wl+(pstar-pl)*nz)/(S1-Sstar)
      U5S = ((S1-unl)*dltot*E1-pl*unl+pstar*Sstar)/(S1-Sstar)

      F_r = F_r + S1*(U1S - dltot)
      F_u = F_u + S1*(U2S - dltot*ul)
      F_v = F_v + S1*(U3S - dltot*vl)
      F_w = F_w + S1*(U4S - dltot*wl)
      F_E = F_E + S1*(U5S - dltot*E1)
    case(3)
      call fluxes(pr,drtot,ur,vr,wr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hr)
      !E4  = E0(pr,dr,sqrt(ur**2+vr**2+wr**2))
      E4 = hr - pr/drtot + 0.5d0*(ur**2+vr**2+wr**2)
      U1S = drtot*(S4-unr)/(S4-Sstar)
      U2S = ((S4-unr)*drtot*ur+(pstar-pr)*nx)/(S4-Sstar)
      U3S = ((S4-unr)*drtot*vr+(pstar-pr)*ny)/(S4-Sstar)
      U4S = ((S4-unr)*drtot*wr+(pstar-pr)*nz)/(S4-Sstar)
      U5S = ((S4-unr)*drtot*E4-pr*unr+pstar*Sstar)/(S4-Sstar)

      F_r = F_r + S4*(U1S - drtot)
      F_u = F_u + S4*(U2S - drtot*ur)
      F_v = F_v + S4*(U3S - drtot*vr)
      F_w = F_w + S4*(U4S - drtot*wr)
      F_E = F_E + S4*(U5S - drtot*E4)
    case(4)
      call fluxes(pr,drtot,ur,vr,wr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hr)
    endselect

  end subroutine riemann_HLLC


  ! HLLC Batten precondizionato secondo Luo–Baum–Löhner (2005)
  subroutine riemann_HLLCprec(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E )

    use ARES_Global_m, only: Uref
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Preconditioning

    implicit none

    real(kind=8), intent(in)  :: pl,ul,vl,wl,hl
    real(kind=8), intent(in)  :: pr,ur,vr,wr,hr
    real(kind=8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    real(kind=8), intent(in)  :: nx, ny, nz

    ! common
    real(kind=8) :: dltot,drtot,al,ar,dl,dr
    integer s
    ! specific HLLC
    real(kind=8) :: S1, S4, Sstar, pstar
    real(kind=8) :: E1, E4, U1S, U2S, U3S, U4S, U5S
    real(kind=8) :: unl,unr
    real(kind=8) :: drtot_ROE,dltot_ROE,u_ROE,v_ROE,w_ROE,dtot_ROEsum, un_ROE
    ! precondizionamento nei flussi (Luo–Baum–Löhner)
    real(kind=8) :: VrL,VrR
    real(kind=8) :: alphaL,alphaR
    real(kind=8) :: vL_star,vR_star,vROE_star
    real(kind=8) :: cL_star,cR_star,cROE_star
    real(kind=8) :: ModVelL,ModVelR
    real(kind=8) :: vnL,vnR
    real(kind=8) :: dx

    !-------------------------------------------
    ! Stati L e R
    !-------------------------------------------

    dltot = ph2vars ( pl, hl, rho_tab )
    drtot = ph2vars ( pr, hr, rho_tab )

    dl = dltot
    dr = drtot
    
    if ( Uref < 0.0d0 ) then
      al = ph2vars ( pl, hl, sound_tab )
      ar = ph2vars ( pr, hr, sound_tab )
    else
      al = Uref
      ar = Uref
    end if

    unl=ul*nx+vl*ny+wl*nz
    unr=ur*nx+vr*ny+wr*nz

    !-------------------------------------------
    ! Medie di Roe 
    !-------------------------------------------
    drtot_ROE   = dsqrt(drtot)
    dltot_ROE   = dsqrt(dltot)
    dtot_ROEsum = drtot_ROE + dltot_ROE

    !-------------------------------------------
    ! Preconditioning HLLC (Luo–Baum–Löhner, eq. 14–19, 27–28)
    !-------------------------------------------

    ! Velocità locali
    ModVelL  = dsqrt( ul**2 + vl**2 + wl**2 )
    ModVelR  = dsqrt( ur**2 + vr**2 + wr**2 )

    vnL   = unl
    vnR   = unr

    ! Reference velocity Vr 
    VrL   = comp_Ur( vel=ModVelL,  sound=al )
    VrR   = comp_Ur( vel=ModVelR,  sound=ar )

    ! Beta = rho_p + rho_T/(rho Cp)
    alphaL    = 0.5d0*( 1.d0 - 1.d0 / al**2 * VrL**2  )
    alphaR    = 0.5d0*( 1.d0 - 1.d0 / ar**2 * VrR**2  )

    ! v* e c* (preconditioned eigenvalues)
    vL_star    = vnL   * (1.d0 - alphaL)
    vR_star    = vnR   * (1.d0 - alphaR)
    vROE_star  = (drtot_ROE*vR_star + dltot_ROE*vL_star)/dtot_ROEsum

    cL_star    = dsqrt( alphaL  **2 * vnL  **2 + VrL  **2 )
    cR_star    = dsqrt( alphaR  **2 * vnR  **2 + VrR  **2 )
    cROE_star  = (drtot_ROE*cR_star + dltot_ROE*cL_star)/dtot_ROEsum

    ! Signal velocities precondizionate (eq. 27–28)
    S1 = min( 0.d0, vL_star   - cL_star,   vROE_star - cROE_star )
    S4 = max( 0.d0, vR_star   + cR_star,   vROE_star + cROE_star )

    !-------------------------------------------
    ! Stato star e pressione star (come Batten)
    !-------------------------------------------
    Sstar = ( pr - pl + dltot*unl*(S1-unl) - drtot*unr*(S4-unr) ) &
          /( dltot*(S1-unl)       - drtot*(S4-unr) )
    pstar = dltot*(unl - S1)*(unl - Sstar) + pl

    !-------------------------------------------
    ! Flussi HLLC con S1,S4 precondizionati
    !-------------------------------------------
    select case( minloc( [-S1, S1*Sstar, Sstar*S4, S4], dim=1 ) )

    case(1)

      call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hl)

    case(2)

      call fluxes(pl,dl,ul,vl,wl,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hl)
      E1 = hl - pl/dltot + 0.5d0*(ul**2+vl**2+wl**2)
      U1S = dltot*(S1-unl)/(S1-Sstar)
      U2S = ((S1-unl)*dltot*ul + (pstar-pl)*nx)/(S1-Sstar)
      U3S = ((S1-unl)*dltot*vl + (pstar-pl)*ny)/(S1-Sstar)
      U4S = ((S1-unl)*dltot*wl + (pstar-pl)*nz)/(S1-Sstar)
      U5S = ((S1-unl)*dltot*E1 - pl*unl + pstar*Sstar)/(S1-Sstar)

      F_r = F_r + S1*(U1S - dltot)
      F_u = F_u + S1*(U2S - dltot*ul)
      F_v = F_v + S1*(U3S - dltot*vl)
      F_w = F_w + S1*(U4S - dltot*wl)
      F_E = F_E + S1*(U5S - dltot*E1)
      
    case(3)

      call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hr)
      E4 = hr - pr/drtot + 0.5d0*(ur**2+vr**2+wr**2)
      U1S = drtot*(S4-unr)/(S4-Sstar)
      U2S = ((S4-unr)*drtot*ur + (pstar-pr)*nx)/(S4-Sstar)
      U3S = ((S4-unr)*drtot*vr + (pstar-pr)*ny)/(S4-Sstar)
      U4S = ((S4-unr)*drtot*wr + (pstar-pr)*nz)/(S4-Sstar)
      U5S = ((S4-unr)*drtot*E4 - pr*unr + pstar*Sstar)/(S4-Sstar)

      F_r = F_r + S4*(U1S - drtot)
      F_u = F_u + S4*(U2S - drtot*ur)
      F_v = F_v + S4*(U3S - drtot*vr)
      F_w = F_w + S4*(U4S - drtot*wr)
      F_E = F_E + S4*(U5S - drtot*E4)

    case(4)

      call fluxes(pr,dr,ur,vr,wr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E,hr)

    end select

  end subroutine riemann_HLLCprec


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


  subroutine riemann_HLLCHLLE(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E )
    use FLINT_Lib_Thermodynamic

    implicit none
    real(kind=8), intent(in)  :: hl,ul,vl,wl,pl ! : density(s), velocity, pressure and sound velocity of left state
    real(kind=8), intent(in)  :: hr,ur,vr,wr,pr ! : density(s), velocity, pressure and sound velocity of right state
    real(kind=8), intent(out) :: F_r, F_u, F_v, F_w, F_E
    real(kind=8), intent(in)  :: nx, ny, nz
    ! specific
    real(kind=8) :: nx1, ny1, nz1, nx2, ny2, nz2, dltot,drtot
    real(kind=8) :: alfa1, alfa2
    real(kind=8) :: F_rHLLE, F_uHLLE, F_vHLLE,F_wHLLE, F_EHLLE,F_rHLLC, F_uHLLC, F_vHLLC, F_wHLLC, F_EHLLC
    real(kind=8) :: abs_dv, dum
    

    ! calcolo del modulo del vettore differenza di velocità
    abs_dv=sqrt((ur-ul)**2+(vr-vl)**2+(wr-wl)**2)
    dum=0.5d0*(sqrt(ul**2+vl**2+wl**2)+sqrt(ur**2+vr**2+wr**2))
    
    if (abs_dv/(dum+1.0d-10)>1.0d-12) then        

    ! calcolo di n1: versore della differenza di velocità
    nx1=(ur-ul)/abs_dv
    ny1=(vr-vl)/abs_dv
    nz1=(wr-wl)/abs_dv

    ! calcolo di alfa1: proiezione di n su n1
    alfa1=nx*nx1+ny*ny1+nz*nz1

    ! rendo alfa1 sempre positivo checkando verso di n1
    nx1=sign(1.d0,alfa1)*nx1
    ny1=sign(1.d0,alfa1)*ny1
    nz1=sign(1.d0,alfa1)*nz1
    alfa1=sign(1.d0,alfa1)*alfa1

    nx2=-ny1
    ny2=nx1
    nz2=nz1
    ! calcolo di alfa2: proiezione di n su n2
    alfa2=nx*nx2+ny*ny2+nz*nz2

    ! rendo alfa2 sempre positivo checkando verso di n2
    nx2=sign(1.d0,alfa2)*nx2
    ny2=sign(1.d0,alfa2)*ny2
    nz2=sign(1.d0,alfa2)*nz2
    alfa2=sign(1.d0,alfa2)*alfa2
    
    if (nz>0.1d0) then
      alfa1=0.d0
      nx2=nx
      ny2=ny
      nz2=nz
      alfa2=1.d0
    endif
      
    ! chiamata a HLLE in direzione n1
    call riemann_HLLE(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx1,ny1,nz1,F_rHLLE,F_uHLLE,F_vHLLE,F_wHLLE,F_EHLLE )

    ! chiamata a HLLCb in direzione n2
    call riemann_HLLC(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx2,ny2,nz2,F_rHLLC,F_uHLLC,F_vHLLC,F_wHLLC,F_EHLLC )

    ! media pesata dei flussi per ottenere flusso rotato
    F_r=F_rHLLE*alfa1+F_rHLLC*alfa2
    F_u=F_uHLLE*alfa1+F_uHLLC*alfa2
    F_v=F_vHLLE*alfa1+F_vHLLC*alfa2
    F_w=F_wHLLE*alfa1+F_wHLLC*alfa2
    if (F_w /= 0.d0 .and. abs(nz) <= 1.d-10) F_w = 0.d0
    F_E=F_EHLLE*alfa1+F_EHLLC*alfa2
    else
    ! chiamata a HLLC in direzione n2, assumendo n1 tangente alla faccia
    call riemann_HLLC(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_rHLLC,F_uHLLC,F_vHLLC,F_wHLLC,F_EHLLC )
    
    F_r=F_rHLLC
    F_u=F_uHLLC
    F_v=F_vHLLC
    F_w=F_wHLLC
    F_E=F_EHLLC
    endif
  end subroutine riemann_HLLCHLLE

end module ARES_Lib_Riemann_HLL