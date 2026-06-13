module ARES_Lib_BC_Fluxes_Inflow
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public

contains

  subroutine BC_Inlet_MassFlux_T ( Bm, Im, Jm, Km, Fm, Blk, BC_T, BC_g, &
                                    BC_rel_fac, BC_alpha, BC_beta, BC_RANS, error )

    use FLINT_Lib_Thermodynamic
    use ARES_Global_m

    implicit none
    type(ARES_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in)    :: BC_T, BC_g, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_RANS(1:)
    integer,  intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign, Fmass, Face_Un, entalpy
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Sound, alpha, beta
    real(R8) :: E1, E2, rho3, p3, Un3, Ut3, Ub3
    real(R8) :: Face_Vel(3), Face_rho, Flux(nprim)
    real(R8) :: b_Vec(3), b_Mod, XA, XB, XC, XD, XE, XF, check

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Sound )

    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    call Compute_Inflow_Direction (Normal, BC_alpha, BC_beta, alpha, beta)

    ! Vector product n x t (why?)
    b_Vec(1) = Normal(2)*t_Vec(3) - Normal(3)*t_Vec(2)
    b_Vec(2) = Normal(3)*t_Vec(1) - Normal(1)*t_Vec(3)
    b_Vec(3) = Normal(1)*t_Vec(2) - Normal(2)*t_Vec(1)
    b_Mod = norm2 ( b_Vec )
    b_Vec = b_Vec/b_Mod

    ! Obscure esotheric stuff
    XA = Normal(2) - Normal(1) * tan(alpha)
    XB =  t_Vec(2) -  t_Vec(1) * tan(alpha)
    XC =  b_Vec(2) -  b_Vec(1) * tan(alpha)
    XD = Normal(3) - Normal(1) * tan(beta)
    XE =  t_Vec(3) -  t_Vec(1) * tan(beta)
    XF =  b_Vec(3) -  b_vec(1) * tan(beta)
    check = ( XC * XE - XF * XB )
    E1 = ( XF * XA - XC * XD ) / check ! ??
    E2 = ( XD * XB - XE * XA ) / check ! ??

    p3      = Bound_Prim(np)      ! pressure extrapolated from interior
    entalpy = pT2h( p3, BC_T )
    rho3    = ph2vars( p3, entalpy, rho_tab ) 
    Un3     = BC_g / rho3
    Ut3     = Un3 * E1
    Ub3     = Un3 * E2

    ! Flux computation (directly imposed state at interface)
    Face_Vel = Un3 * Normal + Ut3 * t_Vec + Ub3 * b_Vec
    Face_Un  = dot_product( Face_Vel, Normal )   ! true normal component (robust on warped faces)
    Fmass    = rho3 * Face_Un * Area
    Flux(np) = Fmass
    Flux(nu:nw) = Fmass * Face_Vel + p3 * Area * Normal
    Flux(nh) = Fmass * ( 0.5d0 * sum( Face_Vel**2 ) + entalpy )
    if (model==2) Flux(nt:nprim) = BC_RANS * Fmass

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux


  end subroutine BC_Inlet_MassFlux_T


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 406: Supersonic inlet — Mach + static temperature T and static pressure p.
  ! Note: BC_Supersonic_Inflow treats its BC_p0/BC_T0 arguments as static values,
  ! so the BC 406 static inputs map directly.  The full 3-cell Riemann stencil is
  ! preserved exactly as for BC 405.
  subroutine BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Km, Fm, Blk, BC_Mach, BC_T, BC_p, &
                                           BC_rel_fac, BC_alpha, BC_beta, BC_RANS, BC_mdot, error )
    implicit none
    type(ARES_block_type), intent(inout) :: Blk
    integer,  intent(in)    :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in)    :: BC_Mach, BC_T, BC_p, BC_rel_fac, BC_alpha, BC_beta
    real(R8), intent(in)    :: BC_RANS(1:)
    real(R8), intent(out)   :: BC_mdot
    integer,  intent(inout) :: error
    ! Local
    integer  :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8) :: Un, Bound_rho, Bound_Sound
    real(R8) :: Flux(nprim)

    error = 0

    call Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm, modfm, modfm1, modfm2, modfm3, &
                                  Int_i, Int_j, Int_k, Normal, Area, t_Vec, BC_Sign,   &
                                  Bound_Prim, Int_Prim, Un, Bound_rho, Bound_Sound )

    if (BC_Sign * Un >= 0d0) then
      error = 1
      return
    end if

    ! BC_Supersonic_Inflow treats BC_p0/BC_T0 as static p/T — pass through directly.
    ! Blending with the interior and the full Riemann stencil are handled inside.
    Flux = 0d0
    call Supersonic_Inflow ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p, &
                                 BC_T, BC_RANS, BC_rel_fac, BC_alpha, BC_beta,         &
                                 Normal, Area, Flux, BC_mdot )

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_Inlet_Supersonic_Static


  subroutine Supersonic_Inflow ( Prim2, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p, BC_T, & 
                                    BC_RANS, BC_rel_fac, BC_alpha, BC_beta, Normal, Area, Flux, Fmass )
    use ARES_Lib_Limiters, only : rlimiter
    use FLINT_Lib_Thermodynamic
    use ARES_Mod_Riemann
    use ARES_Global_m
    implicit none
    integer, intent(in) :: modfm, modfm1, modfm2
    real(R8), intent(in) :: Prim2(nprim), Int_Prim(nprim), BC_Mach, BC_p, BC_T, BC_rel_fac
    real(R8), intent(in) :: BC_RANS(1:), BC_alpha, BC_beta, Normal(3), Area
    real(R8), intent(out) :: Flux(nprim)
    real(R8), intent(out) :: Fmass
    ! Local
    integer :: s
    real(R8) :: Sup_Prim(nprim), Sup_T, Sup_Rgas, Sup_rho, Sup_Sound, Prim1(nprim), Prim3(nprim)
    real(R8), dimension(nprim) :: Diff32, Diff21, Slope, Face_Prim, Prim_L, Prim_R
    real(R8) :: rho_L, rho_R, F_r, F_u, F_v, F_w, F_E
    real(R8) :: su, Sel_L, Sel_R, Vel_T(3)
    real(R8) :: flow_Sound, flow_Mach, BC_Mach_local, BC_p0_local, Sup_p, Sup_h, entalpy, dx(3)
    real(R8) :: alpha_eff, beta_eff

    ! compute flow properties in the interior cell (sound from the boundary-cell state)
    flow_Sound = ph2vars( Prim2(np), Prim2(nh), sound_tab )
    flow_Mach  = norm2( Prim2(nu:nw) ) / flow_Sound

    ! blended Mach number at the boundary
    BC_Mach_local = BC_rel_fac*BC_Mach + (1-BC_rel_fac)*flow_Mach
    ! blended p static at the boundary
    BC_p0_local = BC_rel_fac*BC_p + (1-BC_rel_fac)*Prim2(np)

    ! resolve inflow direction (handles the 'normal' sentinel from parse_dir_tok)
    call Compute_Inflow_Direction ( Normal, BC_alpha, BC_beta, alpha_eff, beta_eff )

    ! supersonic inflow. BC_Mach enforced at boundary, T0 and p0 are static
    Sup_p = BC_p0_local
    Sup_T = BC_T
    Sup_h = pT2h( Sup_p, Sup_T )
    Sup_rho = ph2vars( Sup_p, Sup_h, rho_tab )
    Sup_Prim(np) = Sup_p
    Sup_Prim(nh) = Sup_h
    Sup_Sound = ph2vars( Sup_p, Sup_h, sound_tab )
    Sup_Prim(nu) = BC_Mach_local * Sup_Sound * cos(alpha_eff) * cos(beta_eff)
    Sup_Prim(nv) = BC_Mach_local * Sup_Sound * sin(alpha_eff) * cos(beta_eff)
    Sup_Prim(nw) = BC_Mach_local * Sup_Sound * sin(beta_eff)
    if (model==2) Sup_Prim(nt:nprim) = Sup_rho * BC_RANS

    ! stencil cell 1 is the ghost cell (odd faces) or the interior cell (even faces)
    Prim1 = modfm*Sup_Prim + modfm1*Int_Prim

    ! stencil cell 3 is the interior cell (odd faces) or the ghost cell (even faces)
    Prim3 = modfm*Int_Prim + modfm1*Sup_Prim
    
    ! The state at the interface from the interior side is reconstructed: (ro,u,p)
    Diff32 = Prim3 - Prim2
    Diff21 = Prim2 - Prim1
    do s = 1, nprim
      Slope(s) = rlimiter ( Diff32(s), Diff21(s) )
    end do
    
    ! Extrapolation at the interface from the boundary cell
    Face_Prim = Prim2 + 0.5d0 * Slope * modfm2

    ! state L of riemann problem is the ghost state for even faces and boundary reconstructed state for odd faces
    Prim_L = modfm*Sup_Prim + modfm1*Face_Prim
    rho_L = ph2vars( Prim_L(np), Prim_L(nh), rho_tab )

    ! state R of riemann problem is the boundary reconstructed state for odd faces and the ghost state for even faces
    Prim_R = modfm*Face_Prim + modfm1*Sup_Prim
    rho_R = ph2vars( Prim_R(np), Prim_R(nh), rho_tab )

    call Riemann ( Prim_L(np), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(nh), &
                    Prim_R(np), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(nh), &
                    Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_E )

    su = sign ( 1d0, F_r )
    Sel_L = ( 1d0+su ) / 2d0  ! state left selector
    Sel_R = ( su-1d0 ) / 2d0  ! state right selector

    ! Fluxes
    Fmass = F_r * Area
    Flux(1:np) = Fmass * ( Sel_L - Sel_R )
    Flux(nu) = F_u * Area
    Flux(nv) = F_v * Area
    Flux(nw) = F_w * Area
    Flux(nh) = F_E * Area

    if (model==2) &
      Flux(nt:nprim) = Fmass * ( Sel_L*Prim_L(nt:nprim)/rho_L - Sel_R*Prim_R(nt:nprim)/rho_R )

  end subroutine Supersonic_Inflow


  !─────────────────────────────────────────────────────────────────────────────
  ! Setup_Inflow_Geometry: compute all geometric/metric quantities common to
  ! every inflow BC.  Called at the start of each per-type routine.
  pure subroutine Setup_Inflow_Geometry ( Blk, Im, Jm, Km, Fm,               &
                                     modfm, modfm1, modfm2, modfm3,      &
                                     Int_i, Int_j, Int_k,                &
                                     Normal, Area, t_Vec, BC_Sign,       &
                                     Bound_Prim, Int_Prim,               &
                                     Un, Bound_rho, Bound_Sound )
    use FLINT_Lib_Thermodynamic
    implicit none
    type(ARES_block_type), intent(in)  :: Blk
    integer,  intent(in)               :: Im, Jm, Km, Fm
    integer,  intent(out)              :: modfm, modfm1, modfm2, modfm3
    integer,  intent(out)              :: Int_i, Int_j, Int_k
    real(R8), intent(out)              :: Normal(3), Area, t_Vec(3), BC_Sign
    real(R8), intent(out)              :: Bound_Prim(nprim), Int_Prim(nprim)
    real(R8), intent(out)              :: Un, Bound_rho, Bound_Sound
    ! Local
    integer  :: Dir, Face_i, Face_j, Face_k
    real(R8) :: t_Mod

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    ! (im,jm,km): boundary cell; (Int_i,Int_j,Int_k) next cell
    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)
    Int_k = Km + guide(Fm,3)

    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! boundary cell primitive variables
    Bound_Prim = Blk % P(:,Im,Jm,Km)
    Int_Prim = Blk % P(:,Int_i,Int_j,Int_k)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! Other obscure metric stuff 
    Select case ( Fm )
      Case ( 1 : 2 )
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i,Face_j-1,Face_k-1) % c
      Case ( 3 : 4 )
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i-1,Face_j,Face_k-1) % c
      Case ( 5 : 6 )
        t_Vec = Blk % node(Face_i,Face_j,Face_k) % c - &
                Blk % node(Face_i-1,Face_j-1,Face_k) % c
    end select

    t_Mod = norm2 ( t_Vec )
    t_Vec = t_Vec/t_Mod


    ! Compute stuff in the boundary cell
    Un  = dot_product ( Bound_Prim(nu:nw), Normal )  ! velocity normal to the interface
    Bound_rho = ph2vars( Bound_Prim(np), Bound_Prim(nh), rho_tab )  
    Bound_Sound = ph2vars( Bound_Prim(np), Bound_Prim(nh), sound_tab )  
    BC_Sign = real(modfm2)   ! -1 for faces 1,3,5 ; +1 for faces 2,4,6

  end subroutine Setup_Inflow_Geometry


  pure subroutine Compute_Inflow_Direction(normal, alpha_in, beta_in, alpha_out, beta_out)
    real(R8), intent(in) :: normal(3), alpha_in, beta_in
    real(R8), intent(out) :: alpha_out, beta_out
    real(R8) :: alpha_rad, beta_rad, cos_alpha, cos_beta

    ! Resolve 'normal' direction sentinel (parse_dir_tok returns huge(R8) for token 'normal').
    if (alpha_in >= 0.5_R8 * huge(1.0_R8)) then
      alpha_out = atan(Normal(2) / (Normal(1) + 1d-20))
    else
      alpha_out = alpha_in
    end if
    if (beta_in >= 0.5_R8 * huge(1.0_R8)) then
      beta_out  = atan(Normal(3) / (Normal(1) + 1d-20))
    else
      beta_out  = beta_in
    end if

  end subroutine Compute_Inflow_Direction


end module ARES_Lib_BC_Fluxes_Inflow