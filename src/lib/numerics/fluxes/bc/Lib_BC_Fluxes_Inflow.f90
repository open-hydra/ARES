module ARES_Lib_BC_Fluxes_Inflow
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public

contains

  !─────────────────────────────────────────────────────────────────────────────
  ! BC 404: Inlet — static temperature T + prescribed mass flux g [kg/m2/s].
  ! Pressure is extrapolated from the interior. g is positive into the domain on every face
  ! (IO_BC.f90, Need_Inward).
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
    real(R8) :: Normal(3), Area, t_Vec(3), BC_Sign, Fmass, Face_Un, entalpy, g
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
    ! The mass flux into the domain, |g|, as its component along Normal, which points toward growing index:
    ! out of the domain on the faces 2, 4, 6 (BC_Sign = +1), into it on the faces 1, 3, 5 (BC_Sign = -1)
    g       = - BC_Sign * abs(BC_g)
    Un3     = g / rho3
    Ut3     = Un3 * E1
    Ub3     = Un3 * E2

    ! Flux computation (directly imposed state at interface)
    Face_Vel = Un3 * Normal + Ut3 * t_Vec + Ub3 * b_Vec
    Face_Un  = dot_product( Face_Vel, Normal )   ! true normal component (robust on warped faces)
    Fmass    = rho3 * Face_Un * Area
    Flux(np) = Fmass
    Flux(nu:nw) = Fmass * Face_Vel + p3 * Area * Normal
    Flux(nh) = Fmass * ( 0.5d0 * sum( Face_Vel**2 ) + entalpy )
    ! specific inputs (nuTilde, k, omega, Rij of ATLAS BCB) times the mass flux: the transported
    ! variables are density-weighted (rhoNuTilde, rhoK, rhoOmega, rho*Rij)
    if (model==2) Flux(nt:nprim) = BC_RANS * Fmass

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux


  end subroutine BC_Inlet_MassFlux_T


  !─────────────────────────────────────────────────────────────────────────────
  ! BC 405: Supersonic inlet — Mach, static temperature T and static pressure p. The flux of the face is that of the
  ! Riemann problem of the imposed state (Supersonic_Inflow: the static pressure and the Mach number blended with those
  ! of the cell by rel_fac, the direction from alpha and beta or along the normal of the face into the domain) with the
  ! boundary cell, through the 3-cell stencil.
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

    ! Supersonic_Inflow takes BC_p and BC_T as static values; the blending with the interior and the Riemann stencil are inside.
    Flux = 0d0
    call Supersonic_Inflow ( Bound_Prim, Int_Prim, modfm, modfm1, modfm2, BC_Mach, BC_p, &
                                 BC_T, BC_RANS, BC_rel_fac, BC_alpha, BC_beta,         &
                                 Normal, Area, Flux, BC_mdot )

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + modfm2 * Flux

  end subroutine BC_Inlet_Supersonic_Static


  !─────────────────────────────────────────────────────────────────────────────
  ! Supersonic_Inflow (low-level, works on primitive arrays directly). The imposed state: the static T (BC_T), the
  ! static p (BC_p) and the Mach number M, both blended with those of the cell by rel_fac, and the velocity M a d. With
  ! alpha and beta both normal, d = n = -modfm2 Normal, the unit normal of the face into the domain (Normal points
  ! toward growing index). Otherwise d = (cos alpha cos beta, sin alpha cos beta, sin beta), where an angle given as
  ! normal is that of n: alpha its azimuth, (cos, sin) = (n_x, n_y)/|(n_x, n_y)|, 0 when n is along z; beta its
  ! elevation, (cos, sin) = (|(n_x, n_y)|, n_z)/|n|.
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
    real(R8) :: n_in(3), h_xy, n_mod, cos_a, sin_a, cos_b, sin_b

    ! compute flow properties in the interior cell (sound from the boundary-cell state)
    flow_Sound = ph2vars( Prim2(np), Prim2(nh), sound_tab )
    flow_Mach  = norm2( Prim2(nu:nw) ) / flow_Sound

    ! blended Mach number at the boundary
    BC_Mach_local = BC_rel_fac*BC_Mach + (1-BC_rel_fac)*flow_Mach
    ! blended p static at the boundary
    BC_p0_local = BC_rel_fac*BC_p + (1-BC_rel_fac)*Prim2(np)

    ! supersonic inflow. BC_Mach enforced at boundary, T0 and p0 are static
    Sup_p = BC_p0_local
    Sup_T = BC_T
    Sup_h = pT2h( Sup_p, Sup_T )
    Sup_rho = ph2vars( Sup_p, Sup_h, rho_tab )
    Sup_Prim(np) = Sup_p
    Sup_Prim(nh) = Sup_h
    Sup_Sound = ph2vars( Sup_p, Sup_h, sound_tab )
    n_in = - modfm2 * Normal
    if ( Is_Normal_Dir ( BC_alpha ) .and. Is_Normal_Dir ( BC_beta ) ) then
      Sup_Prim(nu:nw) = BC_Mach_local * Sup_Sound * n_in
    else
      if ( Is_Normal_Dir ( BC_alpha ) ) then
        h_xy = norm2 ( n_in(1:2) )
        cos_a = 1d0
        sin_a = 0d0
        if ( h_xy > 0d0 ) then
          cos_a = n_in(1) / h_xy
          sin_a = n_in(2) / h_xy
        end if
      else
        cos_a = cos ( BC_alpha )
        sin_a = sin ( BC_alpha )
      end if
      if ( Is_Normal_Dir ( BC_beta ) ) then
        n_mod = norm2 ( n_in )
        cos_b = 1d0
        sin_b = 0d0
        if ( n_mod > 0d0 ) then
          cos_b = norm2 ( n_in(1:2) ) / n_mod
          sin_b = n_in(3) / n_mod
        end if
      else
        cos_b = cos ( BC_beta )
        sin_b = sin ( BC_beta )
      end if
      Sup_Prim(nu) = BC_Mach_local * Sup_Sound * cos_a * cos_b
      Sup_Prim(nv) = BC_Mach_local * Sup_Sound * sin_a * cos_b
      Sup_Prim(nw) = BC_Mach_local * Sup_Sound * sin_b
    end if
    ! specific inputs (nuTilde, k, omega, Rij of ATLAS BCB) times the inlet density
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
    if (Is_Normal_Dir(alpha_in)) then
      alpha_out = atan(Normal(2) / (Normal(1) + 1d-20))
    else
      alpha_out = alpha_in
    end if
    if (Is_Normal_Dir(beta_in)) then
      beta_out  = atan(Normal(3) / (Normal(1) + 1d-20))
    else
      beta_out  = beta_in
    end if

  end subroutine Compute_Inflow_Direction


  !─────────────────────────────────────────────────────────────────────────────
  ! Is_Normal_Dir: the flow angle x of an inlet record is the token normal, which parse_dir_tok (IO_BC.f90) reads as
  ! huge(1.0_R8).
  pure logical function Is_Normal_Dir ( x )
    implicit none
    real(R8), intent(in) :: x

    Is_Normal_Dir = x >= 0.5_R8 * huge(1.0_R8)

  end function Is_Normal_Dir


end module ARES_Lib_BC_Fluxes_Inflow