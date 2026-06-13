module ARES_Lib_BC_Fluxes_Outflow
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public

contains

 subroutine BC_Outflow ( Bm, Im, Jm, Km, Fm, Blk, BC_pexit, BC_rel_fac, BC_mdot, error )
    implicit none
    type(ARES_block_type), intent(inout) :: Blk
    integer,  intent(in) :: Bm, Im, Jm, Km, Fm
    real(R8), intent(in) :: BC_rel_fac, BC_pexit
    real(R8), intent(out) :: BC_mdot
    integer :: error
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Int_i, Int_j, Int_k, Dir, Face_i, Face_j, Face_k
    real(R8) :: Bound_Prim(nprim), Int_Prim(nprim), Normal(3), Area, BC_Sign, Un
    real(R8) :: Flux(nprim)

    error = 0

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    ! (im,jm,km): boundary cell; (Int_i,Int_j,Int_k) next cell
    Int_i = Im + guide(Fm,1)
    Int_j = Jm + guide(Fm,2)
    Int_k = Km + guide(Fm,3)

    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Boundary cell primitive variables
    Bound_Prim = Blk % P(:,Im,Jm,Km)
    Int_Prim = Blk % P(:,Int_i,Int_j,Int_k)

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    BC_Sign = Real ( modfm2 )  ! =-1 for faces 1,3,5 ; =1 for faces 2,4,6
    Un  = dot_product ( Bound_Prim(nu:nw), Normal )  ! velocity normal to the interface

    if (BC_Sign * Un < 0d0) then
      error = 1
      return
    endif

    call Outflow( Bound_Prim, BC_Sign, Un, Normal, Area, BC_pexit, Flux, &
                  Int_Prim, modfm, modfm1, modfm2, BC_mdot )

    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + Modfm2 * Flux
  
  end subroutine BC_Outflow


  subroutine Outflow( Prim, BC_Sign, Un, Normal, Area, BC_pexit, Flux, &
                      Int_Prim, modfm, modfm1, modfm2, Fmass )
    use ARES_Lib_Limiters, only : rlimiter
    use FLINT_Lib_Thermodynamic
    use ARES_Mod_Riemann
    use ARES_Global_m
    implicit none
    real(R8), intent(in)  :: Prim(nprim), BC_Sign, Un, Normal(3), Area, BC_pexit
    real(R8), intent(in)  :: Int_Prim(nprim)
    integer,  intent(in)  :: modfm, modfm1, modfm2
    real(R8), intent(out) :: Flux(nprim), Fmass
    ! Local
    integer :: s
    real(R8) :: Mach_N
    real(R8), dimension(nprim) :: Diff32, Diff21, Slope, Face_Prim, Prim_L, Prim_R, Exit_prim, Prim1, Prim3
    real(R8) :: rho_L, rho_R, F_r, F_u, F_v, F_w, F_E
    real(R8) :: su, Sel_L, Sel_R, Bound_rho, Bound_Sound

    Bound_rho   = ph2vars( Prim(np), Prim(nh), rho_tab )
    Bound_Sound = ph2vars( Prim(np), Prim(nh), sound_tab )
    Mach_N = abs ( Un / Bound_Sound )   ! Mach number (NB normal to BC face)
  
    if ( ( Mach_N < 1d0 ) .and. ( BC_pexit > 0d0 ) ) then

      Exit_prim = Prim
      Exit_prim(np) = BC_pexit ! Impongo pressione dall'esterno

      ! stencil cell 1 is the ghost cell (odd faces) or the interior cell (even faces)
      Prim1 = modfm*Exit_prim + modfm1*Int_Prim
      ! stencil cell 3 is the interior cell (odd faces) or the ghost cell (even faces)
      Prim3 = modfm*Int_Prim + modfm1*Exit_prim

      ! The state at the interface from the interior side is reconstructed: (p,u,h)
      Diff32 = Prim3 - Prim
      Diff21 = Prim - Prim1
      do s = 1, nprim
        Slope(s) = rlimiter ( Diff32(s), Diff21(s) )
      end do
      
      ! Extrapolation at the interface from the boundary cell
      Face_Prim = Prim + 0.5d0 * Slope * modfm2

      ! state L of riemann problem is the ghost state for even faces and boundary reconstructed state for odd faces
      Prim_L = modfm*Exit_prim + modfm1*Face_Prim
      rho_L = ph2vars( Prim_L(np), Prim_L(nh), rho_tab )

      ! state R of riemann problem is the boundary reconstructed state for odd faces and the ghost state for even faces
      Prim_R = modfm*Face_Prim + modfm1*Exit_prim
      rho_R = ph2vars( Prim_R(np), Prim_R(nh), rho_tab )

      call Riemann ( Prim_L(np), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(nh), &
                      Prim_R(np), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(nh), &
                      Normal(1), Normal(2), Normal(3), F_r, F_u, F_v, F_w, F_E )

      su = sign ( 1d0, F_r )
      Sel_L = ( 1d0+su ) / 2d0  ! state left selector
      Sel_R = ( su-1d0 ) / 2d0  ! state right selector

      ! Fluxes
      Fmass = F_r * Area
      Flux(np) = Fmass * ( Sel_L - Sel_R )
      Flux(nu) = F_u * Area
      Flux(nv) = F_v * Area
      Flux(nw) = F_w * Area
      Flux(nh) = F_E * Area

      if (model==2) &
        Flux(nt:nprim) = Fmass * ( Sel_L*Prim_L(nt:nprim)/rho_L - Sel_R*Prim_R(nt:nprim)/rho_R )

    else
      ! supersonic exit or pambient=0: extrapolation of all variables
      Face_Prim = Prim ! Non impongo niente e prendo stato all'interfaccia
      Fmass = Bound_rho * Area * Un
      Flux(np) = Fmass
      Flux(nu:nw) = Face_Prim(nu:nw) * Fmass + Face_Prim(np) * Area * Normal
      Flux(nh) = Fmass * ( 0.5d0 * sum( Face_Prim(nu:nw)**2 ) + Face_Prim(nh) )
      if (model==2) Flux(nt:nprim) = Face_Prim(nt:nprim) / Bound_rho * Fmass
    
    end if

  end subroutine Outflow


end module ARES_Lib_BC_Fluxes_Outflow