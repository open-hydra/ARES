module ARES_Lib_BC_Fluxes_Extrapolation
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index

  implicit none
  public

contains

  subroutine BC_Extrapolation ( Im, Jm, Km, Fm, Blk )
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: Blk
    ! Local
    integer :: Dir, Face_i, Face_j, Face_k, FluxSign
    real(R8) :: Normal(3), Area, Prim(nprim), Un, Fmass, Flux(nprim), rho


    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    
    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! boundary cell primitive variables
    Prim = Blk % P(:,Im,Jm,Km)
    rho = ph2vars( Prim(np), Prim(nh), rho_tab ) 
    ! boundary cell: normal velocity to the boundary face
    Un = dot_product ( Prim(nu:nw), Normal )

    Flux(np) = rho * Un * Area
    Fmass = rho * Un * Area
    Flux(nu:nw) =  Fmass * Prim(nu:nw) + Prim(np) * Area * Normal
    Flux(nh) = Fmass * ( Prim(nh) + 0.5_R8* norm2 ( Prim(nu:nw) )**2 )
    if (model==2) Flux(nt:nprim) = Prim(nt:nprim) * Un * Area
    
    ! Residual update
    FluxSign = 1 - 2 * mod(Fm,2)
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + FluxSign * Flux

  end subroutine BC_Extrapolation

end module ARES_Lib_BC_Fluxes_Extrapolation