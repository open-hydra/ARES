module ARES_Lib_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64

  implicit none
  public

  !> Evaluations of the state of a heat-flux wall (Wall_State_Heat), or of a wall coupled to the solid (103, 104, by the
  !> coupling of hydra), whose wall temperature, an estimate (from the heat flux across the half cell, or from the
  !> balance of the two sides), was brought within the fluid tables: reported at the end of the run with the ghost fills
  !> (Report_Wall_Ghosts)
  integer(I8) :: n_wall_heat_table = 0

contains

  pure subroutine Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    implicit none
    integer, intent(in) :: fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1-modFm  
    modfm2 = 1-2*modFm
    modfm3 = 2*( 2*modFm-1 )

  end subroutine Compute_Modfm


  pure subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
    implicit none
    integer, intent(in) :: face, cell_i, cell_j, cell_k
    integer, intent(out) :: dir, face_i, face_j, face_k

    select case ( face )
      case(1:2)
        dir = 1
        face_i = cell_i - mod ( face, 2 )
        face_j = cell_j
        face_k = cell_k
      case(3:4)
        dir = 2
        face_i = cell_i
        face_j = cell_j - mod ( face, 2 )
        face_k = cell_k
      case(5:6)
        dir = 3
        face_i = cell_i
        face_j = cell_j
        face_k = cell_k - mod ( face, 2 )
    end select

  end subroutine Face_Index


  !> Gas at the face of a wall of prescribed temperature, as BC_Wall_Temperature and Ghost_Wall take it: the pressure of
  !> the boundary cell (boundary layer), the enthalpy, the density and the laminar viscosity and conductivity at T_wall.
  subroutine Wall_State_Temperature ( p, T_wall, h_wall, rho_wall, mil, kl )
    use FLINT_Lib_Thermodynamic, only: ph2vars, pT2h, rho_tab, mi_tab2D, k_tab2D
    implicit none
    real(R8), intent(in)  :: p, T_wall
    real(R8), intent(out) :: h_wall, rho_wall, mil, kl

    h_wall   = pT2h ( p, T_wall )
    rho_wall = ph2vars ( p, h_wall, rho_tab )
    mil      = ph2vars ( p, h_wall, mi_tab2D )
    kl       = ph2vars ( p, h_wall, k_tab2D )

  end subroutine Wall_State_Temperature


  !> Gas at the face of a wall of prescribed heat flux, as BC_Wall_Heat and Ghost_Wall take it: the laminar viscosity and
  !> conductivity at the state (p, h) of the boundary cell; the wall temperature Twall = temp + Heat_Flux / (2 kl Mc(Dir,:).n),
  !> from the heat flux across the half cell with the laminar conductivity alone, also on a rough wall (exact when
  !> adiabatic, and the energy flux is Heat_Flux regardless), Mc the metric of the cell; the enthalpy and the density at
  !> Twall and the pressure of the cell (boundary layer). Twall is an estimate (first order; coarse meshes, coarse
  !> multigrid levels): when it lies outside the temperatures of the tables at that pressure it is brought within them
  !> (Table_Bounds, as the ghost cells of Ghost_Wall), and counted; on tables with no temperature common to the two rows
  !> that bracket p, the wall takes the state of the cell.
  subroutine Wall_State_Heat ( p, h, Heat_Flux, Mc, normal, Dir, temp, mil, kl, Twall, h_wall, rho_wall )
    use FLINT_Lib_Thermodynamic,   only: ph2vars, pT2h, rho_tab, T_tab, mi_tab2D, k_tab2D
    use ARES_Load_ThermoTransport, only: Table_Bounds
    implicit none
    real(R8), intent(in)  :: p, h, Heat_Flux, Mc(3,3), normal(3)
    integer,  intent(in)  :: Dir
    real(R8), intent(out) :: temp, mil, kl, Twall, h_wall, rho_wall
    ! Local
    real(R8) :: pw
    logical  :: inside, ok

    temp = ph2vars ( p, h, T_tab )
    mil  = ph2vars ( p, h, mi_tab2D )
    kl   = ph2vars ( p, h, k_tab2D )
    ! Heat_Flux is positive into the fluid on every face: the wall is hotter than the cell when it heats the fluid
    Twall = temp + Heat_Flux / ( 2d0 * kl * dot_product ( Mc(Dir,:), normal ) )
    pw = p
    call Table_Bounds ( pw, Twall, inside, ok )
    if ( .not. inside ) then
      !$omp atomic
      n_wall_heat_table = n_wall_heat_table + 1
    end if
    if ( ok ) then
      h_wall = pT2h ( pw, Twall )     ! pw = p within the tables (bit-identical there)
    else
      Twall  = temp
      h_wall = h
    end if
    rho_wall = ph2vars ( p, h_wall, rho_tab )

  end subroutine Wall_State_Heat


  !> The properties of a wall face for the wall output, the same on every face: the stress is the traction of the gas on
  !> the wall, tau . n with n the normal into the gas; qw is the heat flux into the gas (out of the solid). The callers
  !> give stress and qw along the face normal of the solver, toward growing index: stress = tau . n, qw = k grad T . n
  !> (the heat flux into the wall times into), where into = +1 when that normal points into the gas (faces 1, 3, 5) and
  !> -1 when it points into the wall (faces 2, 4, 6). A zero is written as +0 (Positive_Zero).
  subroutine Compute_Wall_Properties(stress,pw,into,Tw,rhow,mu,qw,y,exit_array)
    use iso_fortran_env, only: int64
    implicit none
    real(R8), dimension(:), intent(in)    :: stress
    real(R8), intent(in)                  :: pw
    integer, intent(in)                   :: into
    real(R8), dimension(8), intent(inout) :: exit_array
    real(R8), optional :: Tw, qw, rhow, mu, y
    ! Local
    real(R8) :: tauw, ut

    ! y+
    tauw = sqrt( sum ( Stress**2 ) )
    ut = sqrt ( tauw/rhow )
    exit_array(1) = rhow*y*ut/mu
    ! Wall stress: the traction of the gas on the wall
    exit_array(2:4) = Positive_Zero ( into * stress )
    ! Wall pressure
    exit_array(5) = pw
    ! Wall temperature
    if (present(Tw)) exit_array(6) = Tw
    ! Heat flux into the gas
    if (present(qw)) exit_array(7) = Positive_Zero ( - into * qw )

  contains

    !> x, with -0 written as +0. The wall output of a run on several ranks is gathered with a sum, which turns -0 into +0,
    !> so a -0 left here would make it differ, bit for bit, from that of a run on one rank. The test is on the bits: under a
    !> fast floating-point model (as -ffast-math allows) the compilers may drop an arithmetic one, such as x + 0.
    elemental function Positive_Zero ( x ) result ( z )
      real(R8), intent(in) :: x
      real(R8) :: z
      integer(int64) :: b
      b = transfer ( x, b )
      if ( ibclr ( b, 63 ) == 0_int64 ) b = 0_int64
      z = transfer ( b, z )
    end function Positive_Zero

  end subroutine Compute_Wall_Properties

end module ARES_Lib_BC_Fluxes