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
  !> conductivity at the state (p, h) of the boundary cell; the wall temperature Twall = temp + Heat_Flux / (2 (kl+kt)
  !> Mw(Dir,:).n), from the heat flux across the half cell with the laminar and turbulent conductivities (exact when
  !> adiabatic, and the energy flux is Heat_Flux regardless), Mw the metric of the wall gradient (with diffusive-metric =
  !> stencil that of the wall face, rec%wM of Wall_Face_Stencils; otherwise Mc, that of the cell), the same as the wall
  !> flux of BC_Wall_Heat; the enthalpy and the density at
  !> Twall and the pressure of the cell (boundary layer). kt = mit cp / (Prt + dPrt) at the wall face, as archive-ARES
  !> took it (explicit, at the state of the cell): mit the eddy viscosity at the wall face (y = 0) from the wall values of
  !> the turbulence variables of the cell values rans_cell (rho nu~ ...) at the cell density rho, cp of the cell, dPrt
  !> the rough-wall correction at y = 0. It is zero but with SA-rough on a wall with ks > 0: k = 0, R_ij = 0 and nu~ = 0
  !> (smooth SA) at the wall give mit = 0, and then Twall is that of kl alone. The eddy viscosity at the
  !> wall face takes the velocity gradient of the wall face, the velocity u of the cell over the half cell (no-slip,
  !> u = 0 at the wall: the index-space difference u modfm3 along Dir of the face Fm, turned with Mw; the velocity is
  !> zero along the wall, so no tangential difference), as BC_Wall_Heat takes it: it enters the limiters of SST and
  !> Wilcox2006 (moot while k = 0 at the wall), not SA and SSG/LRR. The wall values of the turbulence variables take the
  !> metric of the cell Mc and the wall distance of the cell, as in every caller of RANS_Set_Wall_Values. BC_Wall_Heat
  !> and Ghost_Wall pass the same Mw (from the record), so that they give the same Twall.
  !> Twall is an estimate (first order; coarse meshes, coarse multigrid levels): when it lies outside the temperatures of
  !> the tables at that pressure it is brought within them (Table_Bounds, as the ghost cells of Ghost_Wall), and
  !> counted; on tables with no temperature common to the two rows that bracket p, the wall takes the state of the cell.
  subroutine Wall_State_Heat ( p, h, u, Heat_Flux, Mc, Mw, normal, Fm, Dir, rans_cell, rho, yn, ks, temp, mil, kl, &
                               Twall, h_wall, rho_wall )
    use FLINT_Lib_Thermodynamic,   only: ph2vars, pT2h, rho_tab, T_tab, mi_tab2D, k_tab2D, hT_tab
    use ARES_Load_ThermoTransport, only: Table_Bounds
    use ARES_Global_m,             only: model, nrans
    use ARES_Config_Types_m,       only: obj_rans
    use ARES_Lib_RANS,             only: RANS_Set_Wall_Values, Eddy_Viscosity
    use ARES_Lib_Prt_Correction,   only: delta_Prt
    implicit none
    real(R8), intent(in)  :: p, h, u(3), Heat_Flux, normal(3)   ! u: velocity of the cell
    real(R8), intent(in)  :: Mc(3,3), Mw(3,3)                   ! metric of the cell, of the wall gradient
    integer,  intent(in)  :: Fm, Dir
    real(R8), intent(in)  :: rans_cell(:)   ! turbulence variables of the cell (rho nu~, ...), read with model == 2 only
    real(R8), intent(in)  :: rho, yn, ks    ! cell density, wall distance of the cell, roughness (0 on a smooth wall)
    real(R8), intent(out) :: temp, mil, kl, Twall, h_wall, rho_wall
    ! Local
    real(R8) :: pw, metric, metric_c, kt, mit, cp, dPrt, rans_wall(nrans), gradu(3,3)
    integer  :: modfm, modfm1, modfm2, modfm3
    logical  :: inside, ok

    temp = ph2vars ( p, h, T_tab )
    mil  = ph2vars ( p, h, mi_tab2D )
    kl   = ph2vars ( p, h, k_tab2D )
    metric   = dot_product ( Mw(Dir,:), normal )   ! wall gradient: Twall
    metric_c = dot_product ( Mc(Dir,:), normal )   ! wall values of the turbulence variables

    ! Turbulent conductivity at the wall face (only where a heat flux crosses it: adiabatic, Twall = temp anyway)
    kt = 0d0
    if ( model == 2 .and. Heat_Flux /= 0d0 ) then
      call RANS_Set_Wall_Values ( rans_cell=rans_cell(1:nrans), mil_wall=mil, rans_wall=rans_wall, metric=metric_c, &
                                  dist=yn, ks=ks )
      ! Velocity gradient at the wall face (index space along Dir, then physical): as BC_Wall_Heat
      call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
      gradu = 0d0
      gradu(:,Dir) = u * modfm3
      gradu = matmul ( gradu, Mw )
      call Eddy_Viscosity ( mut=mit, rans_variables=rans_wall, mul=mil, rho=rho, vel_gradient=gradu, &
                            walldist=0d0, ks=ks )   ! wall face: y = 0 (yn stays the cell distance for Set_Wall_Values)
      if ( mit /= 0d0 ) then
        cp = ph2vars ( p, h, hT_tab )
        dPrt = 0d0
        if ( obj_rans%Prt_correction ) dPrt = delta_Prt ( mil, cp, kl, rans_wall(1), ks, 0d0 )
        kt = mit * cp / ( obj_rans%Prt + dPrt )
      end if
    end if

    ! Heat_Flux is positive into the fluid on every face: the wall is hotter than the cell when it heats the fluid
    Twall = temp + Heat_Flux / ( 2d0 * ( kl + kt ) * metric )
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