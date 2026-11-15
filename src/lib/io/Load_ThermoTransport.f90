! FLINT wrapper to load thermodynamic data
module ARES_Load_ThermoTransport
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_ThermoTransport, Outside_Tables, Table_Bounds, Table_Temperatures

contains

  subroutine Load_ThermoTransport()
    use ARES_Config_Types_m,        only: obj_thermo, obj_transport, obj_thermoinversion
    use ARES_Global_m,              only: ARES_phase_prefix
    use FLINT_Load_ThermoTransport!, only: read_idealgas_thermo, read_idealgas_transport
    use FLINT_Lib_Thermodynamic,    only: FLINT_phase_prefix
    use ARES_Parameters_m
    implicit none
    ! Local
    integer :: error

    !! ------------------------------------------------------
    !! Thermo-Transport Tables ------------------------------
    !! ------------------------------------------------------
    obj_thermo%warning_message = 'none'
    obj_thermo%error_message   = 'none'
    obj_thermo%description     = 'none'

    FLINT_phase_prefix = ARES_phase_prefix
    error = read_realfluid_thermo( 'INPUT' )

    select case (error)
      case (0)
        obj_thermo%description = 'Real fluid'
      case (5)
        obj_thermo%error_message = '[ERROR] Too many blocks in thermo file'
      case (1)
        obj_thermo%error_message = '[ERROR] Phase file (phase.txt) not found'
      case (2)
        obj_thermo%error_message = '[ERROR] Phase file (phase.txt) found but could not be read'
      case (3)
        obj_thermo%error_message = '[ERROR] Thermo table file not found'
      case (4)
        obj_thermo%error_message = '[ERROR] Thermo table file found but could not be read'
      case (6)
        obj_thermo%error_message = '[ERROR] Composition file found but could not be read'
      case default
        obj_thermo%error_message = '[ERROR] Unknown error loading thermodynamic data'
    end select

    obj_transport%warning_message = 'none'
    obj_transport%error_message   = 'none'
    obj_transport%description     = 'none'

    error = read_realfluid_transport( 'INPUT' )
    select case (error)
      case (0)
        obj_transport%description = 'Fluid transport properties'
      case (1)
        obj_transport%description   = 'Unavailable'
        obj_transport%error_message = '[ERROR] Transport table file not found'
      case (2)
        obj_transport%error_message = '[ERROR] Transport table file found but could not be read'
      case (3)
        obj_transport%error_message = '[ERROR] Too many blocks in transport file'
      case (4)
        obj_transport%error_message = '[ERROR] Transport data mesh size does not match thermo data'
      case default
        obj_transport%error_message = '[ERROR] Unknown error loading transport data'
    end select

    obj_thermoinversion%warning_message = 'none'
    obj_thermoinversion%error_message   = 'none'
    obj_thermoinversion%description     = 'none'
    error = ph2pT( )
    select case (error)
      case (0)
        obj_thermoinversion%description =   'pressure-entalpy table inverted succesfully'
      case (1)
        obj_thermoinversion%error_message = '[ERROR] pT2h: h_tab2D already allocated'
      case(2)
        obj_thermoinversion%error_message = '[ERROR] pT2h: Table dimensions too small (Nh<1 or Np<0)'
      case(3)
        obj_thermoinversion%error_message = '[ERROR] pT2h: T_tab not allocated'
      case(4)
        obj_thermoinversion%error_message = '[ERROR] pT2h: Interpolation failed'
      case default
        obj_thermoinversion%error_message = '[ERROR] pT2h: Unknown error'
    end select

  end subroutine Load_ThermoTransport


  !> True when ph2vars would index outside the tables at (p, h): more than one cell below their first pressure or enthalpy
  !> (within one cell below, ph2vars extrapolates from the first cell), or at or beyond the last. With strict, also the
  !> cell below (for values given by the user: initial and boundary conditions). A value that is not a number is outside.
  pure logical function Outside_Tables ( p, h, strict )
    use FLINT_Lib_Thermodynamic, only: pmin, deltap, hmin, deltah, rho_tab
    implicit none
    real(R8), intent(in)          :: p, h
    logical, intent(in), optional :: strict
    ! Local
    real(R8) :: xp, xh

    xp = ( p - pmin ) / deltap
    xh = ( h - hmin ) / deltah
    ! ph2vars reads the nodes i and i+1 (j and j+1): i = int(xp) must lie in 0 .. ubound - 1. The test is written so that
    ! a value that is not a number fails it
    Outside_Tables = .not. ( xp > -1d0 .and. xp < real ( ubound(rho_tab,1), R8 ) .and. &
                             xh > -1d0 .and. xh < real ( ubound(rho_tab,2), R8 ) )
    if ( present(strict) ) then
      if ( strict ) Outside_Tables = Outside_Tables .or. xp < 0d0 .or. xh < 0d0
    end if

  end function Outside_Tables


  !> The pressure and the temperature of a state brought within the fluid tables, where h(p, T) (pT2h) and rho(p, h)
  !> (ph2vars) read them without leaving their arrays: p within the pressures of the tables, T within the temperatures of
  !> both rows of the tables that bracket p, from h_min to h_max. inside: they were already within. ok: the two rows have
  !> temperatures in common (false only on tables whose T range jumps from one pressure to the next, where the caller
  !> takes the state of the cell).
  subroutine Table_Bounds ( p, T, inside, ok )
    use FLINT_Lib_Thermodynamic, only: pmin, deltap, Tmin2, deltaT, h_tab2D
    implicit none
    real(R8), intent(inout) :: p, T
    logical,  intent(out)   :: inside, ok
    ! Local
    integer  :: i, ip, ih                ! ip, ih: the last pressure and enthalpy indices of the tables
    real(R8) :: q, lo, hi

    ip = ubound ( h_tab2D, 1 )
    ih = ubound ( h_tab2D, 2 )
    q  = min ( max ( p, pmin ), pmin + ( ip - 1d-6 ) * deltap )
    i  = int ( ( q - pmin ) / deltap )
    lo = max ( Tmin2(i), Tmin2(i+1) )
    hi = min ( Tmin2(i) + ( ih - 1d-6 ) * deltaT(i), Tmin2(i+1) + ( ih - 1d-6 ) * deltaT(i+1) )
    ok = lo <= hi
    inside = ok .and. q == p .and. T >= lo .and. T <= hi
    p = q
    T = min ( max ( T, lo ), hi )

  end subroutine Table_Bounds


  !> The lowest and the highest temperature of the fluid tables over all their pressures, and the range of their
  !> pressures: for the check of the values of the boundary file, which give a temperature or a pressure alone
  subroutine Table_Temperatures ( Tlo, Thi, plo, phi )
    use FLINT_Lib_Thermodynamic, only: pmin, deltap, Tmin2, deltaT, h_tab2D
    implicit none
    real(R8), intent(out) :: Tlo, Thi, plo, phi
    ! Local
    integer :: ip, ih

    ip  = ubound ( h_tab2D, 1 )
    ih  = ubound ( h_tab2D, 2 )
    Tlo = minval ( Tmin2(0:ip) )
    Thi = maxval ( Tmin2(0:ip) + ih * deltaT(0:ip) )
    plo = pmin
    phi = pmin + ip * deltap

  end subroutine Table_Temperatures

end module ARES_Load_ThermoTransport