! FLINT wrapper to load thermodynamic data
module ARES_Load_ThermoTransport
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Load_ThermoTransport

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

end module ARES_Load_ThermoTransport