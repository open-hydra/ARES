module ARES_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ir_precision
  
  implicit none
  private
  public :: Assign_Setup

contains

  subroutine Assign_Setup()
    use ARES_Config_Types_m
    use ARES_Global_m,            only: model, Uref, emin
    use ARES_IO_Solution,         only: Setup_Input_Solution
    use ARES_Mod_Space,           only: Setup_Space_Scheme
    use ARES_Mod_Riemann,         only: Assign_Riemann_Solver
    use ARES_Mod_Newstate,        only: Assign_Integration_Variables
    use ARES_Mod_RANS,            only: Setup_RANS_Model
    !use ARES_Lib_RotatingFrame,   only: Setup_RotatingFrame
    use ARES_Lib_Preconditioning
    implicit none

    ! Setting simulation type
    if (obj_sim_param%simulation_type=='euler') model = 0
    if (obj_sim_param%simulation_type=='laminar') model = 1
    if (obj_sim_param%simulation_type=='turbulent') model = 2

    ! Auto-enable preconditioning when integration variables are preconditioned
    if ( trim(obj_time_scheme%integration_variables) == 'prec' ) then
      obj_prec%enabled = .true.
    else
      obj_prec%enabled = .false.
    end if

    ! Copy preconditioning parameters to global module
    Uref = obj_prec%Uref
    emin = obj_prec%emin

    ! Validate preconditioning parameters: must be < 0 (auto) or > 0 (explicit value)
    if ( obj_prec%enabled ) then
      if ( Uref == 0.0d0 ) then
        obj_prec%error_message = '[ERROR] preconditioning-Uref must be < 0 (auto=sound speed) or > 0 (explicit value). Zero is not allowed.'
      end if
      if ( emin == 0.0d0 ) then
        obj_prec%error_message = '[ERROR] preconditioning-eps-min must be < 0 (auto=default) or > 0 (explicit value). Zero is not allowed.'
      end if
    end if

    ! Setting input solution
    call Setup_Input_Solution()

    ! Space
    call Setup_Space_Scheme()
    call Assign_Riemann_Solver()

    ! Validate Riemann solver for preconditioned integration
    if ( obj_prec%enabled ) then
      if ( trim(obj_riemann%description) /= 'Preconditioned HLLC' ) then
        obj_riemann%error_message = '[ERROR] Preconditioning requires Riemann solver "HLLC Prec". ' // &
          'Current selection: "' // trim(obj_riemann%description) // '" is not compatible.'
      end if
    else
      if ( trim(obj_riemann%description) == 'Preconditioned HLLC' ) then
        obj_riemann%error_message = '[ERROR] Riemann solver "HLLC Prec" requires preconditioned integration variables (integration-variables = prec).'
      end if
    end if

    ! Time
    if (obj_time_scheme%solver_type /= 'euler') then
      read(obj_time_scheme%solver_type(3:3), *) obj_time_scheme%n_rk
    else
      obj_time_scheme%n_rk = 1
    end if
    if (obj_irs%beta>0d0) obj_irs%enabled = .true.
    if (obj_irs%enabled .and. obj_irs%beta<=0d0) then
      obj_irs%warning_message = '[WARNING] IRS enabled but irs-beta is not set (or <= 0). IRS will have no effect.'
      obj_irs%enabled = .false.
    end if
    call Assign_Integration_Variables()

    ! Assign RANS model
    if (model==1) obj_rans%rans_name = 'laminar'
    call Setup_RANS_Model()

    ! Assign Rotating frame
    !call Setup_RotatingFrame()

    ! Set preconditioning description
    if ( obj_prec%enabled ) then
      obj_prec%description = 'Weiss-Smith preconditioning enabled'
    else
      obj_prec%description = 'none'
    end if

    !! Descriptions, warnings and errors

    ! Simulation type
    if (obj_sim_param%simulation_type == 'euler') then
      obj_sim_param%description = 'Euler'
    else if (obj_sim_param%simulation_type == 'laminar') then
      obj_sim_param%description = 'Navier-Stokes'
    else if (obj_sim_param%simulation_type == 'turbulent') then
      obj_sim_param%description = 'Navier-Stokes'
    end if
    ! Time scheme
    if (obj_time_scheme%solver_type == 'euler') then
      obj_time_scheme%description = 'Explicit Euler'
    else if (obj_time_scheme%solver_type == 'RK2') then
      obj_time_scheme%description = 'Second-order Runge-Kutta'
    else if (obj_time_scheme%solver_type == 'RK3') then
      obj_time_scheme%description = 'Third-order Runge-Kutta'
    end if
    if (obj_time_scheme%time_accurate) then
      obj_time_scheme%description = trim(obj_time_scheme%description)//' with time-accurate switch enabled'
    end if
    if (obj_irs%enabled) then
      obj_irs%description = 'Beta set to '//trim(str(.true.,real(obj_irs%beta)))
    end if
    ! Space scheme
    ! ... written in Mod_Space ...
    ! Transport
    if (model>0 .and. obj_transport%description=='Unavailable') &
    write(*,'(A)') '[ERROR] Transport properties are unavailable for the selected phase: cannot run Navier-Stokes simulation'

  end subroutine Assign_Setup

end module ARES_Assign_Setup