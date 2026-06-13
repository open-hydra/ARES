module ARES_Read_Numerics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Config_Types_m
  use ARES_Input_Registry
  
  implicit none
  private
  public :: Register_Numerics

contains

  subroutine Register_Numerics (nmgl)
    use ARES_Parameters_m
    use IR_precision
    implicit none
    ! Local
    integer, intent(in) :: nmgl
    character(len=llen) :: section

    section = trim(codename)//'-Numerics'

    !! ------------------------------------------------------
    !! Time Scheme ------------------------------------------
    !! ------------------------------------------------------
    obj_time_scheme%warning_message = 'none'
    obj_time_scheme%error_message   = 'none'
    obj_time_scheme%description     = 'none'

    ! Solver-type
    call reg%add( trim(section), 'time-scheme', obj_time_scheme%solver_type, 'euler', 'Time integration solver', 'euler, RK2, RK3', .false. )

    ! Stability coefficients and related options
    call reg%add( trim(section), 'cfl', obj_time_scheme%cfl, '0.5', 'CFL number', '> 0', .true. )
    call reg%add( trim(section), 'vnn', obj_time_scheme%vnn, '0.3', 'VNN parameter', '> 0', .false. )
    call reg%add( trim(section), 'cfl-rise-threshold', obj_time_scheme%rampa_cfl_iter, '0', 'CFL rise threshold', '>= 0', .false. )

    ! Time-accurate switch
    call reg%add( trim(section), 'time-accurate', obj_time_scheme%time_accurate, '.false.', 'Time accurate switch', 'logical', .true. )

    ! New state
    call reg%add( trim(section), 'integration-variables', obj_time_scheme%integration_variables, 'prim', 'Integration variables (prim/prec)', 'prim, prec', .false. )

    ! Implicit residual smoothing --------------------------
    obj_irs%description     = 'none'
    obj_irs%warning_message = 'none'
    obj_irs%error_message   = 'none'
    call reg%add( trim(section), 'irs', obj_irs%enabled, '.false.', 'Implicit Residual Smoothing', 'logical', .false. )
    call reg%add( trim(section), 'irs-beta', obj_irs%beta, '0.0', 'IRS beta parameter', '>= 0', .false. )


    ! Preconditioning (enabled automatically when integration-variables = prec)
    obj_prec%description     = 'none'
    obj_prec%warning_message = 'none'
    obj_prec%error_message   = 'none'
    obj_prec%enabled         = .false.
    call reg%add( trim(section), 'preconditioning-Uref', obj_prec%Uref, '-1.0', 'Reference velocity for preconditioning (< 0 = sound speed)', '', .false. )
    call reg%add( trim(section), 'preconditioning-eps-min', obj_prec%emin, '-1.0', 'Minimum epsilon for preconditioning (< 0 = default 0.10)', '', .false. )

    !! ------------------------------------------------------
    !! Space Scheme ------------------------------------------
    !! ------------------------------------------------------

    ! Space Discretization ---------------------------------
    obj_space_scheme%description     = 'none'
    obj_space_scheme%warning_message = 'none'
    obj_space_scheme%error_message   = 'none'
    call reg%add( trim(section), 'space-reconstruction', obj_space_scheme%space_reconstruction, 'none', 'Space reconstruction method', 'MUSCL-SD, MUSCL, first-order, none', .true. )
    call reg%add( trim(section), 'flux-limiter', obj_space_scheme%flux_limiter, 'none', 'Flux limiter for space reconstruction', 'vanalbada, minmod, superbee, vanleer, MC, LIMO3, none', .false. )


    ! Riemann solver ---------------------------------------
    obj_riemann%description     = 'none'
    obj_riemann%warning_message = 'none'
    obj_riemann%error_message   = 'none'
    call reg%add( trim(section), 'riemann-solver', obj_riemann%description, 'HLLC', 'Riemann solver', 'Rusanov, PLLF, HLLE, HLLC, HLLC Prec, HLLC Rotated', .false. )
    call reg%add( trim(section), 'riemann-options-Minf', obj_riemann%Minf, '0.0', 'Mach infinity for AUSM+-up', '>= 0', .false. )


    ! Multigrid levels --------------------------------------
    call Register_Multigrid_Levels(nmgl)

    !! ------------------------------------------------------
    !! ------------------------------------------------------

  end subroutine Register_Numerics


  subroutine Register_Multigrid_Levels(nmgl)
    use ARES_Config_Types_m
    use ARES_Parameters_m
    use IR_precision
    implicit none
    integer, intent(in) :: nmgl
    integer :: m
    character(len=llen) :: option

    ! Allocate iteration threshold array for each level
    obj_multigrid%MGL = max(nmgl, 1)
    allocate( obj_multigrid%iter_threshold(obj_multigrid%MGL) )
    obj_multigrid%iter_threshold = 0

    ! Register per-level iteration parameters (level-2-iter, level-3-iter, ...)
    ! HYDRA coupling registers them under a dedicated section.
    if (.not. obj_sim_param%HYDRA_MG) then
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( trim(codename)//'-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false. )
      enddo
    else
      do m = 1, obj_multigrid%MGL
        write(option,'(A5,I0,A5)') 'level', m, '-iter'
        call reg%add( 'HYDRA-Multigrid', trim(option), obj_multigrid%iter_threshold(m), '0', 'Iterations for multigrid level '//trim(str(.true.,m)), '>= 0', .false. )
      enddo
    endif

  end subroutine Register_Multigrid_Levels  

end module ARES_Read_Numerics