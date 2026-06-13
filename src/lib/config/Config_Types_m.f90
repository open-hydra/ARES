module ARES_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Parameters_m

  implicit none
  private

  !! ------------------------------------------------------
  !! Simulation Parameters --------------------------------
  !! ------------------------------------------------------
  type :: simulation_parameters_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    logical   :: newrun           ! Flag di restart
    real(R8)  :: res_threshold    ! Residuo min. per arresto esecuzione
    real(R8)  :: time_threshold   ! Tempo max per arresto esecuzione 
    integer   :: iter_threshold            ! Numero max iterate per arresto esecuzione
    character(len=clen) :: simulation_type ! Type of simulation (euler, laminar, turbulent) 
    ! Useful variables
    integer         :: iter_general     ! Number of iteration - including all MG levels
    integer         :: iter_from_call
    real(R8)        :: time_from_call
    real(R8), allocatable :: residuotot(:)
    integer         :: nthreads         ! Number of threads for simulation
    real(R8)        :: cputime(2)       ! Simulation time duration
    integer         :: TODO             ! Decide solve and/or postprocess
    logical         :: HYDRA_time_accurate = .false.
    logical         :: HYDRA_postprocess   = .false.
    logical         :: HYDRA_MG = .false.
  end type simulation_parameters_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Input-Output -----------------------------------------
  !! ------------------------------------------------------
  type :: io_t
    character(len=llen)  :: warning_message
    character(len=llen)  :: error_message
    character(len=llen)  :: description
    ! USER-DEFINED INPUTS
    integer              :: sol_diter, res_diter   ! iteration interval to save the solution
    real(R8)             :: sol_dtime              ! time interval to save the solution
    logical              :: sol_overwrite      ! switch to overwrite the solution
    character(len=llen)  :: sol_format, ini_format   ! solution format (native,tecplot,vtk) and formatting (ascii,raw)
    character(len=llen)  :: sol_variables, wall_variables   ! Variables to be printed in solution and wall files
    integer              :: shell_diter  ! Shell update
    integer              :: ini_diter    ! input.ini update
    ! Useful variables
    character(len=llen)  :: nameinit  ! Initial file name
    logical              :: write_thermo, write_transport
    character(len=hlen)  :: Ovarnames, ORANSname
    integer              :: Onvar
    real(R8)             :: IOtime
    logical              :: write_wall, write_wall_mechanical, write_wall_thermal, write_wall_mass
    character(len=hlen)  :: Owallnames
    integer              :: Onwall
    integer              :: unitRES             ! Residuals file unit
    character(len=llen)  :: unitRES_format      ! Format for residual_history_file
  end type io_t

  type :: io_probes_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8)            :: dtime
    integer             :: diter
    character(len=clen) :: file
    character(len=hlen) :: varnames
    integer             :: iloc(4)
    real(R8)            :: loc(3)
    ! Useful variables
    ! ...
  end type io_probes_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! NUMERICAL SCHEME !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Time Scheme ------------------------------------------
  !! ------------------------------------------------------
  type :: time_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    character(len=llen) :: solver_type      ! Solver type (explicit/implicit)
    real(R8)            :: cfl              ! Parametro di stabilita' convettiva
    real(R8)            :: vnn              ! Parametro di stabilita' diffusiva
    integer             :: rampa_cfl_iter   ! Numero di iterazioni per rampa di cfl
    logical             :: time_accurate    ! Flag per integrazione time-accurate
    integer             :: n_RK
    character(len=llen) :: integration_variables ! Integration variables (cons/prim)
    ! Useful variables
    ! ...
  end type time_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Implicit residual smoothing --------------------------
  !! ------------------------------------------------------
  type irs_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8)  :: beta
    ! Useful variables
    logical   :: enabled
  end type irs_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Weiss-Smith Preconditioning --------------------------
  !! ------------------------------------------------------
  type prec_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8)  :: Uref, emin
    ! Useful variables
    logical   :: enabled
  end type prec_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Space Discretization ---------------------------------
  !! ------------------------------------------------------
  type :: space_scheme_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    character(len=llen) :: space_reconstruction ! Space reconstruction method
    character(len=llen) :: flux_limiter         ! Flux limiter for space reconstruction
    ! Useful variables
    logical :: SD
  end type space_scheme_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Riemann solver ---------------------------------------
  !! ------------------------------------------------------
  type :: riemann_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    real(R8) :: Minf
    ! Useful variables
    logical  :: SD
  end type riemann_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Shock Detector ---------------------------------
  !! ------------------------------------------------------
  type :: shock_detector_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
    logical :: SD
  end type shock_detector_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! Multigrid --------------------------------------------
  !! ------------------------------------------------------
  type :: multigrid_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    integer              :: MGL                ! Number of multigrid levels
    integer, allocatable :: iter_threshold(:)  ! Number of iterations for each level
    ! Useful variables
    integer, public :: MG_level    ! Current MG-level being solved
    logical, public :: change_MG   ! Logical variable to change MG-level
  end type multigrid_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------

  !! ------------------------------------------------------
  !! BC ---------------------------------------------------
  !! ------------------------------------------------------
  type :: io_bc_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
    logical, allocatable :: viscous_flag(:,:)
    logical, allocatable :: coupling_flag(:,:)
  end type io_bc_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! PHYSICS !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Turbulence closure -----------------------------------
  !! ------------------------------------------------------
  type :: rans_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! USER-DEFINED INPUTS
    real(R8) :: Sc     ! Schmidt laminare
    real(R8) :: Sct    ! Schmidt turbolento
    real(R8) :: Prt    ! Prandtl turbolento
    character(len=llen) :: rans_name
    logical :: SAcomp, SpalartShur, SAR, rough
    logical :: QCR2000, blowing_corr, k_energy_coupling
    logical :: Prt_correction   ! Turbulent Prandtl correction for wall roughness
    ! Useful variables
    logical :: RSM, SD
  end type rans_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------
  !! ------------------------------------------------------
  
  !! Rotating-frame -----------------------------------
  !! ------------------------------------------------------
  type :: rot_t
    character(len=llen)   :: warning_message
    character(len=llen)   :: error_message
    character(len=llen)   :: description
    ! USER-DEFINED INPUTS
    logical               :: enabled
    real(R8)              :: omega                   !> Angular speed [rad/s]
    real(R8)              :: axis(3)                 !> Unit rotation axis
    real(R8)              :: origin(3)               !> A point on the axis [m]
    real(R8)              :: Omega_vec(3) = 0.0_R8   !> omega * axis  [rad/s]
    character(len=128)    :: axis_str      !> String read from registry, parsed in setup
    character(len=128)    :: origin_str    !> String read from registry, parsed in setup
    !> Number and string descriptors for stationary wall patches (registry-populated).
    !> Each stationary_face_str(k) = 'block_index face_direction', parsed in setup.
    integer :: n_stationary = 0
    character(len=32), allocatable :: stationary_face_str(:)
    !> Stationary wall patches in the lab frame, stored as (2, N) pairs:
    !>   stationary_face(1, k) = block index
    !>   stationary_face(2, k) = face direction (1-6)
    !> When ROTATING_FRAME is active, patches in this list receive
    !>   w_wall = -(Omega x r)  (casing/hub); all other type-5/6 walls use
    !>   w_wall = 0             (blades, default).
    !> If unallocated, all type-5/6 walls are treated as rotating (default).
    integer, allocatable :: stationary_face(:,:)
  end type rot_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! OTHER TYPES !!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! ------------------------------------------------------
  !! Basic  -----------------------------------------------
  !! ------------------------------------------------------
  type :: basic_t
    character(len=llen) :: warning_message
    character(len=llen) :: error_message
    character(len=llen) :: description
    ! USER-DEFINED INPUTS
    ! ...
    ! Useful variables
    ! ...
  end type basic_t
  !! ------------------------------------------------------
  !! ------------------------------------------------------


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(simulation_parameters_t), public :: obj_sim_param
  type(io_t), public                    :: obj_io
  type(io_probes_t), allocatable, public:: obj_io_probes(:)
  type(io_bc_t), public                 :: obj_io_bc
  type(time_scheme_t), public           :: obj_time_scheme
  type(irs_t), public                   :: obj_irs
  type(space_scheme_t), public          :: obj_space_scheme
  type(riemann_t), public               :: obj_riemann
  type(prec_t),public                   :: obj_prec
  type(multigrid_t), public             :: obj_multigrid
  type(rans_t), public                  :: obj_rans
  type(rot_t),  public                  :: obj_rot
  type(basic_t), public                 :: obj_thermo
  type(basic_t), public                 :: obj_transport
  type(basic_t), public                 :: obj_thermoinversion
  !! ------------------------------------------------------
  !! ------------------------------------------------------

end module ARES_Config_Types_m