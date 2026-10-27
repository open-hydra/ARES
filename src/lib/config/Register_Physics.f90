module ARES_Read_Physics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Parameters_m
  use ARES_Config_Types_m, only: obj_rans, obj_rot, obj_sim_param
  use ARES_Input_Registry

  implicit none
  private
  public :: Register_Physics

contains

  subroutine Register_Physics()
    use ARES_Global_m
    use IR_Precision, only: str
    implicit none
    integer :: k
    character(len=:), allocatable :: section

    obj_rans%warning_message = 'none'
    obj_rans%error_message   = 'none'

    !! ------------------------------------------------------
    !! Physical models ------------------------
    !! ------------------------------------------------------
    ! Same keys as [MOSE-Physics]: the equations and the turbulence model
    section = trim(codename)//'-Physics'
    call reg%add(section, 'equations', obj_sim_param%simulation_type, '', 'Gas dynamics equations', 'euler , navier-stokes', .true.)
    call reg%add(section, 'turbulence', obj_rans%model, 'none', 'RANS turbulence model', 'SA, SA-R, SA-RC, SA-QCR2000, SA-rough, SA-rough-QCR2000, SAcomp, SST, Wilcox2006, SSGLRR, SSGLRR-SD, SSGLRR-w2019, SSGLRR-w2019-SD, none', .false.)
    ! Laminar (molecular) transport closure
    call reg%add(section, 'schmidt', obj_sim_param%Sc, '0.0', 'Laminar Schmidt number (unused: ARES solves a single species)', '>= 0', .false.)


    !! ------------------------------------------------------
    !! Turbulence -------------------------------------------
    !! ------------------------------------------------------
    ! The keys of [MOSE-Turbulence], then the options of ARES only
    section = trim(codename)//'-Turbulence'
    call reg%add(section, 'Prt', obj_rans%Prt, '0.90', 'Turbulent Prandtl number', '> 0', .false.)
    call reg%add(section, 'Sct', obj_rans%Sct, '0.90', 'Turbulent Schmidt number', '> 0', .false.)
    call reg%add(section, 'point-implicit', obj_rans%point_implicit, '.false.', 'Point-implicit (Patankar) treatment of turbulence destruction source terms', 'logical', .false.)
    call reg%add(section, 'omega-wall-bc', obj_rans%omega_wall_bc, 'practical', 'Omega wall condition for k-omega models and SSGLRR: Menter practical or asymptotic', 'practical, asymptotic', .false.)
    call reg%add(section, 'k-coupling', obj_rans%k_energy_coupling, '.false.', 'Turbulent kinetic energy coupling', 'logical', .false.)
    call reg%add(section, 'Prt-correction', obj_rans%Prt_correction, '.false.', 'Turbulent Prandtl correction for wall roughness', 'logical', .false.)
    call reg%add(section, 'sst-asymptotic', obj_rans%sst_asymptotic, '.false.', 'Impose the asymptotic solution omega = 6*nu/(beta_1*d^2) on the near-wall cells (SST only)', 'logical', .false.)
    call reg%add(section, 'sst-asymptotic-cells', obj_rans%sst_asymptotic_cells, '1', 'Number of cells off viscous walls where the asymptotic omega is imposed', '>= 1', .false.)

    
    !! ------------------------------------------------------
    !! Rotating Frame ---------------------------------------
    !! ------------------------------------------------------
    !section = trim(codename)//'-Rotating-Frame'
    !call reg%add(section, 'omega',  obj_rot%omega,      '0.0',    'Angular speed [rad/s]', '>= 0', .false.)
    !call reg%add(section, 'axis',   obj_rot%axis_str,   '0.0 0.0 1.0', 'Rotation axis direction (3 components)', '', .false.)
    !call reg%add(section, 'origin', obj_rot%origin_str, '0.0 0.0 0.0', 'Point on the rotation axis [m]', '', .false.)
    ! Stationary face entries (count scanned in Scan_Ini, strings allocated there)
    !do k = 1, obj_rot%n_stationary
    !  call reg%add(section, 'stationary-face-'//trim(str(.true.,k)), obj_rot%stationary_face_str(k), &
    !               '', 'Stationary face: block_index face_direction', '', .false.)
    !end do

    !! ------------------------------------------------------
    !! GSI --------------------------------------------------
    !! ------------------------------------------------------
    !! TODO: Refactor GSI material properties to registry
    !! For now: HDPE, PP, Paraffin, HTPB properties
    !! These require dedicated sections in config types
    !! ======================================================

  end subroutine Register_Physics

end module ARES_Read_Physics