module ARES_Read_Physics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Parameters_m
  use ARES_Config_Types_m, only: obj_rans, obj_rot
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


    !! ------------------------------------------------------
    !! Turbulence -------------------------------------------
    !! ------------------------------------------------------
    section = trim(codename)//'-RANS'
    call reg%add(section, 'turbulence-model', obj_rans%rans_name, '', 'RANS turbulence model', 'SA, SA-R, SA-RC, SA-QCR2000, SA-rough, SA-rough-QCR2000, SAcomp, SST, Wilcox2006, SSGLRR, none', .false.)
    call reg%add(section, 'Prt', obj_rans%Prt, '0.90', 'Turbulent Prandtl number', '> 0', .false.)
    call reg%add(section, 'Sct', obj_rans%Sct, '0.90', 'Turbulent Schmidt number', '> 0', .false.)
    call reg%add(section, 'Sc', obj_rans%Sc, '0.7', 'Schmidt number', '> 0', .false.)
    call reg%add(section, 'k-coupling', obj_rans%k_energy_coupling, '.false.', 'Turbulent kinetic energy coupling', 'logical', .false.)
    call reg%add(section, 'Prt-correction', obj_rans%Prt_correction, '.false.', 'Turbulent Prandtl correction for wall roughness', 'logical', .false.)

    
    !! ------------------------------------------------------
    !! Rotating Frame ---------------------------------------
    !! ------------------------------------------------------
    !section = trim(codename)//'-rotating-frame'
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