module ARES_Read_Sim_Param
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: Register_Sim_Param


contains


  subroutine Register_Sim_Param()
    use ARES_Input_Registry
    use ARES_Config_Types_m
    use ARES_Global_m
    implicit none
    ! Local
    character(len=:), allocatable :: section

    obj_sim_param%warning_message = 'none'
    obj_sim_param%error_message   = 'none'
    obj_sim_param%description     = 'none'

    section = trim(codename)//'-Parameters'

    call reg%add(section,'newrun',obj_sim_param%newrun,'true','Start a new simulation (false = restart)','true , false',.false. )
    call reg%add(section,'res-threshold',obj_sim_param%res_threshold,'1e-10','Residual convergence threshold','> 0',.false.)
    call reg%add( section,'time-threshold',obj_sim_param%time_threshold,'1e30','Maximum simulation time','> 0',.false.)
    call reg%add(section,'iter-threshold',obj_sim_param%iter_threshold,'1000000000','Maximum number of iterations','> 0',.false.)
    call reg%add(section,'simulation-type',obj_sim_param%simulation_type,'euler','Type of simulation (euler, laminar, turbulent)','euler , laminar , turbulent',.true.)

  end subroutine Register_Sim_Param


end module ARES_Read_Sim_Param