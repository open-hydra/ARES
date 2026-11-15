module ARES_Read_IO
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Input_Registry
  use ARES_Config_Types_m
  
  implicit none
  private
  public :: Register_IO_Fields, Register_Probes

contains

  subroutine Register_IO_Fields()
    use ARES_Global_m
    implicit none
    character(len=:), allocatable :: section

    section = trim(codename)//'-IO'

    obj_io%warning_message = 'none'
    obj_io%error_message   = 'none'
    obj_io%description     = 'none'
    
    call reg%add( section, 'ic-format', obj_io%ini_format, 'tecplot ascii', 'Initial condition (INPUT/ic.*) format', 'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )
    call reg%add( section, 'sol-format', obj_io%sol_format, 'tecplot ascii', 'Solution (OUTPUT/field.*) format', 'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )
    call reg%add( section, 'decomposition-map', obj_io%decomposition_map, 'none', &
                  'Block-to-rank map in INPUT/ (none: ARES assigns the blocks)', '', .false. )

    call reg%add( section, 'sol-diter', obj_io%sol_diter, '1000000000', 'Solution output frequency (iterations)', '> 0', .false. )
    call reg%add( section, 'sol-dtime', obj_io%sol_dtime, '1e30', 'Solution output frequency (physical time)', '> 0', .false. )
    call reg%add( section, 'sol-overwrite', obj_io%sol_overwrite, 'true', 'Overwrite the solution file (false = numbered snapshots)', 'true, false', .false. )
    
    ! Define variables to write
    call reg%add( section, 'sol-variables', obj_io%sol_variables, 'thermo', 'Extra solution variable groups to write (thermo, transport)', '', .false. )
    call reg%add( section, 'wall-variables', obj_io%wall_variables, 'mech thermal', 'Wall variable groups to write (mech, thermal)', '', .false. )
    
    ! Residual history file
    call reg%add( section, 'res-diter', obj_io%res_diter, '1', 'Residual-history write frequency (iterations)', '> 0', .false. )

    ! Shell options
    call reg%add( section, 'shell-diter', obj_io%shell_diter, '1', 'Console (shell) update frequency (iterations)', '> 0', .false. )
    call reg%add( section, 'ini-diter', obj_io%ini_diter, '10000', 'input.ini re-read frequency (iterations)', '> 0', .false. )
    call reg%add( section, 'ini-error', obj_io%ini_error, 'warn', &
                  'input.ini read again with an error: warn and keep the values held before the reading, or stop '// &
                  '(at set-up every error stops)', 'warn, stop', .false. )

    ! Wall-clock timing report
    call reg%add( section, 'timer-diter', obj_io%timer_diter, '0', 'Wall-clock timing report frequency (iterations; 0 = off)', '>= 0', .false. )

  end subroutine Register_IO_Fields


  subroutine Register_Probes(n, probes_name)
    use ARES_Global_m
    use ir_precision, only: str
    implicit none
    integer, intent(in) :: n
    character(*), intent(in) :: probes_name(n)
    integer :: p

    allocate(obj_io_probes(n))
    do p = 1, n
      call reg%add( trim(codename)//'-Probes', 'probe'//trim(str(.true.,p)), obj_io_probes(p)%file, 'probe-placeholder', 'Name of the section of the probe, and of its file OUTPUT/<name>.txt (16 characters at most; probe1, probe2, ... one per probe)', '', .false. )
      call Register_One_Probe(p, probes_name(p))
    end do

  end subroutine Register_Probes


  subroutine Register_One_Probe(p, probe_name)
    use ARES_Parameters_m
    use IR_precision
    implicit none
    integer, intent(in) :: p
    character(*), intent(in) :: probe_name

    obj_io_probes(p)%warning_message = 'none'
    obj_io_probes(p)%error_message   = 'none'
    obj_io_probes(p)%description     = 'none'

    ! Vars name
    call reg%add( trim(probe_name), 'variables', obj_io_probes(p)%varnames, 'none', 'Probe variables to write', '', .false. )

    ! Frequency
    call reg%add( trim(probe_name), 'dtime', obj_io_probes(p)%dtime, '1e30', 'Probe output frequency (physical time)', '> 0', .false. )
    call reg%add( trim(probe_name), 'diter', obj_io_probes(p)%diter, '1000000000', 'Probe output frequency (iterations)', '> 0', .false. )

    ! Location
    call reg%add( trim(probe_name), 'index-position', obj_io_probes(p)%iloc, '0 0 0 0', 'Probe location by index: block i j k', '>= 0', .false. )
    call reg%add( trim(probe_name), 'position', obj_io_probes(p)%loc, '0.0 0.0 0.0', 'Probe location by coordinates: x y z (the cell with the nearest centre)', '', .false. )

  end subroutine Register_One_Probe

end module ARES_Read_IO