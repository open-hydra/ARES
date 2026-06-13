module ARES_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: ARES_postprocess

  integer, private :: id_stampa=0
  character(81), private :: RANS_shell_format  =    "('ARES | Iter =', i9, ' | Global iter =', i9, ' | Pressure residual =', E13.6)"
  character(69), private :: URANS_shell_format =    "('ARES | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"
  character(97), private :: RANS_shell_format_MG  = "('ARES Grid Level', i2, ' | Iter =', i9, ' | Global iter =', i9, ' | Pressure residual =', E13.6)"
  character(84), private :: URANS_shell_format_MG = "('ARES Grid Level', i2, ' | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"

contains

  subroutine ARES_postprocess ( simulation )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m
    use ARES_Global_m
    use ARES_Parameters_m
    use ARES_IO_Wall,        only: IOwall, Write_Wall_Solution
    use ARES_IO_Solution,    only: Write_Solution
    use ARES_Mod_Diagnostic, only: Write_Diagnostic
    use ARES_IO_Probes,      only: Write_Probes_Data
    use ARES_Mod_Multigrid,  only: Prolongation
    use ARES_Read_Ini,       only: Read_Inifile_Runtime
    use ARES_Mod_MPI,        only: mpi_is_root
    use IR_precision
    implicit none
    type(ARES_simulation_type), intent(inout) :: simulation
    ! Local
    character(len=llen) :: solfile, dgsfile, wallfile, mgsol, mgwall
    real(R8) :: sim_time
    integer  :: m, level

    level = obj_multigrid%MG_level

    ! Update IOfield solutiontime
    if ( obj_time_scheme%time_accurate ) then
      simulation%IOfield(1)%solutiontime = simulation%domain(1)%time
      IOwall%solutiontime = simulation%domain(1)%time
    else
      simulation%IOfield(1)%solutiontime = -real(obj_sim_param%iter_general,8)
      IOwall%solutiontime = -real(obj_sim_param%iter_general,8)
    endif

    if (level == 1) then
      solfile  = '/'//trim(ARES_phase_prefix)//'field'
      wallfile = '/'//trim(ARES_phase_prefix)//'wall'
      dgsfile  = '/'//trim(ARES_phase_prefix)//'diagnostic'
    else
      solfile  = '/'//trim(ARES_phase_prefix)//'field-level'//trim(str(.true.,level))
      mgsol    = '/'//trim(ARES_phase_prefix)//'field-prolongated'
      mgwall   = '/'//trim(ARES_phase_prefix)//'wall-prolongated'
    endif

    ! END OF SIMULATION
    if ( obj_sim_param%TODO == 3 ) then

      if (mpi_is_root) then
        ! Residual history
        if ( mod (simulation%domain(1) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(1)%time, obj_sim_param%residuotot
        endif
      end if

      !if (obj_rot%enabled .and. level == 1 .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0d0) &
      !  call RF_Torque_Walls( simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time )

      if (mpi_is_root) then
        ! Shell
        if ( obj_time_scheme%time_accurate ) then
          if (obj_multigrid%MGL > 1 ) then
            write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
          else
            write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
          endif
        else
          if (obj_multigrid%MGL > 1 ) then
            write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          else
            write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
          endif
        endif

        ! Calculate time at end of simulation
        call Cpu_Time ( obj_sim_param%cputime(2) )

        sim_time = ( obj_sim_param%cputime(2) - obj_sim_param%cputime(1) ) / obj_sim_param%nthreads
        write(*,*)
        write(*,*) '  Time of operation was', sim_time/60, 'min'
      end if

      ! Write output solution (collective — all ranks participate in gather)
      call Write_Solution ( simulation%domain(1), simulation%IOfield(1), solfile )
      if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), wallfile )
      if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), dgsfile )
    

    ! INTERMEDIATE SOLUTION EVALUATION
    elseif ( obj_sim_param%TODO == 2 ) then

      if (level == 1) then
        id_stampa = id_stampa + 1
        if (.not.obj_io%sol_overwrite) solfile  = trim(solfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%sol_overwrite) wallfile = trim(wallfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%sol_overwrite) dgsfile  = trim(dgsfile)//trim(str(.true.,id_stampa))
      endif

      ! Residuals history
      if (mpi_is_root) then
        if ( mod(simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
        endif
      end if

      ! Probes
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Rotating frame torque
      !if (level == 1 .and. obj_rot%enabled .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0d0) &
      !  call RF_Torque_Walls(simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Shell (root only)
      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (obj_multigrid%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            endif
          endif
        endif
      end if
      
      ! Solution output (collective — all ranks participate in gather)
      ! Iter-based solution
      if ( mod(simulation%domain(level) % iter, obj_io%sol_diter) == 0d0 ) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
        endif
        if (level == 1) then
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
          if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
        endif
      ! Time-based solution
      elseif ( simulation%domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime ) then
        obj_sim_param%time_from_call = simulation%domain(level) % time
        if (mpi_is_root) write(*,*) ' ... writing time-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
        endif
        if (level == 1) then
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
          if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
        endif
      endif

      if (obj_multigrid%change_MG) then
        if (mpi_is_root) write(*,*) ' ... writing solution of grid-level ', level
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_multigrid%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
        endif
      endif

      if (obj_multigrid%change_MG) then
        obj_multigrid%MG_level = obj_multigrid%MG_level - 1
        obj_sim_param%iter_from_call = 0
        obj_sim_param%time_from_call = 0.0
      endif


    ! AUXILIARY OPERATIONS AT EACH ITERATION
    elseif ( obj_sim_param%TODO == 1 ) then

      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
        end if
      end if

      ! Probes
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time )

      ! Rotating frame torque
      !if (level == 1 .and. obj_rot%enabled .and. mod(simulation%domain(1)%iter, obj_io%res_diter) == 0d0) &
      !  call RF_Torque_Walls(simulation%domain(1), obj_sim_param%iter_general, simulation%domain(1)%time)

      ! Shell
      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%shell_diter) == 0d0 ) then
          if (obj_time_scheme%time_accurate) then
            if (obj_multigrid%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write (*,URANS_shell_format) obj_sim_param%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (obj_multigrid%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            else
              write (*,RANS_shell_format) obj_sim_param%iter_from_call, obj_sim_param%iter_general, obj_sim_param%residuotot(1)
            endif
          endif
        endif
      end if

    endif

    ! Update input data
    if ( mod (simulation%domain(level) % iter, obj_io%ini_diter) == 0d0 ) then
      call Read_Inifile_Runtime()
    end if
    

    ! PRINT FOR HYDRA-COUPLING PROBLEMS
    if ( obj_sim_param%HYDRA_postprocess ) then
      if (mpi_is_root) write(*,*) ' ... writing HYDRA-coupling solution'
      call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
      if (obj_multigrid%MGL > 1 .and. level > 1) then
        if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
        do m = obj_multigrid%MG_level, 2, -1
          call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
        enddo
        call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
      endif
      if (obj_io % write_wall .and. level == 1) call Write_Wall_Solution ( simulation%domain(1), wallfile )
    endif

  end subroutine ARES_postprocess

end module ARES_Wrap_Postprocess