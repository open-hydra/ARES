module ARES_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: ARES_postprocess

  logical, private :: irs_runtime_warned = .false.
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
    use ARES_IO_Probes,      only: Write_Probes_Data, Close_Probes
    use ARES_Mod_Multigrid,  only: Prolongation
    use ARES_Read_Ini,       only: Read_Inifile_Runtime
    use ARES_Mod_MPI,        only: mpi_is_root
    use ARES_Mod_Timers,     only: timer_summary
    use ARES_Lib_Newstate,   only: Report_Negative_Turbulence
    use IR_precision
    implicit none
    type(ARES_simulation_type), intent(inout) :: simulation
    ! Local
    character(len=llen) :: solfile, dgsfile, wallfile, mgsol, mgwall
    character(len=llen) :: numsol, numwall, numdgs
    real(R8) :: sim_time
    integer  :: m, level
    logical  :: iter_based

    level = obj_multigrid%MG_level

    ! Update IOfield and IOwall solutiontime, on every grid level: a coarse level writes its own solutions
    ! (field-level<L>, wall-level<L>) with its IOfield and IOwall
    do m = 1, obj_multigrid%MGL
      if ( obj_time_scheme%time_accurate ) then
        simulation%IOfield(m)%solutiontime = simulation%domain(m)%time
        IOwall(m)%solutiontime = simulation%domain(m)%time
      else
        simulation%IOfield(m)%solutiontime = -real(obj_sim_param%iter_general,8)
        IOwall(m)%solutiontime = -real(obj_sim_param%iter_general,8)
      endif
    enddo

    if (level == 1) then
      solfile  = '/'//trim(ARES_phase_prefix)//'field'
      wallfile = '/'//trim(ARES_phase_prefix)//'wall'
      dgsfile  = '/'//trim(ARES_phase_prefix)//'diagnostic'
    else
      solfile  = '/'//trim(ARES_phase_prefix)//'field-level'//trim(str(.true.,level))
      wallfile = '/'//trim(ARES_phase_prefix)//'wall-level'//trim(str(.true.,level))
      dgsfile  = '/'//trim(ARES_phase_prefix)//'diagnostic-level'//trim(str(.true.,level))
      mgsol    = '/'//trim(ARES_phase_prefix)//'field-prolongated'
      mgwall   = '/'//trim(ARES_phase_prefix)//'wall-prolongated'
    endif

    ! END OF SIMULATION
    if ( obj_sim_param%TODO == 3 ) then

      if (mpi_is_root) then
        ! Residual history
        if ( mod (simulation%domain(1) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(1)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
        endif
      end if

      ! Probes: the rows due at the last iteration, as at any other
      if (level == 1) call Write_Probes_Data( obj_sim_param%iter_general, simulation%domain(1)%time )

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

      ! The solution due at the last iteration (sol-diter, else sol-dtime), as at any other output iteration (TODO == 2):
      ! it takes the next output number and, when due by time only, sets the time of the last output, which the final
      ! solution below then holds for a restart. With sol-overwrite = true the final solution is that output, written
      ! once. Collective: every rank evaluates the same condition. A run ends on the finest level only (Explicit_Step
      ! sets TODO = 3 at level 1), whose domain and file names the block takes: the guard keeps it there
      iter_based = ( mod (simulation%domain(1) % iter, obj_io%sol_diter) == 0 )
      if ( level == 1 .and. ( iter_based .or. &
           simulation%domain(1) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime ) ) then
        obj_io%output_number = obj_io%output_number + 1
        if ( .not. iter_based ) obj_sim_param%time_from_call = simulation%domain(1) % time
        if ( .not. obj_io%sol_overwrite ) then
          if (mpi_is_root) then
            if ( iter_based ) then
              write(*,*) ' ... writing iter-based solution'
            else
              write(*,*) ' ... writing time-based solution'
            endif
          endif
          numsol  = trim(solfile)//trim(str(.true.,obj_io%output_number))
          numwall = trim(wallfile)//trim(str(.true.,obj_io%output_number))
          numdgs  = trim(dgsfile)//trim(str(.true.,obj_io%output_number))
          call Write_Solution ( simulation%domain(1), simulation%IOfield(1), numsol )
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), numwall )
          if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), numdgs )
        endif
      endif

      ! Timing summary (collective — every rank contributes its accumulators)
      call timer_summary()

      ! Negative turbulence variables after a stage, over the run (collective)
      call Report_Negative_Turbulence()

      ! Write output solution (collective — all ranks participate in gather)
      call Write_Solution ( simulation%domain(1), simulation%IOfield(1), solfile )
      if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), wallfile )
      if (.not. obj_time_scheme%time_accurate) call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), dgsfile )

      ! The history files hold their last rows: close them before the ranks finalize MPI (a unit left
      ! open is closed by the Fortran run-time library only when the process ends, after MPI_FINALIZE)
      call Close_Probes()
      if (mpi_is_root) close(obj_io%unitRES)
    

    ! INTERMEDIATE SOLUTION EVALUATION
    elseif ( obj_sim_param%TODO == 2 ) then

      if (level == 1) then
        ! The output number is stored with every solution: a restart goes on with it
        obj_io%output_number = obj_io%output_number + 1
        if (.not.obj_io%sol_overwrite) solfile  = trim(solfile)//trim(str(.true.,obj_io%output_number))
        if (.not.obj_io%sol_overwrite) wallfile = trim(wallfile)//trim(str(.true.,obj_io%output_number))
        if (.not.obj_io%sol_overwrite) dgsfile  = trim(dgsfile)//trim(str(.true.,obj_io%output_number))
      endif

      ! Residuals history
      if (mpi_is_root) then
        if ( mod(simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
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
      
      ! Solution output (collective — all ranks participate in gather). Every field write of a level also
      ! writes its wall quantities: wall on the finest level, wall-level<n> on a coarse one
      ! Iter-based solution
      if ( mod(simulation%domain(level) % iter, obj_io%sol_diter) == 0d0 ) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
        if (obj_multigrid%MGL > 1 .and. level > 1 .and. prolongated_due()) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution( simulation%domain(1), mgwall )
          obj_multigrid%prolongated_last = simulation%domain(level) % iter
        endif
        if (level == 1 .and. .not. obj_time_scheme%time_accurate) &
          call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
      ! Time-based solution
      elseif ( simulation%domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime ) then
        obj_sim_param%time_from_call = simulation%domain(level) % time
        if (mpi_is_root) write(*,*) ' ... writing time-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
        if (obj_multigrid%MGL > 1 .and. level > 1 .and. prolongated_due()) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = obj_multigrid%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
          if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
          obj_multigrid%prolongated_last = simulation%domain(level) % iter
        endif
        if (level == 1 .and. .not. obj_time_scheme%time_accurate) &
          call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile )
      endif

      if (obj_multigrid%change_MG) then
        if (mpi_is_root) write(*,*) ' ... writing solution of grid-level ', level
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
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
        obj_multigrid%prolongated_last = 0
        obj_sim_param%iter_from_call = 0
        obj_sim_param%time_from_call = 0.0
      endif


    ! AUXILIARY OPERATIONS AT EACH ITERATION
    elseif ( obj_sim_param%TODO == 1 ) then

      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (obj_io%unitRES,trim(obj_io%unitRES_format)) obj_sim_param%iter_general, simulation%domain(level)%time, obj_sim_param%residuotot
          flush(obj_io%unitRES)
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
      call Read_Inifile_Runtime( simulation%domain(level)%iter )
      ! level<n>-iter (multigrid) and iter-threshold (one level) take effect at runtime: the stop
      ! test (Mod_Explicit) uses each level's itermax, set from them at setup.
      if ( obj_multigrid%MGL > 1 ) then
        do m = 1, obj_multigrid%MGL
          simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
        enddo
      else
        simulation%domain(1)%itermax = obj_sim_param%iter_threshold
      endif
      ! IRS can only be switched on at setup, where its arrays are allocated: irs = true, or irs-beta > 0, added to
      ! a run started without it is ignored (irs = false and irs-beta do act on a run started with IRS).
      if ( ( obj_irs%enabled .or. obj_irs%beta > 0d0 ) .and. .not. obj_irs%at_setup ) then
        if ( mpi_is_root .and. .not. irs_runtime_warned ) write(*,'(A)') &
          ' [WARNING] IRS cannot be switched on at runtime (start with irs-beta > 0): irs and irs-beta ignored'
        irs_runtime_warned = .true.
        obj_irs%enabled = .false.
      endif
    end if
    

    ! PRINT FOR HYDRA-COUPLING PROBLEMS
    if ( obj_sim_param%HYDRA_postprocess ) then
      if (mpi_is_root) write(*,*) ' ... writing HYDRA-coupling solution'
      call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
      if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(level), wallfile )
      if (obj_multigrid%MGL > 1 .and. level > 1) then
        if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
        do m = obj_multigrid%MG_level, 2, -1
          call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
        enddo
        call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgsol )
        if (obj_io % write_wall) call Write_Wall_Solution ( simulation%domain(1), mgwall )
      endif
    endif

  contains

    ! A coarse level writes its prolongated solution (field-prolongated, wall-prolongated) with every field
    ! write when prolongated-diter is 0, else only once prolongated-diter of its iterations have passed since
    ! the last one. The end of a level (change_MG) always writes it.
    logical function prolongated_due ()
      prolongated_due = obj_multigrid%prolongated_diter == 0 .or. &
        simulation%domain(level) % iter - obj_multigrid%prolongated_last >= obj_multigrid%prolongated_diter
    end function prolongated_due

  end subroutine ARES_postprocess

end module ARES_Wrap_Postprocess