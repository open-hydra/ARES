module ARES_Wrap_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: ARES_setup

contains

  subroutine ARES_setup ( simulation )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m
    use ARES_Global_m
    use ARES_Read_Ini,             only: Read_Inifile
    use ARES_Load_ThermoTransport, only: Load_ThermoTransport
    use ARES_Assign_Setup,         only: Assign_Setup
    use ARES_IO_Solution,          only: Read_IC, Setup_Output_Solution
    use ARES_Mod_Allocate_Data,    only: Setup_Data_Structure, deallocate_remote_computation_data
    use ARES_Mod_Multigrid,        only: Setup_Multigrid, Restriction
    use ARES_IO_BC,                only: Setup_BC
    use ARES_Mod_Metrics,          only: Setup_Metrics
    use ARES_IO_Probes,            only: Setup_Probes, Restart_Probes, Drop_Probe_Rows, Drop_Rows_After
    use ARES_IO_Wall,              only: Initialize_Wall_File, IOwall
    use ARES_Lib_Ghost,            only: Fill_Ghost_Cell
    use ARES_Mod_MPI,              only: mpi_is_root, mpi_size_, partition_blocks, mpi_stop_all, mpi_bcast_logical
    use ARES_Mod_Timers,           only: timer_run_begin
    use ARES_Mod_GhostExchange,    only: build_ghost_schedule, build_local_bc_index, &
                                         allocate_exchange_schedules
    implicit none
    type(ARES_simulation_type), intent(inout) :: simulation
    ! Local
    integer :: m, ios
    logical :: resume
    logical :: setup_error   ! an error found by a check of rank 0 (Stop_Together)

    !! ------------------------------------------------------
    !! ------------------------------------------------------
    ! Print Header
    if (mpi_is_root) call Print_Header ()

    call Load_ThermoTransport ()
    setup_error = .false.
    if (mpi_is_root) call Check_Tables ( setup_error )
    call Stop_Together ( setup_error, 'thermodynamic or transport data error (see above)' )

    ! Read input.ini
    call Read_Inifile ()
    setup_error = .false.
    if (mpi_is_root) call Check_Input ( setup_error )
    call Stop_Together ( setup_error, 'input file error (see above)' )

    ! Allocate container objects Domain and IOfield.
    allocate ( simulation%domain  ( obj_multigrid%MGL ) )
    allocate ( simulation%IOfield ( obj_multigrid%MGL ) )
    simulation%domain(1)%itermax = obj_sim_param%iter_threshold

    ! Assign setup
    call Assign_Setup ()

    ! Stop on configuration errors now, before the (possibly long) initial-condition read.
    setup_error = .false.
    if (mpi_is_root) call Stop_Simulation ( setup_error )
    call Stop_Together ( setup_error, 'setup error (see above)' )

    ! Read file for initial solution.
    call Read_IC ( simulation%IOfield(1) )
    ! The ranks decide together on an error of the read (Check_IC)
    call Check_IC ( summary = .true. )
    ! A restart goes on from the time of the solution it reads: the physical time of a time-accurate
    ! solution, minus the global iteration count of a steady one. A solution that holds no time (a vtk
    ! file without TIME, a Tecplot ASCII file without SOLUTIONTIME) cannot tell where to go on: stop
    if ( .not. obj_sim_param%newrun .and. .not. obj_io%restart%time_found ) then
      if (mpi_is_root) then
        write(*,'(A)') ' [ERROR] Restart from '//trim(obj_io%nameinit)//': the solution holds no time '// &
                       '(vtk TIME, Tecplot SOLUTIONTIME), and a restart goes on from the time of its solution.'
        write(*,'(A)') ' [ERROR]   To start from this solution, use it as the initial condition: '// &
                       'INPUT/ic.* with newrun = true.'
      end if
      ! Every rank has read the solution and found no time: the ranks stop together, once rank 0 has written the
      ! message (mpi_stop_all). With MPI_ABORT, another rank could stop rank 0 before it writes.
      call mpi_stop_all ( 'restart from a solution without a time (see above)' )
    end if

    ! Allocation of data structures and copy solution from IOfield.
    call Setup_Data_Structure ( simulation%domain(1), simulation%IOfield(1) )
    ! Every value that the initial condition gives the domain is finite, before the coarse levels are built from it
    call Check_IC_Values ( simulation%domain(1), simulation%IOfield(1) )
    allocate(obj_sim_param%residuotot(nres))

    ! With multigrid, allocate Grid and IOfield 2,...,MGL
    if ( obj_multigrid%MGL > 1 ) then
      call Setup_Multigrid ( simulation )
    end if

    ! Read Boundary Conditions file
    call Setup_BC ( simulation%domain )
    setup_error = .false.
    if (mpi_is_root) call Check_BC ( setup_error )
    call Stop_Together ( setup_error, 'boundary condition error (see above)' )

    ! MPI: partition blocks once on finest-level cell counts (restriction and
    ! prolongation are block-local, so a block must live on the same rank at
    ! every level), then build per-level ghost/chimera communication schedules.
    call partition_blocks(simulation%domain(1)%nb, &
      [(product(simulation%domain(1)%blk(ios)%dim(1:3)), ios=1, simulation%domain(1)%nb)])
    call allocate_exchange_schedules(obj_multigrid%MGL)
    do m = 1, obj_multigrid%MGL
      simulation%domain(m)%mg_level = m
      call build_ghost_schedule(simulation%domain(m))
      call build_local_bc_index(simulation%domain(m))
    end do

    ! Read grid file and setup metrics in each grid level.
    do m = 1, obj_multigrid%MGL
      call Setup_Metrics ( simulation%domain(m) )
    end do

    ! Free heavy arrays on non-local blocks to save memory
    do m = 1, obj_multigrid%MGL
      call deallocate_remote_computation_data(simulation%domain(m))
    end do

    ! Setup probes location, output files, and variables to be printed.
    call Setup_Probes ( simulation%domain(1), obj_sim_param%newrun )
 
    ! Solution setup
    call Setup_Output_Solution ( simulation%IOfield )

    ! Initialize the wall files, one per multigrid level (wall, wall-level<n>)
    allocate ( IOwall(obj_multigrid%MGL) )
    if (mpi_is_root .and. obj_io % write_wall) then
      do m = 1, obj_multigrid%MGL
        call Initialize_Wall_File ( simulation%domain(m), simulation%IOfield(1)%tec%extension)
      enddo
    endif

    ! Print simulation onto the logfile/shell.
    if (mpi_is_root) call Print_Shell_Info ()

    ! Residual history file
    if (mpi_is_root) then
      open(newunit=obj_io%unitRES,file='OUTPUT/'//trim(ARES_phase_prefix)//'residual-history.dat',status='unknown',form='formatted')
      if ( .not. obj_sim_param%newrun ) then
        ios = 0
        do while ( ios == 0 )
          read (obj_io%unitRES, *, iostat=ios)
        enddo
        backspace (obj_io%unitRES)
      endif
    end if

    ! Setting the format for writing residual_history file
    write(obj_io%unitRES_format,'(A11,I0,A7)') '(I8,E20.10,', nres, 'E20.10)'

    ! Initialize run variables depending on scheme
    select case (trim(obj_time_scheme%solver_type))

      case ('euler', 'RK2', 'RK3')

        obj_multigrid%MG_level = obj_multigrid%MGL

        do m = 1, obj_multigrid%MGL
          simulation%domain(m)%iter = 0

          if (obj_sim_param%newrun) then
            if (obj_time_scheme%time_accurate) then
              simulation%domain(m)%time = 0.0d0
            else
              simulation%domain(m)%time = -1.0d0
            endif
            obj_sim_param%iter_general = 0
          else
            ! The restart solution is read on the finest level: its time and iteration count hold for
            ! every level (the coarse IOfields are built from the fine one and hold no time)
            if (obj_time_scheme%time_accurate) then
              if (obj_io%restart%steady) then
                ! A transient started from a steady solution starts at t = 0
                simulation%domain(m)%time = 0.0d0
              else
                simulation%domain(m)%time = simulation%IOfield(1)%solutiontime
              endif
              obj_sim_param%iter_general = 0
            else
              obj_sim_param%iter_general = int(simulation%IOfield(1)%solutiontime)
              simulation%domain(m)%time = -1.d0
            endif
          endif
          obj_sim_param%time_from_call = simulation%domain(m)%time
          obj_sim_param%iter_from_call = 0
        enddo

        ! The counters stored with a solution, in any format: the restart goes on as the uninterrupted run
        ! would (CFL ramp, output cadence and numbering, residual history, multigrid level), if
        ! that run was of the same kind, steady or time-accurate. Otherwise, or with a solution
        ! that holds no counters, the run starts its counters as a new one.
        resume = .false.
        if ( .not. obj_sim_param%newrun ) &
          resume = obj_io%restart%found .and. ( obj_io%restart%steady .neqv. obj_time_scheme%time_accurate )
        if ( resume ) then
          obj_sim_param%iter_general   = obj_io%restart%cycle
          obj_multigrid%MG_level       = max(1, min(obj_multigrid%MGL, obj_io%restart%grid_level))
          ! The solution read is that of the finest level, and so is its iteration count
          simulation%domain(1)%iter    = obj_io%restart%iteration
          obj_sim_param%time_from_call = obj_io%restart%last_output_time
          obj_io%output_number         = obj_io%restart%output_number
          ! The unit of res-threshold, when the solution holds it: the first iteration of the restart keeps it
          obj_sim_param%res_reference  = obj_io%restart%res_reference
          obj_sim_param%res_restored   = obj_io%restart%res_found
          ! A run that went on past the solution read left in residual-history.dat and in the probe files
          ! the rows of the iterations after it, which the restart writes again: drop them
          if (mpi_is_root) then
            close(obj_io%unitRES)
            call Drop_Rows_After ( 'OUTPUT/'//trim(ARES_phase_prefix)//'residual-history.dat', 2 + nres, &
                                   obj_io%restart%cycle, simulation%domain(1)%time )
            open(newunit=obj_io%unitRES,file='OUTPUT/'//trim(ARES_phase_prefix)//'residual-history.dat', &
                 status='old',position='append',form='formatted')
          endif
          call Drop_Probe_Rows ( obj_io%restart%cycle, simulation%domain(1)%time )
        endif
        ! Without it (a solution that holds no counters, or a time-accurate one), a steady restart
        ! takes the pressure residual of its own first iteration, as a new run
        if ( .not. ( obj_sim_param%newrun .or. obj_time_scheme%time_accurate .or. obj_sim_param%res_restored ) ) &
          obj_io%warning_message = '[WARNING] The restart solution holds no residual reference of a steady run: '// &
            'res-threshold is relative to the pressure residual of the first iteration of the restart.'
        ! Time-based probe rows fall on the multiples of their dtime: a restart that goes on in time
        ! writes the next one not reached yet; a new transient starts from the first, as a new run
        if ( .not. obj_sim_param%newrun ) then
          call Restart_Probes ( simulation%domain(1)%time, &
                                resume = .not. ( obj_io%restart%steady .and. obj_time_scheme%time_accurate ) )
        endif
        call Fill_Ghost_Cell ( simulation%domain(1) )

        ! Interpolate initial solution on coarser domains
        if ( obj_multigrid%MGL > 1 ) then
          do m = 2, obj_multigrid%MGL
            call Restriction ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
            call Fill_Ghost_Cell ( simulation%domain(m) )
          end do
        endif

    end select

    ! Print warnings onto the logfile/shell.
    if (mpi_is_root) call Print_Warnings ()

    ! If errors are found, print error messages and stop the simulation.
    setup_error = .false.
    if (mpi_is_root) call Stop_Simulation ( setup_error )
    call Stop_Together ( setup_error, 'setup error (see above)' )

    ! Calculate time at beginning of simulation
    call Cpu_Time ( obj_sim_param%cputime(1) )

    ! Set-up is over: start the solver clock
    call timer_run_begin()

  contains

    subroutine Print_Header()

      write(*,*)
      write(*,'(A82)')'_____/\\\\\\\\\_______/\\\\\\\\\______/\\\\\\\\\\\\\\\_____/\\\\\\\\\\\___        '
      write(*,'(A82)')' ___/\\\\\\\\\\\\\___/\\\///////\\\___\/\\\///////////____/\\\/////////\\\_       '
      write(*,'(A82)')'  __/\\\/////////\\\_\/\\\_____\/\\\___\/\\\______________\//\\\______\///__      '
      write(*,'(A82)')'   _\/\\\_______\/\\\_\/\\\\\\\\\\\/____\/\\\\\\\\\\\_______\////\\\_________     '
      write(*,'(A82)')'    _\/\\\\\\\\\\\\\\\_\/\\\//////\\\____\/\\\///////___________\////\\\______    '
      write(*,'(A82)')'     _\/\\\/////////\\\_\/\\\____\//\\\___\/\\\_____________________\////\\\___   '
      write(*,'(A82)')'      _\/\\\_______\/\\\_\/\\\_____\//\\\__\/\\\______________/\\\______\//\\\__  '
      write(*,'(A82)')'       _\/\\\_______\/\\\_\/\\\______\//\\\_\/\\\\\\\\\\\\\\\_\///\\\\\\\\\\\/___ '
      write(*,'(A82)')'        _\///________\///__\///________\///__\///////////////____\///////////_____'
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,'(A82)')'                                         /\           ((                          '
      write(*,'(A82)')'                                        /  \         )))(                         '
      write(*,'(A82)')'                                       | /\ |       )(((()                        '
      write(*,'(A82)')'                      A                | || |       )))))                         '
      write(*,'(A82)')'                                        \||/         )((                          '
      write(*,'(A82)')'                  Real-fluid             \/          \|/                          '
      write(*,'(A82)')'                                         ||          |=|                          '
      write(*,'(A82)')'                   Equations             ||          | |                          '
      write(*,'(A82)')'                                         ||          | |                          '
      write(*,'(A82)')'                    Solver               ||          | |                          '
      write(*,'(A82)')'                                         ||         /| |\                         '
      write(*,'(A82)')'                                         ||         \| |/                         '
      write(*,'(A82)')'                                         \/          |_|                          '
      write(*,*)
      write(*,*)

                                          
    end subroutine Print_Header


    subroutine Check_Tables ( has_error )
      implicit none
      logical, intent(out) :: has_error

      has_error = .false.

      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Loading'
      write(*,'(A)') ' ========================================================================================='

      ! Thermodynamics
      if (index(obj_thermo%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Thermodynamics', 'FAIL'
        write(*,'(4X,A)') trim(obj_thermo%error_message)
        has_error = .true.
      else
        write(*,'(A,T35,A)') '   Thermodynamics', 'OK'
      endif

      ! Transport
      if (index(obj_transport%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Transport', 'FAIL'
        write(*,'(4X,A)') trim(obj_transport%error_message)
        has_error = .true.
      else
        write(*,'(A,T35,A)') '   Transport', 'OK'
      endif

      ! Thermo inversion table (p,h) -> T
      if (index(obj_thermoinversion%error_message,'ERROR')>0) then
        write(*,'(A,T35,A)') '   Thermo inversion', 'FAIL'
        write(*,'(4X,A)') trim(obj_thermoinversion%error_message)
        has_error = .true.
      else
        write(*,'(A,T35,A)') '   Thermo inversion', 'OK'
      endif

    end subroutine Check_Tables


    !> An error of the read of the initial condition (Read_IC, the reader writes it in obj_io%error_message) stops the run.
    !> Every rank reads the file, and a rank can fail to read what rank 0 reads. So the ranks take the decision together
    !> (the largest of a flag): with summary, rank 0 writes the line of the set-up summary with the verdict of all ranks;
    !> the first rank that found an error writes it (with its rank when the run has more than one), and every rank stops
    !> with mpi_stop_all, exit code 1.
    subroutine Check_IC ( summary )
      use ARES_Mod_MPI, only: mpi_rank_, mpi_allreduce_max_r8, mpi_allreduce_min_r8
      implicit none
      logical, intent(in) :: summary   ! write the line of the set-up summary
      real(R8) :: flag, flag_all, first

      flag = 0.0_R8
      if ( index(obj_io%error_message,'ERROR') > 0 ) flag = 1.0_R8
      call mpi_allreduce_max_r8 ( flag, flag_all )
      if ( summary .and. mpi_is_root ) then
        if ( flag_all > 0.0_R8 ) then
          write(*,'(A,T35,A)') '   Initial conditions', 'FAIL'
        else
          write(*,'(A,T35,A)') '   Initial conditions', 'OK'
        endif
      endif
      if ( flag_all > 0.0_R8 ) then
        call mpi_allreduce_min_r8 ( merge(real(mpi_rank_, R8), huge(1.0_R8), flag > 0.0_R8), first )
        if ( flag > 0.0_R8 .and. nint(first) == mpi_rank_ ) then
          if ( mpi_size_ > 1 ) then
            write(*,'(4X,A,I0,A)') trim(obj_io%error_message)//' (rank ', mpi_rank_, ')'
          else
            write(*,'(4X,A)') trim(obj_io%error_message)
          endif
        endif
        call mpi_stop_all ( 'initial condition error (see above)' )
      endif

    end subroutine Check_IC


    !> The initial condition, or the solution that a restart reads, gives the domain the coordinates of its nodes and the
    !> variables of its cells, and a value there that is not finite (NaN, infinity) stops the run at set-up, with the place
    !> of the value. Not_Finite reads the bits of each value, an operation on integers that raises no floating-point
    !> exception: a build that traps the invalid operations (the DEBUG build; a comparison with a NaN is one) stopped
    !> with SIGFPE at the first such operation on the value, at the end of the set-up, and -ffast-math cannot fold the
    !> test away. Each rank checks the values it holds: the nodes of every block, the cells of the blocks whose state it
    !> holds. The ranks take the decision together (the largest of a flag), the first rank that found a value writes it,
    !> and every rank stops with mpi_stop_all: exit code 1, in every build.
    subroutine Check_IC_Values ( domain, IOfield )
      use Lib_ORION_data,    only: ORION_data
      use ARES_Finite_m,     only: Not_Finite
      use ARES_Mod_MPI,      only: mpi_rank_, mpi_allreduce_max_r8, mpi_allreduce_min_r8
      implicit none
      type(ARES_domain_type), intent(in) :: domain
      type(ORION_data), intent(in)       :: IOfield
      ! Local
      character(len=1), parameter :: axis(3) = ['x', 'y', 'z']
      character(len=hlen) :: what
      character(len=40)   :: name
      integer  :: b, i, j, k, d, off
      real(R8) :: flag, flag_all, first

      what = ''
      blocks: do b = 1, domain%nb
        do k = 0, domain%blk(b)%dim(3)
          do j = 0, domain%blk(b)%dim(2)
            do i = 0, domain%blk(b)%dim(1)
              do d = 1, 3
                if ( Not_Finite(domain%blk(b)%node(i,j,k)%c(d)) ) then
                  write(what,'(3A,I0,A,2(I0,A),I0)') 'coordinate ', axis(d), ' of block ', b, ', node i, j, k = ', &
                    i, ', ', j, ', ', k
                  exit blocks
                end if
              end do
            end do
          end do
        end do
        if ( .not. allocated(domain%blk(b)%P) ) cycle
        do k = 1, domain%blk(b)%dim(3)
          do j = 1, domain%blk(b)%dim(2)
            do i = 1, domain%blk(b)%dim(1)
              do d = 1, nprim
                if ( Not_Finite(domain%blk(b)%P(d,i,j,k)) ) then
                  ! The name of the variable when the file gives the names (Tecplot ASCII), which come after the coordinates
                  name = ''
                  if ( allocated(IOfield%varnames) .and. allocated(IOfield%block(1)%vars) ) then
                    off = size(IOfield%varnames) - size(IOfield%block(1)%vars, 1)
                    if ( off >= 0 .and. off + d <= size(IOfield%varnames) ) &
                      name = ' ('//trim(adjustl(IOfield%varnames(off + d)))//')'
                  end if
                  write(what,'(A,I0,2A,I0,A,2(I0,A),I0)') 'variable ', d, trim(name), ' of block ', b, &
                    ', cell i, j, k = ', i, ', ', j, ', ', k
                  exit blocks
                end if
              end do
            end do
          end do
        end do
      end do blocks

      flag = 0.0_R8
      if ( len_trim(what) > 0 ) flag = 1.0_R8
      call mpi_allreduce_max_r8 ( flag, flag_all )
      if ( flag_all > 0.0_R8 ) then
        ! Written once, by the first rank that found a value (rank 0 when every rank holds it)
        call mpi_allreduce_min_r8 ( merge(real(mpi_rank_, R8), huge(1.0_R8), flag > 0.0_R8), first )
        if ( flag > 0.0_R8 .and. nint(first) == mpi_rank_ ) &
          write(*,'(A)') ' [ERROR] '//trim(obj_io%nameinit)//' holds a value that is not finite (NaN or infinity): '// &
                         trim(what)
        call mpi_stop_all ( 'initial condition error (see above)' )
      end if

    end subroutine Check_IC_Values


    subroutine Check_Input ( has_error )
      use ARES_Input_Registry
      use ARES_Backend_INI, only: Obsolete_Keys
      use Finer,            only: file_ini
      implicit none
      logical, intent(out) :: has_error
      character(len=hlen) :: out
      type(file_ini) :: fini

      ! Keys of the earlier layout first: they also explain a missing required key (equations)
      call fini%load(filename='input.ini')
      out = Obsolete_Keys(fini)
      if (len_trim(out) == 0) out = Validate_Registry()
      has_error = index(out,'ERROR') > 0
      if ( has_error ) then
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(out)
      else
        write(*,'(A,T35,A)') '   Input file', 'OK'
      endif

    end subroutine Check_Input


    subroutine Check_BC ( has_error )
      implicit none
      logical, intent(out) :: has_error

      has_error = index(obj_io_bc%error_message,'ERROR') > 0
      if ( has_error ) then
        write(*,'(A,T35,A)') '   Boundary conditions', 'FAIL'
        write(*,'(4X,A)') trim(obj_io_bc%error_message)
      else
        write(*,'(A,T35,A)') '   Boundary conditions', 'OK'
      endif
      write(*,'(A)') ' ========================================================================================='

    end subroutine Check_BC


    subroutine Print_Warnings
      implicit none

      ! IO warnings
      if (index(obj_io%warning_message,'WARNING')>0) write(*,'(A)') obj_io%warning_message
      ! Physics warnings
      if (index(obj_thermo%warning_message,'WARNING')>0) write(*,'(A)') obj_thermo%warning_message
      if (index(obj_transport%warning_message,'WARNING')>0) write(*,'(A)') obj_transport%warning_message
      if (index(obj_thermoinversion%warning_message,'WARNING')>0) write(*,'(A)') obj_thermoinversion%warning_message
      if (index(obj_rans%warning_message,'WARNING')>0) write(*,'(A)') obj_rans%warning_message
      ! Rotation-curvature correction (-RC, a model the input chooses) in a time-accurate run: a limit of the model
      if (model == 2 .and. ( obj_time_scheme%time_accurate .or. obj_sim_param%HYDRA_time_accurate ) .and. &
          obj_rans%SpalartShur) write(*,'(A)') &
        '[WARNING] The rotation-curvature correction (-RC) omits the time derivative of the strain-rate tensor: '// &
        'it holds for steady flows only.'
      !if (index(obj_rot%warning_message,'WARNING')>0) write(*,'(A)') obj_rot%warning_message
      ! Numerical scheme warnings
      if (index(obj_time_scheme%warning_message,'WARNING')>0) write(*,'(A)') obj_time_scheme%warning_message
      if (index(obj_riemann%warning_message,'WARNING')>0) write(*,'(A)') obj_riemann%warning_message
      if (index(obj_space_scheme%warning_message,'WARNING')>0) write(*,'(A)') obj_space_scheme%warning_message
      if (index(obj_irs%warning_message,'WARNING')>0) write(*,'(A)') obj_irs%warning_message

    end subroutine Print_Warnings


    subroutine Stop_Simulation ( has_error )
      implicit none
      logical, intent(out) :: has_error

      has_error = .false.

      ! Physics errors
      if (index(obj_transport%error_message,'ERROR')>0) then
        write(*,'(A)') obj_transport%error_message;  has_error = .true.
      endif
      if (index(obj_rans%error_message,'ERROR')>0) then
        write(*,'(A)') obj_rans%error_message;       has_error = .true.
      endif
      !if (index(obj_rot%error_message,'ERROR')>0) then
      !  write(*,'(A)') obj_rot%error_message;       has_error = .true.
      !endif
      ! Numerical scheme errors
      if (index(obj_time_scheme%error_message,'ERROR')>0) then
        write(*,'(A)') obj_time_scheme%error_message;  has_error = .true.
      endif
      if (index(obj_riemann%error_message,'ERROR')>0) then
        write(*,'(A)') obj_riemann%error_message;      has_error = .true.
      endif
      if (index(obj_space_scheme%error_message,'ERROR')>0) then
        write(*,'(A)') obj_space_scheme%error_message; has_error = .true.
      endif
      if (index(obj_irs%error_message,'ERROR')>0) then
        write(*,'(A)') obj_irs%error_message;          has_error = .true.
      endif
      if (index(obj_prec%error_message,'ERROR')>0) then
        write(*,'(A)') obj_prec%error_message;         has_error = .true.
      endif

    end subroutine Stop_Simulation


    !> A check of the set-up that runs on rank 0 alone writes what it finds there. Rank 0 broadcasts whether it found an
    !> error, and on an error every rank stops together (mpi_stop_all): rank 0 has written its message before any rank
    !> ends, and the run ends with exit code 1. With MPI_ABORT from rank 0 alone, which makes only a best attempt to stop
    !> the processes (MPI-4.1, section 11.4.2), the launcher can end while rank 0 is writing: with Intel MPI 2021.12,
    !> one run in ten ends so (measured), with exit code 141 (128 + SIGPIPE) and without the line of the abort.
    subroutine Stop_Together ( found, why )
      implicit none
      logical, intent(inout)       :: found
      character(len=*), intent(in) :: why

      call mpi_bcast_logical ( found )
      if ( found ) call mpi_stop_all ( why )

    end subroutine Stop_Together


    subroutine Print_Shell_Info()
      use IR_Precision,   only: str  
      use ARES_IO_BC,     only: Print_BC_Summary
      use ARES_IO_Probes, only: nprobes
      implicit none
      ! Local
      character(llen) :: eosword
      integer :: b, total_cells, k

      ! ----- Domain topology -----
      total_cells = 0
      do b = 1, simulation%domain(1)%nb
        total_cells = total_cells + &
          product(simulation%domain(1)%blk(b)%dim(1:3))
      end do

      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Set-up'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Domain'
      write(*,'(A,T35,I0)') '   Blocks', simulation%domain(1)%nb
      write(*,'(A,T35,I0)') '   Cells', total_cells
      write(*,'(A,T35,I0)') '   Boundary faces', simulation%domain(1)%nbound
      if ( obj_multigrid%MGL > 1 ) &
        write(*,'(A,T35,I0)') '   Multigrid levels', obj_multigrid%MGL

      ! ----- Boundary conditions -----
      call Print_BC_Summary ()

      ! ----- IO -----
      write(*,*)
      write(*,'(A)') ' Input/Output'
      write(*,'(A,T35,A)') '   Initial conditions file', trim(obj_io%nameinit)
      write(*,'(A,T35,A)') '   Solution format', trim(obj_io%sol_format)
      write(*,'(A,T35,A)') '   Probes number', str(.true.,nprobes)
      write(*,'(A)') ' ========================================================================================='

      ! ----- Physical model -----
      eosword = 'Real'

      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Physical model'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Real Fluid model'
      write(*,'(A,T35,A)') '   Equations', trim(obj_sim_param%description)
      if (model==1) then
        write(*,'(A,T35,A)') '   Viscous model', 'Laminar'
      elseif (model==2) then
        write(*,'(A,T35,A)') '   Viscous model', 'Turbulent, '//trim(obj_rans%model)
        if (obj_rans%point_implicit) &
          write(*,'(A,T35,A)') '   Point-implicit', 'ON (Patankar destruction source terms)'
        if (nRANS == 2 .and. trim(obj_rans%omega_wall_bc) == 'asymptotic') &
          write(*,'(A,T35,A)') '   Omega wall condition', 'asymptotic (6*nu/(beta1*y^2))'
        if (obj_rans%Prt_correction) &
          write(*,'(A,T35,A)') '   Prt correction', 'ON (rough-wall turbulent Prandtl)'
        if (obj_rans%sst_asymptotic) &
          write(*,'(A,T35,A,I0,A)') '   Asymptotic wall omega', 'ON, first ', obj_rans%sst_asymptotic_cells, ' cell(s)'
      end if
      write(*,'(A,T35,A)') '   Equation of state', trim(eosword)
      write(*,'(A,T35,A)') '   Thermodynamics', trim(obj_thermo%description)

      !if (obj_rot%enabled) then
      !  write(*,*)
      !  write(*,'(A)')         ' Rotating Frame '
      !  write(*,'(A,F14.4,A)') '   omega  = ', obj_rot%omega,   ' rad/s'
      !  write(*,'(A,3F11.5)')  '   axis   = ', obj_rot%axis
      !  write(*,'(A,3F11.5)')  '   origin = ', obj_rot%origin
      !  if ( allocated(obj_rot%stationary_face) ) then
      !    write(*,'(A,I0,A)') '   stationary patches = ', size(obj_rot%stationary_face,2), ':'
      !    do k = 1, size(obj_rot%stationary_face, 2)
      !      write(*,'(A,I0,A,I0)') '     block ', obj_rot%stationary_face(1,k), &
      !                              '  face-dir ', obj_rot%stationary_face(2,k)
      !    end do
      !  else
      !    write(*,'(A)') '   all type-5/6 walls rotate with the frame'
      !  end if
      !end if

      write(*,'(A)') ' ========================================================================================='

      ! ----- Numerical scheme -----
      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Numerical scheme'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A,T35,A)') ' Space'
      write(*,'(A,T35,A)') '   Reconstruction', trim(obj_space_scheme%description)
      if (len_trim(obj_space_scheme%limiter_name)>0) &
        write(*,'(A,T35,A)') '   Flux limiter',trim(obj_space_scheme%limiter_name)
      write(*,*)
      write(*,'(A,T35,A)') ' Time'
      write(*,'(A,T35,A)') '   Scheme', trim(obj_time_scheme%description)
      write(*,'(A,T35,A)') '   Integration variables', trim(obj_time_scheme%integration_name)
      if (len_trim(obj_irs%description)>0) &
        write(*,'(A,T35,A)') '   Implicit residual smoothing', trim(obj_irs%description)
      if ( obj_prec%enabled ) then
        write(*,'(A,T35,A)') '   Preconditioning', 'Weiss-Smith'
        if ( obj_prec%Uref < 0.0d0 ) then
          write(*,'(A,T35,A)') '     Uref', 'sound speed (from table)'
        else
          write(*,'(A,T35,ES12.4)') '     Uref', obj_prec%Uref
        end if
        if ( obj_prec%emin < 0.0d0 ) then
          write(*,'(A,T35,A)') '     eps_min', 'default (0.10)'
        else
          write(*,'(A,T35,ES12.4)') '     eps_min', obj_prec%emin
        end if
      end if
      write(*,*)
      write(*,'(A,T35,A)') ' Fluxes'
      write(*,'(A,T35,A)') '   Riemann solver', trim(obj_riemann%name)
      write(*,'(A)') ' ========================================================================================='
      write(*,*)

    end subroutine Print_Shell_Info

  end subroutine ARES_setup

end module ARES_Wrap_Setup