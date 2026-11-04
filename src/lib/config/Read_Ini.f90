module ARES_Read_Ini
  
  implicit none

contains

  subroutine Read_Inifile()
    use Finer,               only: file_ini
    use ARES_Read_Sim_Param, only: Register_Sim_Param
    use ARES_Read_IO,        only: Register_IO_Fields, Register_Probes
    use ARES_Read_Numerics,  only: Register_Numerics
    use ARES_Read_Physics,   only: Register_Physics
    use ARES_Backend_INI,    only: Load_Ini, Scan_Ini
    use ARES_Input_Registry
    use ARES_Config_Types_m
    implicit none
    ! Local
    type(file_ini) :: fini
    integer :: nprobes, nmgl
    character(len=16), allocatable :: probes_name(:)

    ! Load input.ini: a file that cannot be read stops the run
    call Load_Input_File(fini, 'Input file')

    ! Scan input.ini for unknown number of probes and multigrid levels
    call Scan_Ini(fini, nprobes, probes_name, nmgl)

    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(nprobes, probes_name)
    call Register_Numerics(nmgl)
    call Register_Physics()

    ! Registry is built, now load the values from the ini file
    call Load_Ini(fini)

    ! The keys of input.ini that the reader does not give
    call Warn_Unread_Keys(fini)

  end subroutine Read_Inifile


  !> Prints a warning for every line "key = value" of input.ini whose key the INI reader (FiNeR) does not give in its
  !> section: that key is not read, and its parameter keeps its default. Comment lines (they start with !, ; or #) are
  !> skipped, as the reader skips them since FiNeR 2.1.0, also when they hold "=" or a bracket.
  subroutine Warn_Unread_Keys(fini)
    use Finer,        only: file_ini
    use ARES_Mod_MPI, only: mpi_is_root
    implicit none
    type(file_ini), intent(in) :: fini
    ! Local
    character(len=1024) :: line, section, key, val
    integer :: unit, ios, n, eq, error

    if (.not. mpi_is_root) return
    open(newunit=unit, file='input.ini', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    section = ''
    n = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      n = n + 1
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (scan(line(1:1), '!;#') == 1) cycle
      if (line(1:1) == '[') then
        if (index(line, ']') > 2) section = adjustl(line(2:index(line, ']')-1))
        cycle
      end if
      eq = index(line, '=')
      if (eq <= 1) cycle
      key = adjustl(line(1:eq-1))
      call fini%get(section_name=trim(section), option_name=trim(key), val=val, error=error)
      if (error /= 0) write(*,'(A,I0,A)') '[WARNING] input.ini, line ', n, ': '//trim(key)//' of ['//trim(section)// &
        '] is not read (the INI reader does not give it): its parameter keeps its default'
    end do
    close(unit)

  end subroutine Warn_Unread_Keys


  !> input.ini read again during the run, every ini-diter iterations (ARES_postprocess, at iteration iter). The values
  !> read again keep to the rules of the registry, as at set-up (Check_Input in Wrap_Setup). An error of the reading (a
  !> file that cannot be read, a value that does not read as its type, is out of its rule or is not finite) stops
  !> the run with ini-error = stop, as at set-up. With ini-error = warn, the default, the reading is not taken: a warning
  !> gives the iteration and the error (Write_Not_Taken), and the run goes on with every value it had before the reading,
  !> those of the keys of the same reading that were right too (Save_Values before Load_Ini, Restore_Values on the
  !> error). The value of ini-error in force when the reading starts decides, also when the file read again writes
  !> another one: Load_Ini writes the value of the file in its variable before the checks. Every rank checks the text
  !> that rank 0 read and passed to it (Load_Input_File), so every rank finds the same values; the ranks take each
  !> decision together (the largest of a flag), so that every rank stops, or every rank goes on with the same values (a
  !> rank that stopped alone would wait for the others in mpi_stop_all, and ranks with other values would call their
  !> collectives in other orders).
  subroutine Read_Inifile_Runtime(iter)
    use Finer,               only: file_ini
    use ARES_Backend_INI,    only: Load_Ini
    use ARES_Input_Registry
    use ARES_Config_Types_m, only: obj_io
    use ARES_Mod_MPI,        only: mpi_rank_, mpi_allreduce_max_r8, mpi_allreduce_min_r8, mpi_stop_all
    use iso_fortran_env,     only: R8 => real64
    implicit none
    integer, intent(in) :: iter   !< The iteration at which input.ini is read again
    ! Local
    type(file_ini) :: fini
    character(len=1024) :: out
    real(R8) :: flag, flag_all, first
    logical  :: warn, taken

    ! The value of ini-error in force when the reading starts
    warn = ( trim(obj_io%ini_error) == 'warn' )

    ! Load input.ini: a file that cannot be read stops the run, or with warn leaves the reading not taken
    call Load_Input_File(fini, 'Input file, read again', warn, iter, taken)
    if (.not. taken) return

    ! The values before the reading, which a reading that is not taken puts back; then the values of the file
    call Save_Values()
    call Load_Ini(fini)

    ! Validate registry: the values read again keep to the rules of the registry, as at setup. The rule broken is
    ! written once, by the first rank that found it (rank 0: every rank checks the same text)
    out = Validate_Registry()
    flag = 0.0_R8
    if (index(out,'ERROR') > 0) flag = 1.0_R8
    call mpi_allreduce_max_r8(flag, flag_all)
    if (flag_all > 0.0_R8) then
      call mpi_allreduce_min_r8(merge(real(mpi_rank_, R8), huge(1.0_R8), flag > 0.0_R8), first)
      if (warn) then
        if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) call Write_Not_Taken('Input file, read again', iter, out)
        call Restore_Values()
        return
      end if
      if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) then
        write(*,'(A,T35,A)') '   Input file, read again', 'FAIL'
        write(*,'(4X,A)') trim(out)
      end if
      if (flag > 0.0_R8) call mpi_stop_all('input file error, input.ini read again during the run (see above)')
      call mpi_stop_all('stopped with the ranks that found an error in input.ini read again')
    end if

  end subroutine Read_Inifile_Runtime


  !> Loads input.ini into fini, and stops the run if it cannot be read: a file that cannot be opened (missing, or without
  !> read permission), a file with no line that is not blank, a file whose read ends with an I/O error, or an error of
  !> the INI reader (FiNeR load gives an error code other than 0: a file with no section, a line that it cannot parse;
  !> FiNeR names such a line on stderr, without its number). Without this stop, the parameters would keep their defaults
  !> at set-up and the values of the previous reading during the run, and an empty file would end the run in the reader.
  !> During the run with ini-error = warn (warn present and true) such a file leaves the reading not taken instead: a
  !> warning (Write_Not_Taken) and taken false, so that the caller keeps every value it had (Read_Inifile_Runtime).
  !> Rank 0 alone reads the file, as FiNeR reads a file (string%read_file of StringiFor), and passes its text to every
  !> rank (mpi_bcast_chars); every rank, rank 0 too, gives that text to FiNeR (load from a source), so every rank takes
  !> the same values. If each rank read input.ini on its own, a file rewritten while the ranks read it during the run
  !> could give them different values, and a value that decides a collective call (sol-diter, res-diter, ini-diter)
  !> would make them call their collectives in different orders (forbidden by MPI-4.1, section 6.14), and hang. A file
  !> that rank 0 cannot read stops every rank, or leaves the reading not taken on every rank: rank 0 writes why and
  !> passes the outcome (mpi_bcast_integer), then every rank stops with mpi_stop_all (exit code 1), or returns; an error
  !> of the INI reader, the same on every rank, is decided together too (the largest of a flag), and the first rank that
  !> found it writes it.
  subroutine Load_Input_File(fini, label, warn, iter, taken)
    use Finer,           only: file_ini
    use stringifor,      only: string
    use ARES_Mod_MPI,    only: mpi_is_root, mpi_rank_, mpi_bcast_integer, mpi_bcast_chars, mpi_allreduce_max_r8, &
                               mpi_allreduce_min_r8, mpi_stop_all
    use iso_fortran_env, only: R8 => real64
    implicit none
    type(file_ini),    intent(inout) :: fini
    character(len=*),  intent(in)    :: label   !< Line of the set-up summary: Input file, or Input file, read again
    logical, optional, intent(in)    :: warn    !< A reading during the run with ini-error = warn: an error, not taken
    integer, optional, intent(in)    :: iter    !< The iteration of a reading during the run, for the warning
    logical, optional, intent(out)   :: taken   !< The file is loaded (false: a reading that is not taken, with warn)
    ! Local
    character(len=512)  :: reason, msg
    character(len=1024) :: line
    character(len=:), allocatable :: text      !< The text of input.ini that rank 0 reads
    type(string) :: source
    integer  :: unit, ios, error, n, state, it
    real(R8) :: flag, flag_all, first
    logical  :: not_fatal

    not_fatal = .false.
    if (present(warn)) not_fatal = warn
    it = 0
    if (present(iter)) it = iter
    if (present(taken)) taken = .true.

    ! Rank 0 reads input.ini. state: 0 read, 1 cannot be opened, 2 empty, 3 the read of the text failed
    reason = ''
    state = 0
    if (mpi_is_root) then
      open(newunit=unit, file='input.ini', status='old', action='read', iostat=ios, iomsg=msg)
      if (ios /= 0) then
        state = 1
        reason = 'cannot be opened: '//trim(msg)
      else
        n = 0
        do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          if (len_trim(line) > 0) n = n + 1
        end do
        close(unit)
        if (n == 0) then
          state = 2
          reason = 'is empty'
        end if
      end if
      if (state == 0) then
        ! As FiNeR reads a file (file_ini%load with a file name): read_file, then the parse of its text, below. An error
        ! of read_file is an I/O error of the read (StringiFor), not an error of the INI reader
        call source%read_file(file='input.ini', iostat=ios)
        if (ios > 0) then
          state = 3
          write(reason, '(A,I0)') 'cannot be read: I/O error ', ios
        else
          text = source%chars()
        end if
      end if
    end if
    call mpi_bcast_integer(state)
    if (state /= 0) then
      if (not_fatal) then
        if (mpi_is_root) call Write_Not_Taken(label, it, 'input.ini '//trim(reason))
        if (present(taken)) taken = .false.
        return
      end if
      if (mpi_is_root) then
        write(*,'(A,T35,A)') '   '//trim(label), 'FAIL'
        write(*,'(4X,A)') '[ERROR] input.ini '//trim(reason)
      end if
      call mpi_stop_all('input.ini cannot be read (see above)')
    end if

    ! Every rank parses the text that rank 0 read
    call mpi_bcast_chars(text)
    call fini%load(source=text, error=error)
    if (error /= 0) write(reason, '(A,I0,A)') 'cannot be read by the INI reader (FiNeR error ', error, &
      '; a line that FiNeR could not parse is named above, if any)'

    flag = 0.0_R8
    if (error /= 0) flag = 1.0_R8
    call mpi_allreduce_max_r8(flag, flag_all)
    if (flag_all > 0.0_R8) then
      call mpi_allreduce_min_r8(merge(real(mpi_rank_, R8), huge(1.0_R8), flag > 0.0_R8), first)
      if (not_fatal) then
        if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) call Write_Not_Taken(label, it, 'input.ini '//trim(reason))
        if (present(taken)) taken = .false.
        return
      end if
      if (flag > 0.0_R8 .and. nint(first) == mpi_rank_) then
        write(*,'(A,T35,A)') '   '//trim(label), 'FAIL'
        write(*,'(4X,A)') '[ERROR] input.ini '//trim(reason)
      end if
      call mpi_stop_all('input.ini cannot be read (see above)')
    end if

  end subroutine Load_Input_File


  !> The warning of a reading of input.ini during the run that is not taken (ini-error = warn): the line of the set-up
  !> summary with NOT TAKEN, the iteration and the error, without the [ERROR] in front of it (the run goes on), and that
  !> no value of the reading is taken. Written out at once (flush), as the message of a stop is (mpi_stop_all).
  subroutine Write_Not_Taken(label, iter, error)
    use iso_fortran_env, only: output_unit
    implicit none
    character(len=*), intent(in) :: label   !< Line of the set-up summary: Input file, read again
    integer,          intent(in) :: iter    !< The iteration of the reading
    character(len=*), intent(in) :: error   !< The error, as a stop gives it
    ! Local
    character(len=:), allocatable :: reason

    reason = trim(adjustl(error))
    if (index(reason, '[ERROR] ') == 1) reason = reason(9:)
    write(*,'(A,T35,A)') '   '//trim(label), 'NOT TAKEN'
    write(*,'(4X,A,I0,A)') '[WARNING] input.ini read again at iteration ', iter, ' is not taken: '//reason
    write(*,'(4X,A)') 'No value of this reading is taken: the run goes on with the values it had before the reading.'
    flush(output_unit)

  end subroutine Write_Not_Taken


end module ARES_Read_Ini