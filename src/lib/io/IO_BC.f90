module ARES_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! BC type counters — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, nio, nsym, nper, nchimera, ncoupled, ncoupled_ks, next

  ! Empty fields of a property line. A null value of list-directed input (two commas with only blanks between) leaves
  ! its variable as it was: undefined for a field without a default value. Read_BCfile sets those fields to EMPTY
  ! (IEMPTY for integers, a blank token for a direction) before the read, and Need_R, Need_I and Need_Inlet stop the run
  ! on a field still EMPTY after it; a field with a default value keeps it (docs/user/boundary-conditions.md). Every
  ! number of a property line, required or with a default value, must be finite (Need_Finite, before the test of EMPTY):
  ! every comparison with a NaN is false, so a NaN passes the test of EMPTY, and with -ffast-math the compiler can take
  ! it for EMPTY.
  real(R8), parameter :: EMPTY  = -huge(1.0_R8)
  integer,  parameter :: IEMPTY = -huge(1)

  ! The name of the boundary file that Read_BCfile reads, for the messages that name a line of it: bc.txt, or
  ! bc<level>.txt on a coarse multigrid level, after the phase prefix
  character(len=32) :: bc_name = 'bc.txt'

contains

  subroutine Setup_BC ( domain )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_multigrid, obj_io_bc
    use ARES_Global_m, only: model
    implicit none
    type(ARES_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    ! Local
    integer :: m, error

    do m = 1, obj_multigrid%MGL

      !! Phase 1: Allocate and check
      call Allocate_BC ( domain(m) )

      if (allocated(obj_io_bc%viscous_flag)) deallocate(obj_io_bc%viscous_flag)
      allocate ( obj_io_bc%viscous_flag( 6 * domain(m) % nb, 6 ) )
      obj_io_bc%viscous_flag = .false.

      if (allocated(obj_io_bc%coupling_flag)) deallocate(obj_io_bc%coupling_flag)
      allocate ( obj_io_bc%coupling_flag( 6 * domain(m) % nb, 6 ) )
      obj_io_bc%coupling_flag = .false.

      error = Check_BC ( domain(m) % nbound, m )
      if (error /= 0) cycle

      !! Phase 2: Read
      call Read_BCfile ( domain(m) % bc, domain(m) % n_bf, m )

      ! Every record on a cell of its block, next to the face it names, one record per boundary face (before the maps
      ! of the pairing and of type 200, which index with the cell of each record)
      if ( index(obj_io_bc%error_message,'ERROR') == 0 ) call Check_BC_Records ( domain(m) )

      call Pair_Connections ( domain(m), m )

      ! A face of type 200 takes what lies beyond it from the other k face of its block
      call Check_Rotational_Faces ( domain(m) )

      !if (m == 1) call Setup_Q2D_BC_Data ( domain(m) )

      if ((any(obj_io_bc%viscous_flag)) .and. (model==0)) then
        obj_io_bc%error_message = '[ERROR] Viscous wall BCs assigned with Eulerian model'
      end if

    end do

  end subroutine Setup_BC


  subroutine Allocate_BC ( domain )
    use ARES_Advanced_Types_m
    implicit none
    class(ARES_domain_type) :: domain
    ! Local
    integer :: b, n, ni, nj, nk

    n = 0
    do b = 1, domain % nb
      ni = domain % blk(b) % dim(1)
      nj = domain % blk(b) % dim(2)
      nk = domain % blk(b) % dim(3)
      n = n + 2*nj*nk + 2*ni*nk + 2*nj*ni
    enddo

    domain % nbound = n
    allocate ( domain % bc ( domain % nbound ) )
    allocate ( domain % n_bf ( domain % nb, 6 ) )
    
  end subroutine Allocate_BC


  function Check_BC (n, level) result(ios)
    use ARES_Config_Types_m, only: obj_io_bc
    use ARES_Global_m,       only: ARES_phase_prefix
    use IR_Precision,        only: str
    implicit none
    integer, intent(in)  :: n, level
    integer              :: ios
    integer              :: di(5), ti, unitfile, n_proof
    integer              :: c, ci, cii

    ios = 0

    ! Open file
    if (level == 1) then
      open(newunit=unitfile,file='INPUT/'//trim(ARES_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(ARES_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt',status='old',iostat=ios,action='read')
    endif
    if (ios/=0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    endif

    ! Cheak BC file consistency
    ios = 0; n_proof = -1
    do while (ios==0)
      read( unitfile,*,iostat=ios ) di(1), di(2), di(3), di(4), di(5), ti
      select case(ti)
      case(101, 103, 201, 301:302, 404:406, 501)
        read( unitfile,*,iostat=ios )
      case(102, 104) ! chimera, coupled chimera: donor counts of the two ghost rows, then one line per donor
        read( unitfile,*,iostat=ios ) ci, cii
        do c = 1, ci+cii
          read( unitfile,*,iostat=ios )
        enddo
      case(200, 300, 400) ! no property line
      case default
        ! A type ARES does not have: its property line, if any, would be counted as a record. At the end of the file
        ! (ios /= 0) ti is the type of the last record, already counted
        if ( ios == 0 ) then
          obj_io_bc%error_message = '[ERROR] Boundary file: record '//trim(str(.true.,n_proof+2))//' has type '// &
            trim(str(.true.,ti))//', which is not a boundary type of ARES (101 102 103 104 200 201 300 301 302 400 '// &
            '404 405 406 501)'
          close(unitfile)
          ios = 1
          return
        end if
      end select
      n_proof = n_proof + 1
    enddo

    if (n_proof /= n) then
      obj_io_bc%error_message = '[ERROR] Boundary conditions number ('//str(.true.,n_proof)//') is different than the one of the initial conditions ('//str(.true.,n)//')'
      close(unitfile)
      return
    endif

    ! Validation passed: reset ios (non-zero from EOF) to signal success
    ios = 0
    close(unitfile)

  end function Check_BC


  subroutine Read_BCfile ( bc, n_bf, level )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_io, obj_io_bc, obj_rans
    use ARES_Global_m
    use ARES_Mod_MPI,        only: mpi_abort_all, mpi_stop_all
    use ARES_Load_ThermoTransport, only: Table_Bounds, Table_Temperatures
    use IR_Precision
    implicit none
    type(ARES_bc_type), dimension(:), intent(inout) :: bc
    integer, intent(in)                             :: level
    integer, dimension(1:,1:), intent(inout)        :: n_bf
    ! Local
    integer :: cc, i, s
    integer :: unitfile, ios, cios, ip, ios_ks, idum(9)
    integer :: nline                         ! lines of the file read so far (Next_Line)
    integer :: nzero_chim
    real(R8) :: vf1, vf2
    real(R8) :: Tlo, Thi, plo, phi, pq, Tq_lo, Tq_hi   ! ranges of the fluid tables (Need_In_Tables)
    logical  :: inside, ok
    character(len=256) :: msg
    character(len=32) :: p0file
    character(len=32) :: alpha_tok, beta_tok
    character(len=:), allocatable :: pline   ! a line of the file, read whole

    cios = 0
    
    ! Open file
    if (level == 1) then
      open(newunit=unitfile,file='INPUT/'//trim(ARES_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(ARES_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt',status='old',iostat=ios,action='read')
    endif
    if (ios/=0) then
      obj_io_bc%error_message = '[ERROR] Boundary condition file not found for grid '//trim(str(.true.,level))
      return
    endif
    if (level == 1) then
      bc_name = trim(ARES_phase_prefix)//'bc.txt'
    else
      bc_name = trim(ARES_phase_prefix)//'bc'//trim(str(.true.,level))//'.txt'
    endif

    ! Counters for specific BC types
    if (level == 1) then
      nconnect = 0
      nwall = 0
      nio = 0
      nsym = 0
      nper = 0
      nchimera = 0
      ncoupled = 0
      ncoupled_ks = 0
      next = 0
    endif

    ! Counter for number of cells per face in each block
    n_bf = 0
    nzero_chim = 0

    ! The temperatures and the pressures of the fluid tables (loaded before the boundary file), against which the
    ! values of the property lines are checked
    call Table_Temperatures ( Tlo, Thi, plo, phi )

    ! Read file
    nline = 0
    do i = 1, size(bc)
       
      ! ── First line: block, ijk, face, ATLAS BC ID ─────────────────────────
      ! Every line is read whole (Next_Line, which counts the lines of the file and skips
      ! the blank ones as a list-directed read does), then as its values.
      call Next_Line ( unitfile, pline, ios, nline )
      if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%type
      if ( ios /= 0 ) call BC_Line_Error ( i, -1, nline, pline )

      ! Block and face within the mesh, before they index anything (every rank reads every record: collective stop)
      if ( bc(i)%b < 1 .or. bc(i)%b > size(n_bf,1) .or. bc(i)%f < 1 .or. bc(i)%f > 6 ) then
        write(msg,'(3A,I0,A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', i, &
          '): block ', bc(i)%b, ' or face ', bc(i)%f, ' out of range (blocks 1 to ', size(n_bf,1), ', faces 1 to 6)'
        call mpi_stop_all ( trim(msg) )
      end if

      ! n_bf update
      n_bf( bc(i) % b, bc(i) % f ) = n_bf( bc(i) % b, bc(i) % f ) + 1

      ! Second line is BC type-dependent
      select case( bc(i)%type )

        ! ─────────────────────────────────────────────────────────────────────
        ! Connection and periodic BCs
        case(101, 201)
          if (level == 1) nconnect = nconnect + 1
          bc(i)%bs = IEMPTY ; bc(i)%is = IEMPTY ; bc(i)%js = IEMPTY ; bc(i)%ks = IEMPTY ; bc(i)%fs = IEMPTY
          bc(i)%d11 = IEMPTY ; bc(i)%d12 = IEMPTY ; bc(i)%d21 = IEMPTY ; bc(i)%d22 = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, [ bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
                                        bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22 ], 'bs is js ks fs d11 d12 d21 d22' )
          allocate ( bc(i) % Pg (nprim, 6) )
          allocate ( bc(i) % Fc (nprim, 2) )
          bc(i) % Fc = 0d0

        ! ─────────────────────────────────────────────────────────────────────
        ! Chimera overlap BC (102) and coupled multi-solver chimera (104)
        ! Second line: n1, n2 = number of donors of the 1st and 2nd ghost cell
        ! [, roughness_ks] (104 only, optional and read as for 103),
        ! then n1+n2 lines: donor block, i, j, k, volume fraction.
        ! The 104 donors are FUSS cells: hydra-AF only uses them to find the
        ! FUSS faces the 104 face is coupled with, nothing is blended from them
        case(102, 104)
          if (level == 1) nchimera = nchimera + 1
          bc(i)%ni = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) (bc(i)%ni(cc),cc=1,2)
          if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
          call Need_I ( i, bc(i)%type, nline, bc(i)%ni, 'n1 n2' )
          allocate(bc(i)%donorID(1:sum(bc(i)%ni),1:4))
          allocate(bc(i)%volume_fraction(1:sum(bc(i)%ni)))
          bc(i)%donorID = IEMPTY
          bc(i)%volume_fraction = EMPTY
          if (bc(i)%type == 104) then
            read( pline,*,iostat=ios_ks ) idum(1:2), bc(i)%k_rough
          endif
          do s = 1, sum(bc(i)%ni)
            call Next_Line ( unitfile, pline, ios, nline )
            if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%donorID(s,1:4), bc(i)%volume_fraction(s)
            if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )
            call Need_I ( i, bc(i)%type, nline, bc(i)%donorID(s,1:4), 'bD iD jD kD' )
            call Need_Finite ( i, bc(i)%type, nline, bc(i)%volume_fraction(s:s), 'vol_frac' )
            call Need_R ( i, bc(i)%type, nline, bc(i)%volume_fraction(s:s), 'vol_frac' )
          enddo
          ! The ghost cells take the weighted mean of the donor primitives, so the
          ! fractions of each row must sum to 1: ATLAS writes them normalized, the
          ! division only removes round-off (a row without donors sums to 0)
          vf1 = sum(bc(i)%volume_fraction(1:bc(i)%ni(1)))
          vf2 = sum(bc(i)%volume_fraction(bc(i)%ni(1)+1:sum(bc(i)%ni)))
          if (vf1 < 0.5d0 .or. vf2 < 0.5d0) then
            nzero_chim = nzero_chim + 1
          else
            bc(i)%volume_fraction(1:bc(i)%ni(1)) = bc(i)%volume_fraction(1:bc(i)%ni(1)) / vf1
            bc(i)%volume_fraction(bc(i)%ni(1)+1:sum(bc(i)%ni)) = bc(i)%volume_fraction(bc(i)%ni(1)+1:sum(bc(i)%ni)) / vf2
          endif
          if (bc(i)%type == 102) then
            allocate ( bc(i) % Pg (nprim, 6) )
          else
            ! 104: coupled wall as for 103 (solid T stencil, coupling flux, wall output),
            ! counted with the 103 walls and with their roughness
            if (level == 1) ncoupled = ncoupled + 1
            if (ios_ks == 0) then
              if (level == 1) ncoupled_ks = ncoupled_ks + 1
            else
              bc(i)%k_rough = 0.0_R8 ! smooth
            endif
            obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
            allocate ( bc(i) % Pg (1, 6) )
            allocate ( bc(i) % ext_flux (nprim) )
            bc(i) % ext_flux = 0.0
          endif

        ! ─────────────────────────────────────────────────────────────────────
        ! Coupled multi-solver wall
        ! Second line: bs, is, js, ks, fs, d11, d12, d21, d22 [, roughness_ks]
        ! The trailing roughness is optional (older files have none): the line
        ! is read whole, so that a missing value is not taken from the next record.
        case(103)
          if (level == 1) ncoupled = ncoupled + 1
          obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%bs = IEMPTY ; bc(i)%is = IEMPTY ; bc(i)%js = IEMPTY ; bc(i)%ks = IEMPTY ; bc(i)%fs = IEMPTY
          bc(i)%d11 = IEMPTY ; bc(i)%d12 = IEMPTY ; bc(i)%d21 = IEMPTY ; bc(i)%d22 = IEMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, [ bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, &
                                        bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22 ], 'bs is js ks fs d11 d12 d21 d22' )
          if ( ios == 0 ) then
            read( pline,*,iostat=ios_ks ) idum, bc(i)%k_rough
            if (ios_ks == 0) then
              if (level == 1) ncoupled_ks = ncoupled_ks + 1
            else
              bc(i)%k_rough = 0.0_R8 ! smooth
            endif
          endif
          allocate(bc(i)%ext_flux(nprim))
          bc(i)%ext_flux = 0.0
          allocate ( bc(i) % Pg (1, 6) )

        ! ─────────────────────────────────────────────────────────────────────
        ! Axisymmetry
        case(200)
          if (level == 1) nper = nper + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Euler symmetry
        case(300)
          if (level == 1) nsym = nsym + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed heat flux
        ! Second line: q, roughness_ks, emissivity_eps  (comma-separated reals)
        case(301)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%qw = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%qw, bc(i)%k_rough, bc(i)%eps_wall
          if ( ios == 0 ) call Need_Finite ( i, bc(i)%type, nline, [ bc(i)%qw, bc(i)%k_rough, bc(i)%eps_wall ], &
                                             'q_wall k_rough eps_wall' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%qw ], 'q_wall' )

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed temperature
        ! Second line: T, roughness_ks, emissivity_eps
        case(302)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          bc(i)%Tw = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%Tw, bc(i)%k_rough, bc(i)%eps_wall
          if ( ios == 0 ) call Need_Finite ( i, bc(i)%type, nline, [ bc(i)%Tw, bc(i)%k_rough, bc(i)%eps_wall ], &
                                             'T_wall k_rough eps_wall' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%Tw ], 'T_wall' )
          if ( ios == 0 ) call Need_In_Tables ( i, bc(i)%type, nline, bc(i)%Tw, 'T_wall', Tlo, Thi, 'temperatures' )
          
        ! ─────────────────────────────────────────────────────────────────────
        ! Extrapolation
        case(400)
          if (level == 1) next = next + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, mass-flux g + T (static)
        ! Second line: T, g, alpha, beta, rel_fac, turb
        ! turb = the nrans turbulence values of the boundary, read by position as the specific
        ! ones ATLAS BCB writes (nuTilde, k, omega, Rij): the inlet weights them with the density
        ! of the boundary state (see Lib_BC_Fluxes_Inflow)
        case(404)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : np+nrans) )
          bc(i)%T0 = EMPTY ; bc(i)%mdot = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci(1:nrans) = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%T0, bc(i)%mdot ], 'T g', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )
          if ( ios == 0 ) call Need_Inward ( i, bc(i)%type, nline, bc(i)%f, bc(i)%mdot, 'g' )
          if ( ios == 0 ) call Need_In_Tables ( i, bc(i)%type, nline, bc(i)%T0, 'T', Tlo, Thi, 'temperatures' )

        ! ─────────────────────────────────────────────────────────────────────
        ! Supersonic inlet, M + T (static) + p (static)
        ! Second line: mach, T, p, alpha, beta, rel_fac, turb (specific, as for BC 404)
        case(405)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : np+nrans) )
          bc(i)%mach = EMPTY ; bc(i)%T0 = EMPTY ; bc(i)%pamb = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%ci(1:nrans) = EMPTY
          alpha_tok = ' ' ; beta_tok = ' '
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) &
            bc(i)%mach, bc(i)%T0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = 1, nrans)
          if ( ios == 0 ) call Need_Inlet ( i, nline, bc(i), [ bc(i)%mach, bc(i)%T0, bc(i)%pamb ], 'mach T p', &
                                            alpha_tok, beta_tok )
          if ( ios == 0 ) call Need_Direction ( i, bc(i)%type, nline, alpha_tok, beta_tok, bc(i)%alpha, bc(i)%beta )
          if ( ios == 0 ) then
            ! The pressure within the tables, then the temperature within those of the tables at that pressure
            ! (Table_Bounds brings a temperature far below and far above to the two ends of that range)
            call Need_In_Tables ( i, bc(i)%type, nline, bc(i)%pamb, 'p', plo, phi, 'pressures' )
            pq = bc(i)%pamb ; Tq_lo = -huge(1.0_R8)
            call Table_Bounds ( pq, Tq_lo, inside, ok )
            pq = bc(i)%pamb ; Tq_hi = huge(1.0_R8)
            call Table_Bounds ( pq, Tq_hi, inside, ok )
            call Need_In_Tables ( i, bc(i)%type, nline, bc(i)%T0, 'T', Tq_lo, Tq_hi, &
                                  'temperatures (at the p of the line)' )
          end if

        ! ─────────────────────────────────────────────────────────────────────
        ! Outlet, prescribed back pressure (if zero, extrapolation)
        ! Second line: p_back, rel_fac
        case(406)
          if (level == 1) nio = nio + 1
          bc(i)%pamb = EMPTY ; bc(i)%rel_fac = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%pamb, bc(i)%rel_fac
          if ( ios == 0 ) call Need_Finite ( i, bc(i)%type, nline, [ bc(i)%pamb, bc(i)%rel_fac ], 'p_back rel_fac' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%pamb, bc(i)%rel_fac ], 'p_back rel_fac' )
          ! p_back = 0: the pressure is extrapolated from the cell
          if ( ios == 0 .and. bc(i)%pamb /= 0d0 ) &
            call Need_In_Tables ( i, bc(i)%type, nline, bc(i)%pamb, 'p_back', plo, phi, 'pressures' )

        ! ─────────────────────────────────────────────────────────────────────
        ! Manifold (501): the faces of a group joined by a lumped plenum (Lib_BC_Fluxes_Manifold)
        ! Second line: gid, role, K, rel_fac, mdot_ratio, p0_ratio, dh0, as the boundary-condition builder of ATLAS writes
        !   them: the group id and the role of the face (1 outlet: the gas leaves the domain there; 2 inlet: it comes
        !   back), two integers; the loss K on the dynamic pressure of the outlet faces; the relaxation factor of the
        !   faces; the mass flow, total pressure and total enthalpy of the inlet faces against those of the outlet faces.
        !   The groups (the same parameters on every record, an outlet and an inlet face in each) are checked once every
        !   record of the level is read (Manifold_Setup).
        case(501)
          if (level == 1) nio = nio + 1
          bc(i)%man_gid = IEMPTY ; bc(i)%man_role = IEMPTY
          bc(i)%man_K = EMPTY ; bc(i)%rel_fac = EMPTY ; bc(i)%man_mratio = EMPTY ; bc(i)%man_p0ratio = EMPTY
          bc(i)%man_dh0 = EMPTY
          call Next_Line ( unitfile, pline, ios, nline )
          if ( ios == 0 ) read( pline,*,iostat=ios ) bc(i)%man_gid, bc(i)%man_role, bc(i)%man_K, bc(i)%rel_fac, &
                                                       bc(i)%man_mratio, bc(i)%man_p0ratio, bc(i)%man_dh0
          if ( ios == 0 ) call Need_I ( i, bc(i)%type, nline, [ bc(i)%man_gid, bc(i)%man_role ], 'gid role' )
          if ( ios == 0 ) call Need_Finite ( i, bc(i)%type, nline, [ bc(i)%man_K, bc(i)%rel_fac, bc(i)%man_mratio, &
                                             bc(i)%man_p0ratio, bc(i)%man_dh0 ], 'K rel_fac mdot_ratio p0_ratio dh0' )
          if ( ios == 0 ) call Need_R ( i, bc(i)%type, nline, [ bc(i)%man_K, bc(i)%rel_fac, bc(i)%man_mratio, &
                                        bc(i)%man_p0ratio, bc(i)%man_dh0 ], 'K rel_fac mdot_ratio p0_ratio dh0' )
          if ( ios == 0 ) call Need_Manifold ( i, nline, bc(i) )

        ! ─────────────────────────────────────────────────────────────────────
        ! A type ARES does not have: without a property line it would run with no flux on the face (Check_BC stops the
        ! types with one first, on the count of the records)
        case default
          write(msg,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', i, '): type ', &
            bc(i)%type, ' is not a boundary type of ARES (101 102 103 104 200 201 300 301 302 400 404 405 406 501)'
          call mpi_stop_all ( trim(msg) )

      end select

      ! A property line that cannot be read as the values of its type stops the run
      if ( ios /= 0 ) call BC_Line_Error ( i, bc(i)%type, nline, pline )

    enddo
       
    if (nwall > 0) obj_io % write_wall = .true.

    close( unitfile )

    if (nzero_chim > 0) then
      write(msg,'(I0,A,I0,A)') nzero_chim, ' chimera (102/104) entries in the level-', level, &
        ' BC file have donor volume fractions summing to ~0 (no donors found by the BC '// &
        'builder). Their ghost states (102) or coupling partners (104) would be garbage; '// &
        'fix the BC file before running.'
      call mpi_abort_all(trim(msg))
    endif

  end subroutine Read_BCfile


  subroutine Print_BC_Summary ()
    use ARES_Config_Types_m, only: obj_rans
    implicit none

    write(*,*)
    write(*,'(A)') ' Boundary conditions'
    if (nconnect > 0) write(*,'(A,T35,I0)') '   Connection', nconnect
    if (nwall > 0) write(*,'(A,T35,I0)') '   Viscous wall', nwall
    if (nio > 0) write(*,'(A,T35,I0)') '   Inflow/outflow', nio
    if (nsym > 0) write(*,'(A,T35,I0)') '   Symmetry', nsym
    if (nper > 0) write(*,'(A,T35,I0)') '   Periodicity', nper
    if (next > 0) write(*,'(A,T35,I0)') '   Extrapolation', next
    if (nchimera > 0) write(*,'(A,T35,I0)') '   Chimera', nchimera
    if (ncoupled > 0) write(*,'(A,T35,I0)') '   Coupled wall', ncoupled
    if (ncoupled_ks > 0) write(*,'(A,T35,I0)') '   Coupled wall with roughness', ncoupled_ks
    ! Coupled walls (103, 104) without the trailing roughness in bc.txt enter the
    ! nearest-wall roughness map with ks = 0, i.e. as smooth walls.
    if (ncoupled_ks < ncoupled .and. obj_rans%rough) &
      write(*,'(A,I0,A)') '   [WARNING] ', ncoupled - ncoupled_ks, &
        ' coupled walls carry no roughness: treated as smooth (ks = 0) by the rough-wall model'

  end subroutine Print_BC_Summary


  !> Pairs the two records of each connection face (101, 201): the record of the
  !> same face seen from the other block (partner), and which of the two computes
  !> the flux of the face for both cells (owner, BC_Fluxes). The owner is the
  !> record whose face is the high-index one of its cell (2, 4, 6) when the
  !> other is a low-index one, as on an interior face, where the flux between
  !> cells n and n+1 is computed once. When the two faces have the same parity
  !> (blocks that meet with their index directions permuted), the owner is the
  !> record with the smaller (block, face, k, j, i): a choice that does not
  !> depend on the order of the records in bc.txt. A record whose partner does
  !> not point back at it, or that points at itself, is an error.
  !> The search goes through a map, per block face, from the two indices of the
  !> cell along the face to its connection record: one pass over the records,
  !> in any order.
  !> The gc ghost cells of a connection face are the first gc cells of the
  !> block across it (Ghost_Connection), and their metrics come from its first
  !> gc node planes past the face (BC_Connect_Metrics). That block must have at
  !> least gc cells across its face, at every multigrid level (level): with
  !> fewer, a ghost cell would be a ghost cell of that block, which another
  !> record fills at the same time, and a node plane would lie outside the
  !> block. A record whose block across has fewer cells is an error; the
  !> directions along the face, and those of blocks with no connection, need
  !> nothing.
  subroutine Pair_Connections ( domain, level )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_io_bc
    use ARES_Global_m,       only: gc
    use IR_Precision, only: str
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    integer, intent(in)                   :: level
    ! Local
    type :: face_records_type
      integer, allocatable :: r(:,:)
    end type face_records_type
    type(face_records_type), allocatable :: rec(:,:)
    integer :: b, f, i, p, nt(2), t(2), n
    character(len=32) :: at_level

    allocate ( rec(domain % nb, 6) )

    ! Map of the connection records of each block face
    do i = 1, domain % nbound
      if ( domain % bc(i) % type /= 101 .and. domain % bc(i) % type /= 201 ) cycle
      b = domain % bc(i) % b
      f = domain % bc(i) % f
      if ( .not. allocated(rec(b,f) % r) ) then
        nt = Tangential ( f, domain % blk(b) % dim(1), domain % blk(b) % dim(2), domain % blk(b) % dim(3) )
        allocate ( rec(b,f) % r(nt(1), nt(2)) )
        rec(b,f) % r = 0
      end if
      t = Tangential ( f, domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k )
      rec(b,f) % r(t(1), t(2)) = i
    end do

    ! Partner and owner of each connection record
    do i = 1, domain % nbound
      if ( domain % bc(i) % type /= 101 .and. domain % bc(i) % type /= 201 ) cycle
      b = domain % bc(i) % bs
      f = domain % bc(i) % fs
      p = 0
      if ( b >= 1 .and. b <= domain % nb .and. f >= 1 .and. f <= 6 ) then
        if ( allocated(rec(b,f) % r) ) then
          t = Tangential ( f, domain % bc(i) % is, domain % bc(i) % js, domain % bc(i) % ks )
          if ( all(t >= 1) .and. all(t <= shape(rec(b,f) % r)) ) p = rec(b,f) % r(t(1), t(2))
        end if
      end if
      if ( p == 0 ) then
        obj_io_bc % error_message = '[ERROR] Connection record '//trim(str(.true.,i))// &
                                    ': no connection record at the cell and face it points at'
        exit
      end if
      if ( p == i ) then
        obj_io_bc % error_message = '[ERROR] Connection record '//trim(str(.true.,i))// &
                                    ' points at its own cell and face'
        exit
      end if
      if ( domain % bc(p) % type /= domain % bc(i) % type .or. domain % bc(p) % bs /= domain % bc(i) % b .or. &
           domain % bc(p) % fs /= domain % bc(i) % f .or. domain % bc(p) % is /= domain % bc(i) % i .or. &
           domain % bc(p) % js /= domain % bc(i) % j .or. domain % bc(p) % ks /= domain % bc(i) % k ) then
        obj_io_bc % error_message = '[ERROR] Connection records '//trim(str(.true.,i))//' and '//trim(str(.true.,p))// &
                                    ' do not point at each other'
        exit
      end if
      n = domain % blk(b) % dim((f + 1) / 2)
      if ( n < gc ) then
        at_level = ''
        if ( level > 1 ) at_level = ' at multigrid level '//trim(str(.true.,level))
        obj_io_bc % error_message = '[ERROR] Connection record '//trim(str(.true.,i))//' (block '// &
          trim(str(.true.,domain % bc(i) % b))//', face '//trim(str(.true.,domain % bc(i) % f))//'): block '// &
          trim(str(.true.,b))//' across it has '//trim(str(.true.,n))//' cell across its face '//trim(str(.true.,f))// &
          trim(at_level)//', fewer than the '//trim(str(.true.,gc))//' a connection needs'
        exit
      end if
      domain % bc(i) % partner = p
      if ( mod(domain % bc(i) % f, 2) == 0 .and. mod(domain % bc(p) % f, 2) == 1 ) then
        domain % bc(i) % owner = .true.
      else if ( mod(domain % bc(i) % f, 2) == 1 .and. mod(domain % bc(p) % f, 2) == 0 ) then
        domain % bc(i) % owner = .false.
      else
        domain % bc(i) % owner = Key_Less ( domain % bc(i) % b, domain % bc(i) % f, domain % bc(i) % k, &
                                            domain % bc(i) % j, domain % bc(i) % i, &
                                            domain % bc(p) % b, domain % bc(p) % f, domain % bc(p) % k, &
                                            domain % bc(p) % j, domain % bc(p) % i )
      end if
    end do

    contains

      !> Whether the record (b1, f1, k1, j1, i1) comes before (b2, f2, k2, j2, i2):
      !> block, then face, then k, j, i. Every record has its own key, so the two
      !> records of a face always tell which of them is first, whatever their
      !> order in bc.txt.
      pure logical function Key_Less ( b1, f1, k1, j1, i1, b2, f2, k2, j2, i2 )
        implicit none
        integer, intent(in) :: b1, f1, k1, j1, i1, b2, f2, k2, j2, i2
        integer :: key1(5), key2(5), n

        key1 = [ b1, f1, k1, j1, i1 ]
        key2 = [ b2, f2, k2, j2, i2 ]
        Key_Less = .false.
        do n = 1, 5
          if ( key1(n) /= key2(n) ) then
            Key_Less = key1(n) < key2(n)
            return
          end if
        end do

      end function Key_Less

      !> The two indices of a cell (or the two dimensions of a block) along the face f
      pure function Tangential ( f, i, j, k ) result ( t )
        implicit none
        integer, intent(in) :: f, i, j, k
        integer :: t(2)

        select case ( f )
          case (1:2)
            t = [ j, k ]
          case (3:4)
            t = [ i, k ]
          case default
            t = [ i, j ]
        end select

      end function Tangential

  end subroutine Pair_Connections


  !> Checks the records of type 200 against the rule of the type, which joins the two k faces of one block. Beyond a
  !> face of type 200 lies the cell next to the other k face of the same block, in the same column (i, j), turned by
  !> the angle of the wedge or of the sector: on a 3-D sector the ghost cell and the flux of the face take the cell
  !> k' = dim(3) - k + 1 of that block, and the record of face 6 computes the face once, for its cell and for that one
  !> (Ghost_Axysimmetry, BC_Rotational_Periodic_Eul); on a 2-D wedge, one cell in k, type 200 on faces 5 and
  !> 6 is what makes the mesh a wedge (docs/user/boundary-conditions.md). The record names no other block. Without
  !> this check a sector cut along k into blocks, a ring of several blocks closed by type 200 between the last block
  !> and the first, type 200 on one k face of a wedge or on a face 1 to 4 run with wrong fluxes, in silence (on a 3-D
  !> mesh a face 1 to 4 stops at the first flux, on the rank of its block alone). Every record of type
  !> 200 on face f of block b, cell (i, j, k), must lie on face 5 or 6 (R1), in a block of the mesh and inside it
  !> (R0), and have a record of type 200 on face 11 - f of block b in the same column (i, j) (R2); the first record
  !> of the file that does not writes the error, and the run stops at set-up on every rank (each rank reads every
  !> record: Check_BC and Stop_Together in Wrap_Setup). An error already written stays. The map holds the records of
  !> faces 5 and 6 of the blocks with type 200 there, of every type, inside their block: nothing is read or written
  !> outside an array. The numbers of the messages are written with I0, which keeps the sign of a negative index read
  !> from a wrong file.
  subroutine Check_Rotational_Faces ( domain )
    use ARES_Advanced_Types_m, only: ARES_domain_type
    use ARES_Config_Types_m,   only: obj_io_bc
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    ! Local
    type :: kface_records_type
      integer, allocatable :: r(:,:)
    end type kface_records_type
    type(kface_records_type), allocatable :: rec(:,:)
    logical, allocatable :: closed(:)
    integer :: n, b, f, i, j, k, p
    character(len=512) :: msg
    ! The remedy of a record without its partner: where the other k face of its column is a face between blocks (101,
    ! 201), the block goes on in k beyond it, a sector or a ring cut along k; elsewhere (another type, no record) the
    ! block holds the wedge or the sector, and its other k face takes type 200 too
    character(len=*), parameter :: rule = ': type 200 joins the k faces of one block: '
    character(len=*), parameter :: cut = 'cut sectors along i or j; join rings by 101 at each junction'

    if ( index(obj_io_bc % error_message, 'ERROR') > 0 ) return

    ! The blocks of the mesh with a record of type 200 on face 5 or 6
    allocate ( closed(domain % nb) )
    closed = .false.
    do n = 1, domain % nbound
      b = domain % bc(n) % b
      f = domain % bc(n) % f
      if ( domain % bc(n) % type == 200 .and. ( f == 5 .or. f == 6 ) .and. b >= 1 .and. b <= domain % nb ) &
        closed(b) = .true.
    end do

    ! Map of the records of faces 5 and 6 of those blocks, from the indices (i, j) of their cell
    allocate ( rec(domain % nb, 5:6) )
    do b = 1, domain % nb
      if ( .not. closed(b) ) cycle
      do f = 5, 6
        allocate ( rec(b,f) % r(domain % blk(b) % dim(1), domain % blk(b) % dim(2)) )
        rec(b,f) % r = 0
      end do
    end do
    do n = 1, domain % nbound
      b = domain % bc(n) % b
      f = domain % bc(n) % f
      if ( b < 1 .or. b > domain % nb ) cycle
      if ( .not. closed(b) .or. ( f /= 5 .and. f /= 6 ) ) cycle
      i = domain % bc(n) % i
      j = domain % bc(n) % j
      if ( i < 1 .or. i > domain % blk(b) % dim(1) .or. j < 1 .or. j > domain % blk(b) % dim(2) ) cycle
      rec(b,f) % r(i,j) = n
    end do

    ! Every record of type 200, in the order of the file: R1, R0, R2
    do n = 1, domain % nbound
      if ( domain % bc(n) % type /= 200 ) cycle
      b = domain % bc(n) % b
      f = domain % bc(n) % f
      i = domain % bc(n) % i
      j = domain % bc(n) % j
      k = domain % bc(n) % k
      msg = ''
      if ( f /= 5 .and. f /= 6 ) then
        write(msg,'(3A,I0,A,I0,A,I0,A,I0,1X,I0,1X,I0,A)') '[ERROR] ', trim(bc_name), ', record ', n, &
          ': type 200 on face ', f, ' of block ', b, ' at i, j, k = ', i, j, k, &
          ': type 200 is for the k faces (5 and 6) of a wedge or a sector only; declare symmetry (300) on the axis'
      else if ( b < 1 .or. b > domain % nb ) then
        ! Reached in builds without bound checks only: Read_BCfile counts the records of a block face (n_bf) first
        write(msg,'(3A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ', record ', n, ': type 200 on face ', f, &
          ' of block ', b, ', but the mesh has ', domain % nb, ' blocks'
      else if ( any ( [ i, j, k ] < 1 ) .or. any ( [ i, j, k ] > domain % blk(b) % dim ) ) then
        write(msg,'(3A,I0,A,I0,A,I0,A,I0,1X,I0,1X,I0,A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ', record ', n, &
          ': type 200 on face ', f, ' of block ', b, ' at i, j, k = ', i, j, k, ': the cell is outside the block, of ', &
          domain % blk(b) % dim(1), ' x ', domain % blk(b) % dim(2), ' x ', domain % blk(b) % dim(3), ' cells'
      else
        p = rec(b,11-f) % r(i,j)
        if ( p == 0 ) then
          write(msg,'(3A,I0,A,I0,A,I0,A,I0,1X,I0,A,I0,3A,I0,A)') '[ERROR] ', trim(bc_name), ', record ', n, &
            ': type 200 on face ', f, ' of block ', b, ' at i, j = ', i, j, ', but face ', 11 - f, ' there has no record', &
            rule, 'give face ', 11 - f, ' there type 200 too'
        else if ( domain % bc(p) % type == 101 .or. domain % bc(p) % type == 201 ) then
          write(msg,'(3A,I0,A,I0,A,I0,A,I0,1X,I0,A,I0,A,I0,2A)') '[ERROR] ', trim(bc_name), ', record ', n, &
            ': type 200 on face ', f, ' of block ', b, ' at i, j = ', i, j, ', but face ', 11 - f, ' there is of type ', &
            domain % bc(p) % type, rule, cut
        else if ( domain % bc(p) % type /= 200 ) then
          write(msg,'(3A,I0,A,I0,A,I0,A,I0,1X,I0,A,I0,A,I0,2A,I0,A)') '[ERROR] ', trim(bc_name), ', record ', n, &
            ': type 200 on face ', f, ' of block ', b, ' at i, j = ', i, j, ', but face ', 11 - f, ' there is of type ', &
            domain % bc(p) % type, rule, 'give face ', 11 - f, ' there type 200 too'
        end if
      end if
      if ( len_trim(msg) > 0 ) then
        obj_io_bc % error_message = msg
        exit
      end if
    end do

  end subroutine Check_Rotational_Faces

  !> Checks every record of the boundary file against the mesh: its cell lies inside its block and next to the face it
  !> names (face 1: i = 1, face 2: i = dim(1); 3 and 4: j; 5 and 6: k), and no two records name the same cell and face.
  !> Check_BC has found as many records as boundary faces, so a face named twice leaves another without a record, which
  !> would run with no flux and with ghost cells that nothing fills. Every rank reads every record: the first record
  !> that fails stops the run on every rank (mpi_stop_all), with its number and its cell.
  subroutine Check_BC_Records ( domain )
    use ARES_Advanced_Types_m, only: ARES_domain_type
    use ARES_Mod_MPI,          only: mpi_stop_all
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    ! Local
    type :: face_mask_type
      logical, allocatable :: m(:,:)
    end type face_mask_type
    type(face_mask_type), allocatable :: seen(:,:)
    integer :: n, b, f, d, c(3), dims(3), t(2)
    character(len=160) :: why

    allocate ( seen(domain % nb, 6) )
    do b = 1, domain % nb
      dims = domain % blk(b) % dim(1:3)
      do f = 1, 6
        t = Along ( f, dims )
        allocate ( seen(b,f) % m(t(1), t(2)) )
        seen(b,f) % m = .false.
      end do
    end do

    do n = 1, domain % nbound
      b = domain % bc(n) % b          ! within the mesh (Read_BCfile)
      f = domain % bc(n) % f
      c = [ domain % bc(n) % i, domain % bc(n) % j, domain % bc(n) % k ]
      dims = domain % blk(b) % dim(1:3)
      if ( any ( c < 1 ) .or. any ( c > dims ) ) then
        write(why,'(A,3(I0,A))') 'the cell is outside the block, of ', dims(1), ' x ', dims(2), ' x ', dims(3), ' cells'
        call Fail ( n, why )
      end if
      d = ( f + 1 ) / 2
      if ( c(d) /= merge ( 1, dims(d), mod(f, 2) == 1 ) ) &
        call Fail ( n, 'the cell is not next to the face it names (faces 1, 3, 5: index 1; faces 2, 4, 6: the last '// &
                       'cell of the block in that direction)' )
      t = Along ( f, c )
      if ( seen(b,f) % m(t(1), t(2)) ) &
        call Fail ( n, 'an earlier record names the same cell and face: each boundary face has one record' )
      seen(b,f) % m(t(1), t(2)) = .true.
    end do

    contains

      !> The two indices of a cell (or the two dimensions of a block) along the face f
      pure function Along ( f, c ) result ( t )
        implicit none
        integer, intent(in) :: f, c(3)
        integer :: t(2)

        select case ( f )
          case (1:2)
            t = [ c(2), c(3) ]
          case (3:4)
            t = [ c(1), c(3) ]
          case default
            t = [ c(1), c(2) ]
        end select

      end function Along

      subroutine Fail ( n, why )
        implicit none
        integer, intent(in)          :: n
        character(len=*), intent(in) :: why
        character(len=400) :: msg

        write(msg,'(3A,I0,A,I0,A,I0,A,3(1X,I0),A,I0,2A)') '[ERROR] ', trim(bc_name), ', record ', n, ' (block ', &
          domain % bc(n) % b, ', face ', domain % bc(n) % f, ', cell i, j, k =', domain % bc(n) % i, domain % bc(n) % j, &
          domain % bc(n) % k, ', type ', domain % bc(n) % type, '): ', trim(why)
        call mpi_stop_all ( trim(msg) )

      end subroutine Fail

  end subroutine Check_BC_Records


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_In_Tables: stop the run when the value v (named name) of record n, type btype, line nline of the file lies
  !   outside the range lo to hi of the fluid tables (what: the quantity of the range, 'temperatures' or 'pressures'):
  !   the state of the boundary would be read outside their arrays. v is finite (Need_Finite).
  subroutine Need_In_Tables ( n, btype, nline, v, name, lo, hi, what )
    use ARES_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    real(R8), intent(in)         :: v, lo, hi
    character(len=*), intent(in) :: name, what
    ! Local
    character(len=320) :: msg

    if ( v >= lo .and. v <= hi ) return
    write(msg,'(3A,I0,A,I0,A,I0,3A,ES12.5,3A,ES12.5,A,ES12.5,A)') '[ERROR] ', trim(bc_name), ' line ', nline, &
      ' (record ', n, ', type ', btype, '): ', name, ' = ', v, ' is outside the ', what, ' of the fluid tables (', &
      lo, ' to ', hi, ')'
    call mpi_stop_all ( trim(msg) )

  end subroutine Need_In_Tables


  !─────────────────────────────────────────────────────────────────────────────
  ! Next_Line: the next line of the boundary file that is not blank (Read_Line), as a list-directed read skips
  !   blank lines; nline counts every line read, blank ones included, so it is the line number in the file.
  subroutine Next_Line ( unit, line, ios, nline )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out)                       :: ios
    integer, intent(inout)                     :: nline

    do
      call Read_Line ( unit, line, ios )
      if ( ios /= 0 ) exit
      nline = nline + 1
      if ( len_trim(line) > 0 ) exit
    end do

  end subroutine Next_Line


  !─────────────────────────────────────────────────────────────────────────────
  ! BC_Line_Error: stop the run on line nline of the boundary file, of record n, that cannot be read as the values it
  !   must hold (short, or a value that is not a number): its record line (btype < 0) or a property line of type btype;
  !   an empty line is the end of the file.
  subroutine BC_Line_Error ( n, btype, nline, line )
    use ARES_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: line
    ! Local
    character(len=160) :: head

    if ( len_trim(line) == 0 ) then
      write(head,'(3A,I0,A,I0)') '[ERROR] ', trim(bc_name), ' ends at line ', nline, ', before the end of record ', n
      call mpi_stop_all ( trim(head) )
    end if
    if ( btype < 0 ) then
      write(head,'(3A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, &
        '): the line is not block, i, j, k, face and type'
    else
      write(head,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, &
        '): the line is short or holds a value that is not a number'
    end if
    call mpi_stop_all ( trim(head)//': "'//trim(line(1:min(len(line),200)))//'"' )

  end subroutine BC_Line_Error


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_R, Need_I: stop the run (Empty_Field) on the first of the values v of record n, type btype, line nline of the
  !   file, still EMPTY after the read of its property line. names: their names, blank-separated, in the order of the
  !   line; the last name stands for the values after it (a list of mass fractions, say).
  subroutine Need_R ( n, btype, nline, v, names )
    implicit none
    integer, intent(in)          :: n, btype, nline
    real(R8), intent(in)         :: v(:)
    character(len=*), intent(in) :: names
    ! Local
    integer :: k

    do k = 1, size(v)
      if ( v(k) == EMPTY ) call Empty_Field ( n, btype, nline, Name_Of ( names, k ) )
    end do

  end subroutine Need_R


  subroutine Need_I ( n, btype, nline, v, names )
    implicit none
    integer, intent(in)          :: n, btype, nline
    integer, intent(in)          :: v(:)
    character(len=*), intent(in) :: names
    ! Local
    integer :: k

    do k = 1, size(v)
      if ( v(k) == IEMPTY ) call Empty_Field ( n, btype, nline, Name_Of ( names, k ) )
    end do

  end subroutine Need_I


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_Inlet: the fields of an inlet line, in its order: the leading values lead (names), alpha, beta, rel_fac and the
  !   turbulence values: each one finite (Need_Finite; alpha and beta are read and checked by Need_Direction), then
  !   given (not empty).
  subroutine Need_Inlet ( n, nline, bc, lead, names, alpha_tok, beta_tok )
    use ARES_Advanced_Types_m, only: ARES_bc_type
    use ARES_Global_m,         only: nrans
    implicit none
    integer, intent(in)            :: n, nline
    type(ARES_bc_type), intent(in) :: bc
    real(R8), intent(in)           :: lead(:)
    character(len=*), intent(in)   :: names, alpha_tok, beta_tok

    call Need_Finite ( n, bc%type, nline, lead, names )
    call Need_Finite ( n, bc%type, nline, [ bc%rel_fac ], 'rel_fac' )
    call Need_Finite ( n, bc%type, nline, bc%ci(1:nrans), 'turb' )
    call Need_R ( n, bc%type, nline, lead, names )
    if ( len_trim(alpha_tok) == 0 ) call Empty_Field ( n, bc%type, nline, 'alpha' )
    if ( len_trim(beta_tok) == 0 ) call Empty_Field ( n, bc%type, nline, 'beta' )
    call Need_R ( n, bc%type, nline, [ bc%rel_fac ], 'rel_fac' )
    call Need_R ( n, bc%type, nline, bc%ci(1:nrans), 'turb' )

  end subroutine Need_Inlet

  !─────────────────────────────────────────────────────────────────────────────
  ! Need_Manifold: the values of a manifold record (501), record n at line nline: the group id from 1, the role 1 (outlet)
  !   or 2 (inlet), K >= 0, 0 < rel_fac <= 1, mdot_ratio > 0 and p0_ratio > 0 (dh0 takes any value). The groups
  !   themselves -- the same parameters on every record, an outlet and an inlet face in each -- are checked once every
  !   record of the level is read (Manifold_Setup, Lib_BC_Fluxes_Manifold).
  subroutine Need_Manifold ( n, nline, bc )
    use ARES_Advanced_Types_m, only: ARES_bc_type
    use ARES_Mod_MPI,          only: mpi_stop_all
    implicit none
    integer,            intent(in) :: n, nline
    type(ARES_bc_type), intent(in) :: bc
    ! Local
    character(len=400) :: msg
    character(len=120) :: head

    write(head,'(3A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type 501):'
    if ( bc%man_gid < 1 ) then
      write(msg,'(A,I0,A)') trim(head)//' the group id is ', bc%man_gid, '; the ids of the manifolds start at 1'
    else if ( bc%man_role /= 1 .and. bc%man_role /= 2 ) then
      write(msg,'(A,I0,A)') trim(head)//' the role is ', bc%man_role, &
        '; it is 1 (outlet: the gas leaves the domain through the face) or 2 (inlet: it comes back)'
    else if ( bc%man_K < 0d0 ) then
      write(msg,'(A,ES12.5,A)') trim(head)//' K = ', bc%man_K, &
        ' is negative; the loss on the dynamic pressure of the outlet faces is 0 or positive'
    else if ( bc%rel_fac <= 0d0 .or. bc%rel_fac > 1d0 ) then
      write(msg,'(A,ES12.5,A)') trim(head)//' rel_fac = ', bc%rel_fac, ' is outside (0, 1]'
    else if ( bc%man_mratio <= 0d0 ) then
      write(msg,'(A,ES12.5,A)') trim(head)//' mdot_ratio = ', bc%man_mratio, &
        ' is not positive; it is the mass flow of the inlet faces over that of the outlet faces'
    else if ( bc%man_p0ratio <= 0d0 ) then
      write(msg,'(A,ES12.5,A)') trim(head)//' p0_ratio = ', bc%man_p0ratio, &
        ' is not positive; it is the total pressure of the inlet faces over that of the outlet faces'
    else
      return
    end if
    call mpi_stop_all ( trim(msg) )

  end subroutine Need_Manifold


  !─────────────────────────────────────────────────────────────────────────────
  ! Need_Inward: the mass flux g of an inlet 404 (v, named name), of record n, type btype, line nline of the file, on
  !   face face of its block, already finite and given (Need_Inlet). v is positive into the domain on every face: an
  !   inlet only lets gas in, so zero, or a negative value on a face 1, 3 or 5, stops the run. On a face 2, 4 or 6 a
  !   negative value is taken as its magnitude: it is the component of a gas that enters along the face normal of the
  !   solver, which points toward growing index, out of the domain on those faces. The flux routine takes -BC_Sign |v|
  !   along that normal.
  subroutine Need_Inward ( n, btype, nline, face, v, name )
    use ARES_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline, face
    real(R8), intent(in)         :: v
    character(len=*), intent(in) :: name
    ! Local
    character(len=400) :: msg
    character(len=16)  :: what, val

    what = 'a mass flux'
    if ( v == 0d0 ) then
      write(msg,'(3A,I0,A,I0,A,I0,3A,I0,5A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', &
        btype, '): ', name, ' = 0 on face ', face, ': an inlet needs ', trim(what), ' into the domain (', name, &
        ' > 0 on every face)'
      call mpi_stop_all ( trim(msg) )
    else if ( v < 0d0 .and. mod(face, 2) == 1 ) then
      write(val,'(ES12.5)') v
      write(msg,'(3A,I0,A,I0,A,I0,5A,I0,3A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', &
        btype, '): ', name, ' = ', trim(adjustl(val)), ' on face ', face, &
        ' means gas leaving the domain: an inlet only lets gas in (', name, &
        ' > 0 into the domain on every face); give an outlet type (406, or 400) to a face where the gas leaves'
      call mpi_stop_all ( trim(msg) )
    end if

  end subroutine Need_Inward

  !─────────────────────────────────────────────────────────────────────────────
  ! Need_Finite: stop the run on the first of the values v of record n, type btype, line nline of the file that is not
  !   finite (NaN or infinity): the numbers of a property line, required or with a default value. names as for Need_R.
  !   The test reads the bits of each value (Not_Finite), and comes before the test of EMPTY (Need_R): every comparison
  !   with a NaN is false, so that test lets a NaN through, and with -ffast-math the compiler can take a NaN for EMPTY
  !   (observed with ifx); a build that traps the invalid operations (the DEBUG build) would stop with SIGFPE at the
  !   first operation on the value.
  subroutine Need_Finite ( n, btype, nline, v, names )
    use ARES_Finite_m, only: Not_Finite
    use ARES_Mod_MPI,  only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    real(R8), intent(in)         :: v(:)
    character(len=*), intent(in) :: names
    ! Local
    integer :: k
    character(len=240) :: msg

    do k = 1, size(v)
      if ( Not_Finite ( v(k) ) ) then
        write(msg,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, &
          '): the value of '//trim(Name_Of ( names, k ))//' is not finite (NaN or infinity)'
        call mpi_stop_all ( trim(msg) )
      end if
    end do

  end subroutine Need_Finite


  !─────────────────────────────────────────────────────────────────────────────
  ! Empty_Field: stop the run on the field what of the property line of record n, type btype, line nline of the file,
  !   left empty.
  subroutine Empty_Field ( n, btype, nline, what )
    use ARES_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: what
    ! Local
    character(len=240) :: msg

    write(msg,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, &
      '): the value of '//trim(what)//' is empty (docs/user/boundary-conditions.md)'
    call mpi_stop_all ( trim(msg) )

  end subroutine Empty_Field


  !─────────────────────────────────────────────────────────────────────────────
  ! Name_Of: the k-th blank-separated word of names, or the last one when there are fewer.
  pure function Name_Of ( names, k ) result ( w )
    implicit none
    character(len=*), intent(in) :: names
    integer, intent(in)          :: k
    character(len=len(names))    :: w
    ! Local
    integer :: p, q, m

    w = ' '
    p = 1
    m = 0
    do
      do while ( p <= len_trim(names) )
        if ( names(p:p) /= ' ' ) exit
        p = p + 1
      end do
      if ( p > len_trim(names) ) exit
      q = index ( names(p:), ' ' )
      if ( q == 0 ) then
        q = len(names) + 1
      else
        q = p + q - 1
      end if
      m = m + 1
      w = names(p:q-1)
      if ( m == k ) exit
      p = q
    end do

  end function Name_Of


  !─────────────────────────────────────────────────────────────────────────────
  ! Read_Line: the next line of a formatted file, whole, whatever its length (read in
  !   pieces without advancing). ios < 0 at the end of the file; a last line without its
  !   end-of-line mark is still a line.
  subroutine Read_Line ( unit, line, ios )
    implicit none
    integer, intent(in)                        :: unit
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out)                       :: ios
    ! Local
    character(len=256) :: piece
    integer :: nread

    line = ''
    do
      read ( unit, '(A)', advance='no', iostat=ios, size=nread ) piece
      line = line // piece(1:nread)
      if ( ios /= 0 ) exit
    end do
    if ( is_iostat_eor(ios) ) ios = 0
    if ( is_iostat_end(ios) .and. len(line) > 0 ) ios = 0

  end subroutine Read_Line

  !─────────────────────────────────────────────────────────────────────────────
  ! Need_Direction: the flow directions alpha and beta of an inlet line (record n, type btype, line nline of the file)
  !   from their tokens; a token that is neither a number nor normal, or a number that is not finite (Need_Finite), stops
  !   the run.
  subroutine Need_Direction ( n, btype, nline, alpha_tok, beta_tok, alpha, beta )
    use ARES_Mod_MPI, only: mpi_stop_all
    implicit none
    integer, intent(in)          :: n, btype, nline
    character(len=*), intent(in) :: alpha_tok, beta_tok
    real(R8), intent(out)        :: alpha, beta
    ! Local
    integer            :: ios
    character(len=200) :: msg

    call parse_dir_tok ( alpha_tok, alpha, ios )
    if ( ios /= 0 ) then
      write(msg,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, &
        "): alpha is neither a number nor normal: '"//trim(adjustl(alpha_tok))//"'"
      call mpi_stop_all ( trim(msg) )
    end if
    call parse_dir_tok ( beta_tok, beta, ios )
    if ( ios /= 0 ) then
      write(msg,'(3A,I0,A,I0,A,I0,A)') '[ERROR] ', trim(bc_name), ' line ', nline, ' (record ', n, ', type ', btype, &
        "): beta is neither a number nor normal: '"//trim(adjustl(beta_tok))//"'"
      call mpi_stop_all ( trim(msg) )
    end if
    call Need_Finite ( n, btype, nline, [ alpha, beta ], 'alpha beta' )

  end subroutine Need_Direction


  !─────────────────────────────────────────────────────────────────────────────
  ! parse_dir_tok: convert a direction token (read as character) to real.
  !   If the token contains 'normal', return huge (face-normal flag).
  !   Otherwise parse it as a floating-point number; ios /= 0 when it is not one (the caller then stops the run).
  pure subroutine parse_dir_tok ( tok, val, ios )
    implicit none
    character(len=*), intent(in) :: tok
    real(R8), intent(out)        :: val
    integer, intent(out)         :: ios
    character(len=len(tok)) :: tok_
    tok_ = adjustl(tok)
    ios = 0
    if (index(trim(tok_), 'normal') > 0) then
      val = huge(1.0_R8)
    else
      read(tok_, *, iostat=ios) val
      if (ios /= 0) val = 0.0_R8
    endif
  end subroutine parse_dir_tok
end module ARES_IO_BC
