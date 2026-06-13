module ARES_IO_BC
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_BC
  public :: Print_BC_Summary

  ! BC type counters — populated by Read_BCfile, consumed by Print_BC_Summary
  integer :: nconnect, nwall, nio, nsym, nper, ncoupled, next
  logical :: has_tdep_bc

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
      case(101, 103, 201, 301:302, 404:406)
        read( unitfile,*,iostat=ios )
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
    use IR_Precision
    implicit none
    type(ARES_bc_type), dimension(:), intent(inout) :: bc
    integer, intent(in)                             :: level
    integer, dimension(1:,1:), intent(inout)        :: n_bf
    ! Local
    integer :: cc, i, s
    integer :: unitfile, ios, cios, ip
    character(len=32) :: p0file
    character(len=32) :: alpha_tok, beta_tok

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

    ! Counters for specific BC types
    if (level == 1) then
      nconnect = 0
      nwall = 0
      nio = 0
      nsym = 0
      nper = 0
      ncoupled = 0
      next = 0
      has_tdep_bc = .false.
    endif

    ! Counter for number of cells per face in each block
    n_bf = 0

    ! Read file
    do i = 1, size(bc)
       
      ! ── First line: block, ijk, face, ATLAS BC ID ─────────────────────────
      read( unitfile,*,iostat=ios ) bc(i)%b, bc(i)%i, bc(i)%j, bc(i)%k, bc(i)%f, bc(i)%type
      if (ios/=0) write(*,'(A)') '  Error in BC file'

      ! n_bf update
      n_bf( bc(i) % b, bc(i) % f ) = n_bf( bc(i) % b, bc(i) % f ) + 1

      ! Second line is BC type-dependent
      select case( bc(i)%type )

        ! ─────────────────────────────────────────────────────────────────────
        ! Connection and periodic BCs
        case(101, 201)
          if (level == 1) nconnect = nconnect + 1
          read( unitfile,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
          allocate ( bc(i) % Pg (nprim, 6) )

        ! ─────────────────────────────────────────────────────────────────────
        ! Coupled multi-solver wall
        case(103)
          if (level == 1) ncoupled = ncoupled + 1
          obj_io_bc%coupling_flag( bc(i)%b , bc(i)%f ) = .true.
          read( unitfile,*,iostat=ios ) &
            bc(i)%bs, bc(i)%is, bc(i)%js, bc(i)%ks, bc(i)%fs, bc(i)%d11, bc(i)%d12, bc(i)%d21, bc(i)%d22
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
          read(unitfile,*,iostat=ios) bc(i)%qw, bc(i)%k_rough, bc(i)%eps_wall

        ! ─────────────────────────────────────────────────────────────────────
        ! Wall, prescribed temperature
        ! Second line: T, roughness_ks, emissivity_eps
        case(302)
          if (level == 1) nwall = nwall + 1
          obj_io_bc%viscous_flag( bc(i)%b , bc(i)%f ) = .true.
          read(unitfile,*,iostat=ios) bc(i)%Tw, bc(i)%k_rough, bc(i)%eps_wall
          
        ! ─────────────────────────────────────────────────────────────────────
        ! Extrapolation
        case(400)
          if (level == 1) next = next + 1

        ! ─────────────────────────────────────────────────────────────────────
        ! Inlet, mass-flux g + T (static)
        ! Second line: T, g, alpha, beta, rel_fac, massf, turb
        case(404)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : np+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%T0, bc(i)%mdot, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Supersonic inlet, M + T (static) + p (static)
        ! Second line: mach, T, p, alpha, beta, rel_fac, massf, turb
        case(405)
          if (level == 1) nio = nio + 1
          allocate( bc(i) % ci(1 : np+nrans) )
          read( unitfile,*,iostat=ios ) &
            bc(i)%mach, bc(i)%T0, bc(i)%pamb, alpha_tok, beta_tok, bc(i)%rel_fac, (bc(i)%ci(s), s = 1, nrans)
          bc(i)%alpha = parse_dir_tok(alpha_tok)
          bc(i)%beta  = parse_dir_tok(beta_tok)

        ! ─────────────────────────────────────────────────────────────────────
        ! Outlet, prescribed back pressure (if zero, extrapolation)
        ! Second line: p_back, rel_fac
        case(406)
          if (level == 1) nio = nio + 1
          read( unitfile,*,iostat=ios ) bc(i)%pamb, bc(i)%rel_fac
          
      end select

    enddo
       
    if (nwall > 0) obj_io % write_wall = .true.

    close( unitfile )

  end subroutine Read_BCfile


  subroutine Print_BC_Summary ()
    implicit none

    write(*,*)
    write(*,'(A)') ' Boundary conditions'
    if (nconnect > 0) write(*,'(A,T35,I0)') '   Connection', nconnect
    if (nwall > 0) write(*,'(A,T35,I0)') '   Viscous wall', nwall
    if (nio > 0) write(*,'(A,T35,I0)') '   Inflow/outflow', nio
    if (nsym > 0) write(*,'(A,T35,I0)') '   Symmetry', nsym
    if (nper > 0) write(*,'(A,T35,I0)') '   Periodicity', nper
    if (next > 0) write(*,'(A,T35,I0)') '   Extrapolation', next
    if (ncoupled > 0) write(*,'(A,T35,I0)') '   Coupled wall', ncoupled
    if (has_tdep_bc) write(*,'(A)') '   Time-dependent BC detected'

  end subroutine Print_BC_Summary


  !─────────────────────────────────────────────────────────────────────────────
  ! parse_dir_tok: convert a direction token (read as character) to real.
  !   If the token contains 'normal', return 1.0e30 (face-normal flag).
  !   Otherwise parse as a floating-point number.
  pure function parse_dir_tok(tok) result(val)
    implicit none
    character(len=*), intent(in) :: tok
    real(R8) :: val
    integer  :: ios_loc
    character(len=len(tok)) :: tok_
    tok_ = adjustl(tok)
    if (index(trim(tok_), 'normal') > 0) then
      val = huge(1.0_R8)
    else
      read(tok_, *, iostat=ios_loc) val
      if (ios_loc /= 0) val = 0.0_R8
    endif
  end function parse_dir_tok
end module ARES_IO_BC
