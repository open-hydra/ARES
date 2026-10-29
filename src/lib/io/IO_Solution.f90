module ARES_IO_Solution

  implicit none

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(r_solution_if), pointer, public :: Read_IC
  procedure(w_solution_if), pointer, public :: Write_solution
  
  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
    subroutine r_solution_if ( IOfield )
      use Lib_ORION_data
      implicit none
      type(ORION_data), intent(inout) :: IOfield
    end subroutine r_solution_if

    subroutine w_solution_if ( domain, IOfield, file )
      use ARES_Advanced_Types_m
      use ARES_Config_Types_m, only: obj_io
      use ARES_Global_m
      use ARES_Parameters_m
      use ARES_Lib_RANS
      use IR_Precision
      use Lib_ORION_data
      use Lib_VTK
      use Lib_Tecplot
      use FLINT_Lib_Thermodynamic
      implicit none
      type(ARES_domain_type), intent(inout) :: domain
      type(ORION_data), intent(inout)    :: IOfield
      character(llen), intent(in)        :: file
    end subroutine w_solution_if
  end interface

contains

  subroutine Setup_Input_Solution()
    use ARES_Config_Types_m, only: obj_sim_param, obj_io 
    use ARES_Global_m
    use ARES_Parameters_m
    use IR_Precision
    implicit none
    logical         :: present
    integer         :: i
    character(llen) :: try
    character(6)    :: extension

    if (obj_sim_param%newrun) then
      if (index(obj_io%ini_format,'vtk')>0) then
        extension = '.vtm'
        Read_IC  => Read_vtk_tec
      else
        if (index(obj_io%ini_format,'ascii')>0) then
          extension = '.tec'
        else
          extension = '.szplt'
        endif
        Read_IC  => Read_vtk_tec
      endif
      obj_io%nameinit = 'INPUT/'//trim(ARES_phase_prefix)//'ic'//extension
      inquire(file=obj_io%nameinit, exist=present)

    else

      if (index(obj_io%sol_format,'vtk')>0) then
        extension = '.vtm'
        Read_IC => Read_vtk_tec
      else
        if (index(obj_io%sol_format,'ascii')>0) then
          extension = '.tec'
        else
          extension = '.szplt'
        endif
        Read_IC => Read_vtk_tec
      end if
      obj_io%nameinit = 'OUTPUT/'//trim(ARES_phase_prefix)//'field'//extension
      inquire(file=obj_io%nameinit, exist=present)

      if (.not.present) then
        i = 0
        do
          i = i+1
          try = 'OUTPUT/'//trim(ARES_phase_prefix)//'field'//trim(str(.true.,i))//extension
          inquire(file=try,exist=present)
          if (present) then
            obj_io%nameinit = try
          else
            exit
          endif
        enddo
      endif

    endif

  end subroutine Setup_Input_Solution

  subroutine Read_vtk_tec ( IOfield )
    use ARES_Config_Types_m, only: obj_sim_param, obj_io
    use ARES_Parameters_m
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use strings, only: parse
    implicit none
    type(ORION_data), intent(inout) :: IOfield
    ! Local
    integer         :: error
    character(clen) :: format(2)

    error = 0

    if ( obj_sim_param%newrun ) then
      call parse(obj_io%ini_format,' ', format)
    else
      call parse(obj_io%sol_format,' ', format)
    endif

    obj_io%IOtime = 0.d0

    select case(trim(format(1)))
    case('tecplot')
      IOfield%tec%format = trim(format(2))
      error = tec_read_structured_multiblock(orion=IOfield,filename=trim(obj_io%nameinit))
    case('vtk')
      IOfield%tec%format = trim(format(2))
      error = vtk_read_structured_multiblock(orion=IOfield,vtmpath=obj_io%nameinit(1:len(trim(obj_io%nameinit))-4),vtspath='INPUT/vtk/field',time=obj_io%IOtime)
    end select

    if (error/=0) obj_io%error_message = "[ERROR] reading input file "//trim(obj_io%nameinit)

  end subroutine Read_vtk_tec


  !> Output setup
  subroutine Setup_Output_Solution ( IOfield )
    use IR_Precision
    use ARES_Config_Types_m, only: obj_multigrid, obj_io, obj_prec
    use ARES_Global_m
    use ARES_Parameters_m
    use Lib_ORION_data
    use strings, only: parse
    implicit none
    type(ORION_data), intent(inout) :: IOfield(obj_multigrid%MGL)
    ! Local
    integer :: b, i, m

    ! IO Variables specification
    obj_io%Ovarnames=' '
    obj_io%ORANSname=' '
    ! The turbulence variables are the density-weighted ones ARES solves for, under the names of
    ! the density-weighted keys of ATLAS (those of the NASA Turbulence Modeling Resource), as in
    ! the initial condition written by ICB: rhoNuTilde = rho*nu_tilde (the Spalart-Allmaras
    ! variable, not the eddy viscosity, which is the auxiliary variable mut below),
    ! rhoK = rho*k, rhoOmega = rho*omega
    if ( nrans == 1 ) then
      obj_io%ORANSname=' "rhoNuTilde"'
    elseif ( nrans == 2 ) then
      obj_io%ORANSname=' "rhoK" "rhoOmega"'
    elseif ( nrans == 7 ) then
      obj_io%ORANSname=" ru'u' rv'v' rw'w' ru'v' ru'w' rv'w' rhoOmega"
    end if
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "p" "u" "v" "w" "h"'
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//obj_io%ORANSname

    ! Auxiliary variables
    ! T, rho, sound are printed by default
    obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "T" "rho" "sound"'
    obj_io%Onvar = nprim + 3 ! T, rho, sound
    ! Thermodynamic properties
    obj_io%write_thermo = .false.
    if (index(obj_io%sol_variables,'thermo')>0) then
      obj_io%write_thermo = .true.
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "cp"'
      obj_io%Onvar = obj_io%Onvar + 1
      if (obj_prec%enabled) then
        obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "Ur"'
        obj_io%Onvar = obj_io%Onvar + 1
      endif
    endif

    ! Transport properties
    obj_io%write_transport = .false.
    if (index(obj_io%sol_variables,'transport')>0) then
      obj_io%write_transport = .true.
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "mil" "kl"'
      obj_io%Onvar = obj_io%Onvar + 2
    endif
    if (obj_io%write_transport .and. (nrans>0)) then
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "mut"'   ! eddy viscosity
      obj_io%Onvar = obj_io%Onvar + 1
    endif

    do m = 1, obj_multigrid%MGL
      IOfield(m)%vtk%node = .false.
      IOfield(m)%tec%node = .false.
      IOfield(m)%tec%bc = .false.
      if (index(obj_io%sol_format,'ascii')>0) then
        IOfield(m)%tec%extension = '.tec'
        IOfield(m)%tec%format = 'ascii'
        IOfield(m)%vtk%format = 'ascii'
      else
        IOfield(m)%tec%extension = '.szplt'
        IOfield(m)%tec%format = 'binary'
        IOfield(m)%vtk%format = 'binary'
      endif
    end do

    ! Concretize the sol subroutine
    Write_solution => Write_vtk_tec

    do m = 1, obj_multigrid%MGL
      do b = 1, Size(IOfield(m)%block)
        IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
      enddo
    enddo

    ! Reallocate IOfield vars if the number of variables read into the backup file is different from the solution one
    ! The reallocation is performed after reading the ICs and before the first solution is written!
    do m = 1, obj_multigrid%MGL
      if ( obj_io%Onvar /= Size(IOfield(m)%block(1)%vars,1) ) then
        do b = 1, Size ( IOfield(m)%block )
          IOfield(m)%block(b)%name = 'Block'//trim(str(.true.,b))
          deallocate ( IOfield(m)%block(b)%vars )
          allocate( IOfield(m)%block(b)%vars(1:obj_io%Onvar, &
                                             1:IOfield(m)%block(b)%Ni, &
                                             1:IOfield(m)%block(b)%Nj, &
                                             1:IOfield(m)%block(b)%Nk) )
        enddo
      endif
    enddo

  end subroutine Setup_Output_Solution

  subroutine Write_vtk_tec ( domain, IOfield, file )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_io
    use ARES_Global_m
    use ARES_Parameters_m
    use IR_Precision
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use ARES_Mod_MPI, only: mpi_is_root, mpi_size_, is_local_block, block_owner, check_mpi_error
    use ARES_Mod_GhostExchange, only: gather_P_to_root, mpi_io_barrier
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    character(llen), intent(in)           :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    integer             :: E_IO, b

    ! Gather P from all ranks to root (collective — all ranks must call)
    call gather_P_to_root(domain)

    ! Output variables are computed by the rank that owns the block: the auxiliary
    ! ones need the metrics, the wall distance and the ghost cells, which root
    ! does not hold for remote blocks (geometry-only, interior P from the gather).
    if (mpi_is_root) then
      do b = 1, domain%nb
        if (is_local_block(b)) call Fill_Output_Vars ( domain%blk(b), domain%blk(b)%dim, IOfield%block(b)%vars )
      enddo
    endif

#ifdef USE_MPI
    if (mpi_size_ > 1) then
      block
        use mpi
        use ARES_Mod_GhostExchange, only: TAG_IO
        integer :: ierr, nreq, n, offset
        integer, allocatable :: reqs(:), stats(:,:)
        real(R8P), allocatable :: buf(:)

        allocate(reqs(domain%nb), stats(MPI_STATUS_SIZE, domain%nb))
        nreq = 0
        if (mpi_is_root) then
          ! Root: receive the remote blocks straight into IOfield
          do b = 1, domain%nb
            if (is_local_block(b)) cycle
            n = obj_io%Onvar * product(domain%blk(b)%dim)
            nreq = nreq + 1
            call MPI_IRECV(IOfield%block(b)%vars, n, MPI_DOUBLE_PRECISION, &
                           block_owner(b), TAG_IO, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
          enddo
        else
          ! Non-root: fill the local blocks in one contiguous buffer and send them
          n = 0
          do b = 1, domain%nb
            if (is_local_block(b)) n = n + obj_io%Onvar * product(domain%blk(b)%dim)
          enddo
          allocate(buf(max(n, 1)))
          offset = 0
          do b = 1, domain%nb
            if (.not. is_local_block(b)) cycle
            n = obj_io%Onvar * product(domain%blk(b)%dim)
            call Fill_Output_Vars ( domain%blk(b), domain%blk(b)%dim, buf(offset+1:offset+n) )
            nreq = nreq + 1
            call MPI_ISEND(buf(offset+1), n, MPI_DOUBLE_PRECISION, &
                           0, TAG_IO, MPI_COMM_WORLD, reqs(nreq), ierr)
            call check_mpi_error(ierr)
            offset = offset + n
          enddo
        endif
        if (nreq > 0) then
          call MPI_WAITALL(nreq, reqs(1:nreq), stats(:,1:nreq), ierr)
          call check_mpi_error(ierr)
        endif
      end block
    endif
#endif

    if (mpi_is_root) then
      path = 'OUTPUT/'

      ! Write the IOfield accordingly to the solution format
      if (index(obj_io%sol_format,'vtk')>0) then
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                              vtmpath=trim(path)//trim(file),varnames=obj_io%Ovarnames,time=domain%time)
      else
        E_IO = tec_write_structured_multiblock(orion=IOfield,varnames=obj_io%Ovarnames,filename=trim(path)//trim(file)//trim(IOfield%tec%extension))
      end if
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_vtk_tec


  !> Output variables of one block: primitives, then T, rho, the speed of sound and the optional thermo,
  !> transport and eddy-viscosity ones, in the order of obj_io%Ovarnames.
  !> Needs the full block (P with ghosts, metrics, wall distance): call it on the owner.
  subroutine Fill_Output_Vars ( blk, nijk, out )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_io, obj_prec
    use ARES_Global_m
    use ARES_Parameters_m
    use ARES_Lib_RANS
    use IR_Precision
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Preconditioning
    implicit none
    type(ARES_block_type), intent(in) :: blk
    integer, intent(in)               :: nijk(3)
    real(R8P), intent(out)            :: out(obj_io%Onvar, nijk(1), nijk(2), nijk(3))
    ! Local
    integer :: i, j, k, pv
    real(8) :: p, rho, T, mu, kl, mut, entalpy, sound
    real(8) :: vel_gradient(3,3)

    out(1:nprim,:,:,:) = blk%P(:,1:nijk(1),1:nijk(2),1:nijk(3))
    ! Auxiliary variables
    do k = 1, nijk(3) ; do j = 1, nijk(2) ; do i = 1, nijk(1)
      entalpy = blk%P(nh,i,j,k)
      p    = blk%P(np,i,j,k)
      T = ph2vars( p, entalpy, T_tab )
      rho = ph2vars( p, entalpy, rho_tab )
      sound = ph2vars( p, entalpy, sound_tab )
      pv = 3; out(nprim+1,i,j,k) = T; out(nprim+2,i,j,k) = rho
      out(nprim+3,i,j,k) = sound
      if ( Uref > 0.d0 ) out(nprim+3,i,j,k) = Uref
      if (obj_io%write_thermo) then
        pv = pv + 1; out(nprim+pv,i,j,k) = ph2vars( p, entalpy, hT_tab ) ! cp = dh/dT
        if (obj_prec%enabled) then
          pv = pv + 1; out(nprim+pv,i,j,k) = Comp_Ur(Norm2(blk%P(nu:nw,i,j,k)), sound )
        endif
      endif
      if (obj_io%write_transport) then
        mu = ph2vars( p, entalpy, mi_tab2D )
        kl = ph2vars( p, entalpy, k_tab2D )
        pv = pv +1; out(nprim+pv,i,j,k) = mu
        pv = pv +1; out(nprim+pv,i,j,k) = kl
        if (nrans>0) then
          vel_gradient(:,1) = ( blk%P(nu:nw,i+1,j,k) - blk%P(nu:nw,i-1,j,k) )/2d0
          vel_gradient(:,2) = ( blk%P(nu:nw,i,j+1,k) - blk%P(nu:nw,i,j-1,k) )/2d0
          vel_gradient(:,3) = ( blk%P(nu:nw,i,j,k+1) - blk%P(nu:nw,i,j,k-1) )/2d0
          vel_gradient = matmul ( vel_gradient + 1d-40, blk%m(i,j,k)%c )
          call Eddy_Viscosity ( mut=mut, rans_variables=blk%p(nt:,i,j,k), &
                                mul=mu, rho=rho, vel_gradient=vel_gradient, &
                                walldist=blk%yn(i,j,k), ks=blk%ks(i,j,k))
          pv = pv +1; out(nprim+pv,i,j,k) = mut
        end if
      endif
    enddo; enddo; enddo

  end subroutine Fill_Output_Vars

end module ARES_IO_Solution