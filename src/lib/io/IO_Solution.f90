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

    obj_io%ini_format=obj_io%sol_format

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
    if ( nrans == 1 ) then
      obj_io%ORANSname=' "mi_tilde"'
    elseif ( nrans == 2 ) then
      obj_io%ORANSname=' "kappa" "omega"'
    elseif ( nrans == 7 ) then
      obj_io%ORANSname=" ru'u' rv'v' rw'w' ru'v' ru'w' rv'w' omega"
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
      obj_io%Ovarnames = trim(obj_io%Ovarnames)//' "mit"'
      obj_io%Onvar = obj_io%Onvar + 1
    endif
    ! Composition
    ! write_composition = .false.
    ! if (index(sol_variables,'compos')>0) then
    !   write_composition = .true.
    !   Ovarnames = trim(Ovarnames)//' "mu" "kl"'
    !   Onvar = Onvar + 2
    ! endif

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
    use ARES_Config_Types_m, only: obj_io, obj_prec
    use ARES_Global_m
    use ARES_Parameters_m
    use ARES_Lib_RANS
    use IR_Precision
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Preconditioning
    use ARES_Mod_MPI, only: mpi_is_root
    use ARES_Mod_GhostExchange, only: gather_P_to_root, mpi_io_barrier
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    character(llen), intent(in)           :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    integer             :: E_IO, b, i, j, k, pv
    real(8)             :: p, rho, T, mu, kl, mut, entalpy, sound
    real(8)             :: vel_gradient(3,3)

    ! Gather P from all ranks to root (collective — all ranks must call)
    call gather_P_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      ! Update IOfield variables with domain primitives and other variables
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1:nprim,:,:,:) = domain%blk(b)%P(:,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        ! Auxiliary variables
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
          entalpy = domain%blk(b)%P(nh,i,j,k)
          p    = domain%blk(b)%P(np,i,j,k)
          T = ph2vars( p, entalpy, T_tab )
          rho = ph2vars( p, entalpy, rho_tab )
          sound = ph2vars( p, entalpy, sound_tab )
          pv = 3; IOfield%block(b)%vars(nprim+1,i,j,k) = T; IOfield%block(b)%vars(nprim+2,i,j,k) = rho
          IOfield%block(b)%vars(nprim+3,i,j,k) = sound
          if ( Uref > 0.d0 ) IOfield%block(b)%vars(nprim+3,i,j,k) = Uref
          if (obj_io%write_thermo) then
            pv = pv + 1; IOfield%block(b)%vars(nprim+pv,i,j,k) = ph2vars( p, entalpy, hT_tab ) ! cp = dh/dT
            if (obj_prec%enabled) then
              pv = pv + 1; IOfield%block(b)%vars(nprim+pv,i,j,k) = Comp_Ur(Norm2(domain%blk(b)%P(nu:nw,i,j,k)), sound )
            endif
          endif
          if (obj_io%write_transport) then
            !call co_k_mi_lam_Wilke(rhoi,rho,T,mu,kl)
            mu = ph2vars( p, entalpy, mi_tab2D )
            kl = ph2vars( p, entalpy, k_tab2D )
            pv = pv +1; IOfield%block(b)%vars(nprim+pv,i,j,k) = mu
            pv = pv +1; IOfield%block(b)%vars(nprim+pv,i,j,k) = kl
            if (nrans>0) then
              vel_gradient(:,1) = ( domain%blk(b)%P(nu:nw,i+1,j,k) - domain%blk(b)%P(nu:nw,i-1,j,k) )/2d0
              vel_gradient(:,2) = ( domain%blk(b)%P(nu:nw,i,j+1,k) - domain%blk(b)%P(nu:nw,i,j-1,k) )/2d0
              vel_gradient(:,3) = ( domain%blk(b)%P(nu:nw,i,j,k+1) - domain%blk(b)%P(nu:nw,i,j,k-1) )/2d0
              vel_gradient = matmul ( vel_gradient + 1d-40, domain%blk(b)%m(i,j,k)%c )
              call Eddy_Viscosity ( mut=mut, rans_variables=domain%blk(b)%p(nt:,i,j,k), &
                                    mul=mu, rho=rho, vel_gradient=vel_gradient, &
                                    walldist=domain%blk(b)%yn(i,j,k), ks=domain%blk(b)%ks(i,j,k))
              pv = pv +1; IOfield%block(b)%vars(nprim+pv,i,j,k) = mut
            end if
          endif
        enddo; enddo; enddo
      enddo


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

end module ARES_IO_Solution