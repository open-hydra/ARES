module ARES_Mod_Diagnostic
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Parameters_m

  implicit none
  private
  public :: Compute_Residual, Write_Diagnostic

  character(len=llen), private :: Dvarnames = '"rho" "rhou" "rhov" "rhow" "rhoe" "dt"'

contains

  subroutine Compute_Residual ( new, old, dt, n, average, total )
    use ARES_Global_m
    implicit none
    integer, intent(in)     :: n(3)
    real(R8), intent(in)    :: new ( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(R8), intent(in)    :: old ( nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc )
    real(R8), intent(in)    :: dt ( n(1), n(2), n(3) )
    real(R8), intent(out)   :: average(nres)
    real(R8), intent(inout) :: total(nres)
    ! Local
    integer :: i, j, k
    real(R8) :: resn(nres), residuo(nres), reslocal ( nres, n(1), n(2), n(3) )

    residuo = 0d0
    
    !$omp parallel
    !$omp do private ( resn ) reduction ( + : residuo ) collapse (3)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      
      resn(1:nres) = abs ( new(np:nprim,i,j,k) - old(np:nprim,i,j,k) ) ! delta(p - u - v - w - h - rans(if))
      reslocal(:,i,j,k) = resn
      residuo = residuo + resn*resn  ! sum ( r(i)**2 )

    enddo; enddo; enddo
    !$omp end parallel

    ! Local max residual (UNUSED)
    ! do i = 1, nres; resmax_(i) = maxval ( reslocal(i,:,:,:) ); enddo

    ! Average residual of this block
    average = sqrt ( residuo / float ( n(1)*n(2)*n(3) ) )
    ! Overall residual
    total = total + residuo

  end subroutine Compute_Residual


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine Write_Diagnostic( domain, IOfield, file )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_io
    use ARES_Global_m
    use Lib_ORION_data
    use Lib_VTK
    use Lib_Tecplot
    use ARES_Mod_MPI, only: mpi_is_root
    use ARES_Mod_GhostExchange, only: gather_diagnostic_to_root, mpi_io_barrier
    use strings, only: parse
    implicit none
    type(ARES_domain_type), intent(inout)  :: domain
    type(ORION_data), intent(inout)        :: IOfield
    character(llen), intent(in)            :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    integer             :: E_IO, b, i, j, k
    character(len=clen) :: format(2)

    ! Gather R, dtlocal, beta from all ranks to root (collective)
    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      call parse(obj_io%sol_format,' ', format)

      ! Update IOfield variables with domain residuals
      do b = 1, size(IOfield%block)
        IOfield%block(b)%vars(1:5,:,:,:) = domain%blk(b)%r(np:nh,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        if ( nrans > 0 ) then ! turbulence variables
          IOfield%block(b)%vars(6:6+nrans-1,:,:,:) = domain%blk(b)%r(nt:nprim,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk)
        end if
        ! Auxiliary variables: dt
        do k = 1, IOfield%block(b)%Nk ; do j = 1, IOfield%block(b)%Nj ; do i = 1, IOfield%block(b)%Ni
              IOfield%block(b)%vars(nres+1,i,j,k) = domain%blk(b)%dtlocal(i,j,k)
        enddo; enddo; enddo
      enddo

      ! Write the IOfield accordingly to the solution format
      select case(trim(format(1)))
      case('vtk')
        IOfield%vtk%format = trim(format(2))
        localpath_vtk = trim(path)//'/vtk'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=Dvarnames,time=domain%time)
      case('tecplot')
        IOfield%tec%format = trim(format(2))
        E_IO = tec_write_structured_multiblock(Nvars=nres+1,orion=IOfield,varnames=Dvarnames,filename=trim(path)//trim(file)//'.tec')
      end select
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Diagnostic

end module ARES_Mod_Diagnostic