module ARES_Mod_Diagnostic
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Parameters_m

  implicit none
  private
  public :: Compute_Residual, Write_Diagnostic

  character(len=llen), private :: Dvarnames = '"dp" "du" "dv" "dw" "dh" "dt"'

contains

  !> The residual of an iteration: for each primitive variable (p, u, v, w, h and the turbulence variables), the sum over
  !> the cells of the square of its change in the iteration (new - old), added to total; average is the root mean square
  !> over the cells of the block. It is a change per iteration, not a time derivative: dt is not used, so the residual
  !> grows with cfl (and cfl-rise-threshold changes the unit of res-threshold, which is taken at the first iteration).
  !> The sum over the cells of the block is taken in a fixed order, the same for any number of threads: each thread sums
  !> whole i lines, and the sums of the lines are added in turn.
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
    real(R8) :: resn(nres), residuo(nres)
    real(R8), allocatable :: line(:,:,:)   ! sum over each i line

    allocate ( line(nres, n(2), n(3)) )

    !$omp parallel do private ( i, resn ) collapse (2)
    do k = 1, n(3)
    do j = 1, n(2)
      line(:,j,k) = 0d0
      do i = 1, n(1)

        resn(1:nres) = abs ( new(np:nprim,i,j,k) - old(np:nprim,i,j,k) ) ! delta(p - u - v - w - h - rans(if))
        line(:,j,k) = line(:,j,k) + resn*resn  ! sum ( r(i)**2 )

      enddo
    enddo; enddo
    !$omp end parallel do

    residuo = 0d0
    do k = 1, n(3)
    do j = 1, n(2)
      residuo = residuo + line(:,j,k)
    enddo; enddo

    ! Average residual of this block
    average = sqrt ( residuo / real ( n(1)*n(2)*n(3), R8 ) )
    ! Overall residual
    total = total + residuo
    deallocate ( line )

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
    character(len=llen) :: vtknames
    integer             :: E_IO, b, i, j, k, ndiag
    character(len=clen) :: format(2)
    real(R8), allocatable :: full(:,:,:,:)

    ! Gather R, dtlocal, beta from all ranks to root (collective)
    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then
      path = 'OUTPUT/'

      call parse(obj_io%sol_format,' ', format)

      ! The residual field holds nres+1 variables (5 mean-flow residuals, the
      ! nrans turbulence residuals, and dt). Build the name list to match that
      ! count; a fixed string omits the RANS residual(s) and makes the binary
      ! writer declare fewer variable names than the data written.
      ! blk%r holds, at output time, the increments of the last stage of the primitive variables ARES advances (times
      ! dt/V, smoothed by IRS when it is on): p, u, v, w, h, then the density-weighted turbulence variables
      ! (rhoNuTilde = rho*nu_tilde, rhoK, rhoOmega, as in the solution file, see IO_Solution). Each is named d<variable>
      Dvarnames = '"dp" "du" "dv" "dw" "dh"'
      if (nrans==1) then
        Dvarnames = trim(Dvarnames)//' "drhoNuTilde"'
      elseif (nrans==2) then
        Dvarnames = trim(Dvarnames)//' "drhoK" "drhoOmega"'
      elseif (nrans==7) then
        Dvarnames = trim(Dvarnames)//' "druu" "drvv" "drww" "druv" "druw" "drvw" "drhoOmega"'
      endif
      Dvarnames = trim(Dvarnames)//' "dt"'

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

      ! The residuals and dt
      ndiag = nres + 1

      ! Write the IOfield accordingly to the solution format
      E_IO = 0
      select case(trim(format(1)))
      case('vtk')
        IOfield%vtk%format = trim(format(2))
        ! The names become XML attributes (Name="..."): the double quotes of the Tecplot list must go
        vtknames = Dvarnames
        do i = 1, len_trim(vtknames)
          if ( vtknames(i:i) == '"' ) vtknames(i:i) = ' '
        enddo
        ! IOfield is the output buffer of the field, with obj_io%Onvar columns, and the vtk
        ! writer writes every column it holds: give it the ndiag diagnostic columns only,
        ! then the field's size back (the field writer fills every column at each write).
        do b = 1, size(IOfield%block)
          call move_alloc(IOfield%block(b)%vars, full)
          allocate(IOfield%block(b)%vars(1:ndiag,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk))
          IOfield%block(b)%vars(:,:,:,:) = full(1:ndiag,:,:,:)
          deallocate(full)
        enddo
        localpath_vtk = trim(path)//'/vtk'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        ! The time is the solution time of the Tecplot output: minus the global iteration count in a steady run
        E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                             vtmpath=trim(path)//trim(file),varnames=vtknames, &
                                                             time=IOfield%solutiontime)
        do b = 1, size(IOfield%block)
          deallocate(IOfield%block(b)%vars)
          allocate(IOfield%block(b)%vars(1:obj_io%Onvar,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk))
        enddo
      case('tecplot')
        IOfield%tec%format = trim(format(2))
        E_IO = tec_write_structured_multiblock(Nvars=ndiag,orion=IOfield,varnames=Dvarnames,filename=trim(path)//trim(file)//trim(IOfield%tec%extension))
      end select
      ! A failed write (disk full, quota, permission) leaves a cut file: say so, and go on (the next output may succeed)
      if ( E_IO /= 0 ) write(*,'(A,I0,A)') ' [WARNING] '//trim(path)//trim(file)//' was not written (status ', E_IO, &
                                           '): disk full, quota or permission?'
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Diagnostic

end module ARES_Mod_Diagnostic