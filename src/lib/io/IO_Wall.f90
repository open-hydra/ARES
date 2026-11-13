module ARES_IO_Wall
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Config_Types_m, only: obj_io, obj_io_bc
  use ARES_Global_m
  use ARES_Parameters_m
  use Lib_ORION_data
  implicit none
  ! Wall output, one container per multigrid level (index = domain%mg_level): allocated on
  ! every rank, filled and written on root only
  type(ORION_data), allocatable :: IOwall(:)

contains

  subroutine Initialize_Wall_File ( domain, extension )
    use ARES_Advanced_Types_m
    use IR_Precision
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    character(len=*), intent(in)       :: extension
    ! Local
    integer :: b, f, bb, i, j, k, lv

    lv = domain%mg_level
    IOwall(lv)%tec%extension = extension

    obj_io%write_wall_mechanical = .false.
    obj_io%write_wall_thermal = .false.
    obj_io%write_wall_mass = .false.

    obj_io%Onwall = 1
    obj_io%Owallnames = '"y+"'
    if (index(obj_io%wall_variables,'mech')>0) then
      obj_io%Onwall = obj_io%Onwall + 4
      obj_io%Owallnames = trim(obj_io%Owallnames)//' "tauX" "tauY" "tauZ" "pw"'
      obj_io%write_wall_mechanical = .true.
    endif
    if (index(obj_io%wall_variables,'thermal')>0) then
      obj_io%Onwall = obj_io%Onwall + 2
      obj_io%Owallnames = trim(obj_io%Owallnames)//' "Tw" "qw"'
      obj_io%write_wall_thermal = .true.
    endif

    ! Sand-grain roughness of the wall cell (always written; 0 if the case is
    ! not rough). Kept as the last wall variable.
    obj_io%Onwall = obj_io%Onwall + 1
    obj_io%Owallnames = trim(obj_io%Owallnames)//' "hs"'

    allocate(IOwall(lv)%block(1:count(obj_io_bc%viscous_flag)))

    bb = 0
    do b = 1, domain % nb; do f = 1, 6
      if (obj_io_bc%viscous_flag( b , f )) then
        bb = bb + 1
        IOwall(lv)%block(bb)%name = 'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))
        select case(f)
        case(1)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall(lv) % block(bb) % mesh(1,0,j,k) = domain % blk(b) % node(0,j,k) % c(1)
              IOwall(lv) % block(bb) % mesh(2,0,j,k) = domain % blk(b) % node(0,j,k) % c(2)
              IOwall(lv) % block(bb) % mesh(3,0,j,k) = domain % blk(b) % node(0,j,k) % c(3)
              ! Section to fool compiler...
              if (IOwall(lv)%block(bb)%mesh(1,0,j,k)==0d0 .and. &
                  IOwall(lv)%block(bb)%mesh(2,0,j,k)==0d0 .and. &
                  IOwall(lv)%block(bb)%mesh(3,0,j,k)==0d0) &
              print*, IOwall(lv) % block(bb) % mesh(:,0,j,k)
          enddo; enddo
        case(2)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:0,0:domain%blk(b)%dim(2),0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1,1:domain%blk(b)%dim(2),1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2)
              IOwall(lv) % block(bb) % mesh(:,0,j,k) = domain % blk(b) % node(domain%blk(b)%dim(1),j,k) % c
          enddo; enddo
        case(3)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,0,k) % c
          enddo; enddo
        case(4)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:0,0:domain%blk(b)%dim(3)))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1,1:domain%blk(b)%dim(3)))
          do k = 0, domain%blk(b)%dim(3); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,0,k) = domain % blk(b) % node(i,domain%blk(b)%dim(2),k) % c
          enddo; enddo
        case(5)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,0) % c
          enddo; enddo
        case(6)
          allocate(IOwall(lv)%block(bb)%mesh(1:3,0:domain%blk(b)%dim(1),0:domain%blk(b)%dim(2),0:0))
          allocate(IOwall(lv)%block(bb)%vars(1:obj_io%Onwall,1:domain%blk(b)%dim(1),1:domain%blk(b)%dim(2),1))
          do j = 0, domain%blk(b)%dim(2); do i = 0, domain%blk(b)%dim(1)
              IOwall(lv) % block(bb) % mesh(:,i,j,0) = domain % blk(b) % node(i,j,domain%blk(b)%dim(3)) % c
          enddo; enddo
        end select
        ! Cells of the face in each direction, 0 across it: the vtk writer of ORION sizes the block
        ! file with them (the Tecplot writer takes the sizes from the mesh array)
        IOwall(lv)%block(bb)%Ni = size(IOwall(lv)%block(bb)%mesh, 2) - 1
        IOwall(lv)%block(bb)%Nj = size(IOwall(lv)%block(bb)%mesh, 3) - 1
        IOwall(lv)%block(bb)%Nk = size(IOwall(lv)%block(bb)%mesh, 4) - 1
      endif
    enddo; enddo

  end subroutine Initialize_Wall_File


  !> Update the orion-field data and write them accordingly to the chosen format (vtk,tecplot)
  subroutine Write_Wall_Solution( domain, file )
    use Lib_VTK
    use Lib_Tecplot
    use ARES_Advanced_Types_m
    use ARES_Parameters_m
    use ARES_Mod_MPI, only: mpi_is_root
    use ARES_Mod_GhostExchange, only: mpi_io_barrier
    use strings, only: parse
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    character(len=llen), intent(in)       :: file
    ! Local
    character(len=llen) :: path
    character(len=llen) :: localpath_vtk
    character(len=llen) :: name
    character(len=len(obj_io%Owallnames)) :: vtknames
    integer             :: E_IO, i, lv
    character(len=clen) :: format(2)

    lv = domain%mg_level

    ! Each rank computes wall properties for its local blocks and reduces to root.
    ! (collective — all ranks must call)
    call Extract_Properties_From_BC( domain )

    if (mpi_is_root) then
      path = 'OUTPUT/'

      call parse(obj_io%sol_format,' ', format)

      ! Write the IOwall accordingly to the solution format
      select case(trim(format(1)))
      case('vtk')
        IOwall(lv)%vtk%format = trim(format(2))
        ! The names become XML attributes (Name="..."): the double quotes of the Tecplot list must go
        vtknames = obj_io%Owallnames
        do i = 1, len_trim(vtknames)
          if ( vtknames(i:i) == '"' ) vtknames(i:i) = ' '
        enddo
        ! As for the field (IO_Solution): the .vtm lists its blocks relative to its folder
        name = file
        if ( name(1:1) == '/' ) name = name(2:)
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        ! The time is the solution time of the Tecplot output: minus the global iteration count in a steady run
        E_IO = vtk_write_structured_multiblock(orion=IOwall(lv),vtspath=trim(localpath_vtk)//trim(name), &
                                               vtmpath=trim(path)//trim(name),varnames=vtknames,time=IOwall(lv)%solutiontime)
      case('tecplot')
        IOwall(lv)%tec%format = trim(format(2))
        E_IO = tec_write_structured_multiblock(orion=IOwall(lv),varnames=obj_io%Owallnames, &
                                               filename=trim(path)//trim(file)//trim(IOwall(lv)%tec%extension))
      end select
    end if

    ! Synchronize all ranks after I/O
    call mpi_io_barrier()

  end subroutine Write_Wall_Solution


  !> Compute wall BC properties on all ranks (local blocks only), reduce to root,
  !> and fill IOwall on root. Collective — all ranks must call. The residual and
  !> the ghost cells are left as they were: an output does not change the state.
  subroutine Extract_Properties_From_BC ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Mod_MPI, only: mpi_is_root, is_local_block, mpi_reduce_sum_r8_array
    use ARES_Lib_BC_Fluxes_Wall_Heat
    use ARES_Lib_BC_Fluxes_Wall_Temperature
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    real(R8), dimension(8) :: Ovar
    integer :: f, lower, upper, i, b, bb, pv, Bm, Im, Jm, Km, Fm, Is, Js, Ks
    integer :: n_wall_total, wall_i, ii, nrec, lv
    real(R8), allocatable :: wall_buf(:)
    integer, allocatable  :: order(:), first(:,:), next(:,:)
    real(R8) :: r_keep(nprim)   ! residual of the wall cell, before the call

    lv = domain%mg_level

    ! Count total BC cells across all viscous faces
    n_wall_total = 0
    do b = 1, domain%nb
      do f = 1, 6
        if (obj_io_bc%viscous_flag(b,f)) n_wall_total = n_wall_total + domain%n_bf(b,f)
      end do
    end do

    allocate(wall_buf(8 * n_wall_total))
    wall_buf = 0d0

    ! The records grouped by block, then face, in the order of the file within a group
    ! (a counting sort over n_bf): the wall properties are computed and laid out face by
    ! face, block by block, whatever the order of bc.txt. The ranges n_bf(b,f) of the
    ! record list itself hold one face of one block only when bc.txt is sorted.
    allocate(order(domain%nbound), first(domain%nb,6), next(domain%nb,6))
    nrec = 0
    do b = 1, domain%nb
      do f = 1, 6
        first(b,f) = nrec + 1
        nrec = nrec + domain%n_bf(b,f)
      end do
    end do
    next = first
    do i = 1, domain%nbound
      b = domain%bc(i)%b
      f = domain%bc(i)%f
      order(next(b,f)) = i
      next(b,f) = next(b,f) + 1
    end do

    ! Each rank computes Ovar only for its local blocks
    wall_i = 0
    do b = 1, domain%nb
      do f = 1, 6
        lower = first(b,f)
        upper = lower + domain%n_bf(b,f) - 1
        if (.not. obj_io_bc%viscous_flag(b,f)) cycle

        do ii = lower, upper
          i = order(ii)
          wall_i = wall_i + 1
          Bm = domain%bc(i)%b
          if (.not. is_local_block(Bm)) cycle

          Im = domain%bc(i)%i
          Jm = domain%bc(i)%j
          Km = domain%bc(i)%k
          Fm = domain%bc(i)%f

          ! The wall flux routines below are called for the wall properties (Ovar)
          ! only, but they also add the wall flux to the residual of the wall cell.
          ! The diagnostic output, written next, prints that residual: keep it and
          ! put it back. (They write no ghost value: Fill_Ghost_Cell does.)
          r_keep = domain % blk(Bm) % r (:,Im,Jm,Km)

          Ovar = 0d0

          ! Roughness of the wall cell (0 for a smooth case; blk%ks is allocated
          ! whenever model>0, i.e. always for viscous walls).
          if (allocated(domain%blk(Bm)%ks)) Ovar(8) = domain%blk(Bm)%ks(Im,Jm,Km)

          select case ( domain%bc(i)%type )

            case (301) ! wall: prescribed heat flux
              call BC_Wall_Heat ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % qw, Ovar, rec = domain % bc(i) )

            case (302) ! wall: prescribed temperature
              call BC_Wall_Temperature ( Im, Jm, Km, Fm, domain % blk(Bm), domain % bc(i) % Tw, Ovar, rec = domain % bc(i) )

          end select

          domain % blk(Bm) % r (:,Im,Jm,Km) = r_keep

          wall_buf(8*(wall_i-1)+1 : 8*wall_i) = Ovar
        enddo

      enddo ! face
    enddo ! block

    ! Reduce wall_buf to root. Each cell is owned by exactly one rank, so SUM gives
    ! correct values. On root, wall_buf is updated in-place by mpi_reduce_sum_r8_array.
    call mpi_reduce_sum_r8_array(wall_buf, 8*n_wall_total)

    ! Root fills IOwall from the gathered buffer
    if (mpi_is_root) then
      wall_i = 0; bb = 0
      do b = 1, domain%nb
        do f = 1, 6
          lower = first(b,f)
          upper = lower + domain%n_bf(b,f) - 1
          if (.not. obj_io_bc%viscous_flag(b,f)) cycle
          bb = bb + 1

          do ii = lower, upper
            i = order(ii)
            wall_i = wall_i + 1
            Ovar = wall_buf(8*(wall_i-1)+1 : 8*wall_i)

            select case(f)
              case(1,2); Is = 1                ; Js = domain%bc(i)%j ; Ks = domain%bc(i)%k
              case(3,4); Is = domain%bc(i)%i   ; Js = 1              ; Ks = domain%bc(i)%k
              case(5,6); Is = domain%bc(i)%i   ; Js = domain%bc(i)%j ; Ks = 1
            end select
            IOwall(lv)%block(bb)%vars(1, Is, Js, Ks) = Ovar(1)
            pv = 1
            if (obj_io%write_wall_mechanical) then
              IOwall(lv)%block(bb)%vars(pv+1:pv+4, Is, Js, Ks) = Ovar(2:5)
              pv = pv + 4
            endif
            if (obj_io%write_wall_thermal) then
              IOwall(lv)%block(bb)%vars(pv+1:pv+2, Is, Js, Ks) = Ovar(6:7)
              pv = pv + 2
            endif
            ! Roughness — always last
            IOwall(lv)%block(bb)%vars(pv+1, Is, Js, Ks) = Ovar(8)
            pv = pv + 1
          enddo

        enddo ! face
      enddo ! block
    end if

    deallocate(wall_buf, order, first, next)

  end subroutine Extract_Properties_From_BC
  
end module ARES_IO_Wall