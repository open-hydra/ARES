module ARES_Mod_Allocate_Data
  use, intrinsic :: iso_fortran_env, only : iostat_end

  implicit none
  private
  public :: Setup_Data_Structure, Allocate_Block, deallocate_remote_computation_data

contains

  subroutine Setup_Data_Structure ( domain, IOfield )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use Lib_ORION_data
    use ARES_Config_Types_m, only: obj_io

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    type(ORION_data), intent(inout)       :: IOfield
    ! Local
    integer :: b, d, nblocks


    nu = np + 1
    nv = nu + 1
    nw = nv + 1
    nh = nw + 1
    nprim = nh + nrans
    if (nrans > 0) then
      nt = nh + 1
    else
      nt = nprim
    endif
    nres = nprim ! Primitive

    ! Domain is the ARES-alias of IOfield
    nblocks = size ( IOfield%block )
    allocate( domain%blk( 1:nblocks ) )
    domain%nb = nblocks

    ! Define the ijk dimensions of each block
    do b = 1, nblocks
      domain%blk(b)%dim(1) = IOfield%block(b)%Ni
      domain%blk(b)%dim(2) = IOfield%block(b)%Nj
      domain%blk(b)%dim(3) = IOfield%block(b)%Nk
    enddo

    ! Check if number of simulation variables in orion-field matches ARES expectation
    if ( size(IOfield%block(1)%vars, 1) < nprim ) then
      write(*,'(A)')         '[ERROR] Number of variables in IOfield does not match ARES expectation.'
      write(*,'(A,I0,A,I0)') '        Expected: ', nprim, ', Found: ', size(IOfield%block(1)%vars, 1)
      stop
    end if

    do b = 1, nblocks
      write(*,'(A9,I2,A3,3I5)') ' - Block ', b, ' = ', domain%blk(b)%dim

      ! Allocate domain block
      call Allocate_Block ( domain%blk(b), domain%blk(b)%dim )
      
      ! Import domain-block nodes from orion-field
      do d = 1, 3
        domain%blk(b)%node(0:IOfield%block(b)%Ni,0:IOfield%block(b)%Nj,0:IOfield%block(b)%Nk)%c(d) &
        = IOfield%block(b)%mesh(d,0:,0:,0:)
      enddo
      ! Import domain-block primitives from orion-field variables
      do d = 1, nprim
        domain%blk(b)%P(d,1:IOfield%block(b)%Ni,1:IOfield%block(b)%Nj,1:IOfield%block(b)%Nk) &
        = IOfield%block(b)%vars(d,1:,1:,1:)
      enddo
      domain%time = obj_io%IOtime

    enddo

  end subroutine Setup_Data_Structure


  subroutine Allocate_Block( blk, nijk )
    use ARES_Global_m
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_rans, obj_irs
    
    implicit none
    integer, intent(in)                  :: nijk(3)
    type(ARES_block_type), intent(inout) :: blk
    ! Local
    integer :: ni, nj, nk


    ni = nijk(1) ; nj = nijk(2) ; nk = nijk(3)

    ! Metrics
    allocate( blk % node ( 0:ni, 0:nj, 0:nk ) )
    allocate( blk % M    ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dl   ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % vol  ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % dir(1) % f (0:ni, nj, nk) )
    allocate( blk % dir(2) % f (ni, 0:nj, nk) )
    allocate( blk % dir(3) % f (ni, nj, 0:nk) )
    if (model>0) then 
      allocate ( blk % yn ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) ) ! Distance from the nearest wall
      allocate ( blk % ks ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) ) ! roughness of the nearest wall
      blk % yn = 0.0_R8
      blk % ks = 0.0_R8
    end if
    
    ! Prim and Residuals
    allocate( blk % P (nprim, 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
    allocate( blk % PO, blk % R, mold = blk % P )

    ! Dt cell center with no ghost cells
    allocate( blk % dtlocal ( 1:ni, 1:nj, 1:nk) )

    ! Velocity gradient and related, no ghost cells
    if ( obj_rans%SpalartShur ) then
      allocate( blk % vel_gradient ( 1-gc:ni+gc, 1-gc:nj+gc, 1-gc:nk+gc ) )
      allocate( blk % rc_term1(ni, nj, nk) )
      allocate( blk % rc_term2(ni, nj, nk) )
    end if

    ! Temp storage for residuals in IRS
    if ( obj_irs%enabled ) then
      allocate( blk % RS1, blk % RS2, mold = blk % R )
    end if

  end subroutine Allocate_Block


  subroutine deallocate_remote_computation_data(domain)
    use ARES_Advanced_Types_m
    use ARES_Mod_MPI, only: is_local_block, mpi_is_root

    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    integer :: b, d, i, c
    logical, allocatable :: needs_remote_P(:)

    ! Build mask of remote blocks whose P (and dir) must be kept:
    !  - chimera (102): donorID(:,1) can reference remote blocks
    !  - manifold (501): bc%bs can be remote
    allocate(needs_remote_P(domain%nb))
    needs_remote_P = .false.

    do b = 1, domain%nb
      if (is_local_block(b)) cycle

      ! Computation arrays — free on all ranks
      if (allocated(domain%blk(b)%PO))           deallocate(domain%blk(b)%PO)
      if (allocated(domain%blk(b)%R))            deallocate(domain%blk(b)%R)
      if (allocated(domain%blk(b)%RS1))          deallocate(domain%blk(b)%RS1)
      if (allocated(domain%blk(b)%RS2))          deallocate(domain%blk(b)%RS2)
      if (allocated(domain%blk(b)%dtlocal))      deallocate(domain%blk(b)%dtlocal)
      if (allocated(domain%blk(b)%vel_gradient)) deallocate(domain%blk(b)%vel_gradient)
      if (allocated(domain%blk(b)%rc_term1))     deallocate(domain%blk(b)%rc_term1)
      if (allocated(domain%blk(b)%rc_term2))     deallocate(domain%blk(b)%rc_term2)

      ! P — free unless this remote block is a chimera donor or manifold source
      if (.not. needs_remote_P(b)) then
        if (allocated(domain%blk(b)%P)) deallocate(domain%blk(b)%P)
      end if

      ! Metrics — free on non-root ranks only (root needs them for wall I/O)
      ! Keep dir on blocks needed for manifold (BC_Manifold reads blk(Bs)%dir%f%A)
      if (.not. mpi_is_root) then
        if (allocated(domain%blk(b)%node)) deallocate(domain%blk(b)%node)
        if (allocated(domain%blk(b)%M))    deallocate(domain%blk(b)%M)
        if (allocated(domain%blk(b)%dl))   deallocate(domain%blk(b)%dl)
        if (allocated(domain%blk(b)%vol))  deallocate(domain%blk(b)%vol)
        if (allocated(domain%blk(b)%yn))   deallocate(domain%blk(b)%yn)
        if (allocated(domain%blk(b)%ks))   deallocate(domain%blk(b)%ks)
        if (.not. needs_remote_P(b)) then
          do d = 1, 3
            if (allocated(domain%blk(b)%dir(d)%f)) deallocate(domain%blk(b)%dir(d)%f)
          end do
        end if
      end if
    end do

    deallocate(needs_remote_P)

  end subroutine deallocate_remote_computation_data

end module ARES_Mod_Allocate_Data