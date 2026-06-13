module ARES_IO_Probes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Parameters_m
  use ARES_Base_Types_m, only: ARES_vector_3D_type
  use ARES_Config_Types_m, only: obj_io_probes
  use ARES_Advanced_Types_m, only: ARES_domain_type

  implicit none
  private
  public :: Setup_Probes, Write_Probes_Data

  type :: real_ptr
    real(R8), pointer :: p
  end type real_ptr

  type :: obj_probe
    type(ARES_vector_3D_type)        :: location
    integer, dimension(4)            :: ilocation
    integer                          :: nvar
    character(len=clen), allocatable :: names(:)
    type(real_ptr), allocatable      :: variables(:)
    real(R8), pointer                :: P(:)
    real(R8)                         :: T, M
    real(R8)                         :: dtime
    integer                          :: ntime
    integer                          :: diter
    integer                          :: unit
  contains
    private 
    procedure, pass(self) :: Place
    procedure, pass(self) :: Update
  end type obj_probe

  integer, public :: nprobes
  type(obj_probe), allocatable, target :: probe(:)

contains

  subroutine Setup_Probes( domain, newrun )
    use IR_precision
    use strings, only: parse
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(in)  :: domain
    logical, intent(in)                 :: newrun
    ! Local
    integer             :: i, ii, error
    character(len=clen) :: string(8)

    nprobes = size(obj_io_probes)
    allocate(probe(1:nprobes))

    do i = 1, nprobes

      obj_io_probes(i)%file = 'OUTPUT/'//trim(obj_io_probes(i)%file)//'.txt'

      ! Vars name
      string = ''
      probe(i)%nvar = 0
      call parse(obj_io_probes(i)%varnames,' ',string)
      do ii = 1, size(string)
        if (string(ii) /= '') then
          probe(i)%nvar = probe(i)%nvar + 1
        endif
      enddo
      allocate(probe(i)%names(1:probe(i)%nvar))
      do ii = 1, probe(i)%nvar
        probe(i)%names(ii) = trim(string(ii))
      enddo

      ! Frequency
      probe(i)%dtime = obj_io_probes(i)%dtime
      probe(i)%diter = obj_io_probes(i)%diter

      ! Location
      probe(i)%location%c = obj_io_probes(i)%loc
      probe(i)%ilocation  = obj_io_probes(i)%iloc

      if (sum(obj_io_probes(i)%iloc)==0) call probe(i)%Place(domain)

      ! Only the rank owning this probe's block sets up pointers and opens the file
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle

      call Assign_Variables(probe(i), domain)

      if (.not.newrun) then
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='OLD',iostat=error)
        if (error/=0) then
          obj_io_probes(i)%error_message = "[ERROR] You restarted from an old solution but the probes files were not found."
          return
        endif
        error = 0
        do while ( error == 0 ); read (probe(i)%unit, *, iostat=error); enddo
        backspace (probe(i)%unit)
      else
        open(newunit=probe(i)%unit,file=trim(obj_io_probes(i)%file),status='REPLACE',iostat=error)
      endif

    enddo

  end subroutine Setup_Probes


  !> Find the cell indexes related to the probe starting from space coords
  subroutine Place (self, domain)
    implicit none
    class(obj_probe)                   :: self
    type(ARES_domain_type), intent(in) :: domain
    ! Local
    integer :: i, j, k, b
    real(8) :: d0, d

    d0 = huge(1d0)
    do b = 1, domain%nb
      do k = 0, domain%blk(b)%dim(3); do j = 0, domain%blk(b)%dim(2); do i = 1, domain%blk(b)%dim(1)
        d = norm2(self%location%c-domain%blk(b)%node(i,j,k)%c)
        if (d<d0) then
          d0 = d
          self%ilocation = [b , i, j, k]
        endif
      enddo; enddo; enddo
    enddo

  end subroutine Place


  subroutine Assign_Variables(probe, domain)
    use IR_precision
    use ARES_Global_m, only: np, nu, nv, nw, nh
    implicit none
    type(obj_probe), intent(inout), target     :: probe
    type(ARES_domain_type), intent(in), target :: domain
    ! Local
    integer :: v,s,b,i,j,k

    b = probe%ilocation(1)
    i = probe%ilocation(2)
    j = probe%ilocation(3)
    k = probe%ilocation(4)

    probe%P => domain%blk(b)%P(:,i,j,k)
    allocate(probe%variables(1:probe%nvar))

    do v = 1, probe%nvar
      if (probe%names(v)=='p') probe%variables(v)%p => probe%P(np)
      if (probe%names(v)=='u') probe%variables(v)%p => probe%P(nu)
      if (probe%names(v)=='v') probe%variables(v)%p => probe%P(nv)
      if (probe%names(v)=='w') probe%variables(v)%p => probe%P(nw)
      if (probe%names(v)=='h') probe%variables(v)%p => probe%P(nh)
      !if (probe%names(v)=='rho') probe%variables(v)%p => probe%rho
      if (probe%names(v)=='T') probe%variables(v)%p => probe%T
      if (probe%names(v)=='M') probe%variables(v)%p => probe%M
    enddo

  end subroutine Assign_Variables


  subroutine Update(self)
    use ARES_Global_m,           only: np, nu, nv, nw, nh
    use FLINT_Lib_Thermodynamic

    implicit none
    class(obj_probe), intent(inout) :: self
    ! Local
    real(8) :: speed_of_sound

    speed_of_sound = ph2vars( self%P(np), self%P(nh), sound_tab )
    self%T = ph2vars( self%P(np), self%P(nh), T_tab )
    self%M = norm2(self%P(nu:nw))/speed_of_sound

  end subroutine Update


  subroutine Write_Probes_Data( iter, time )
    use ARES_Mod_MPI, only: is_local_block

    implicit none
    integer, intent(in) :: iter
    real(8), intent(in) :: time
    ! Local
    integer :: i, v

    do i = 1, nprobes
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
      if (mod(iter,probe(i)%diter)==0) then
        call probe(i)%update
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
      elseif (time >= probe(i)%dtime*probe(i)%ntime) then
        call probe(i)%update
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p,v=1,probe(i)%nvar)
        probe(i)%ntime = probe(i)%ntime+1
      endif
    enddo

  end subroutine Write_Probes_Data


end module ARES_IO_Probes