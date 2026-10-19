module ARES_Backend_INI
  implicit none
  private
  public :: Load_Ini, Scan_Ini, Obsolete_Keys

contains

  subroutine Load_Ini(fini)
    use ARES_Input_Registry, only: reg
    use Finer, only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error

    do i = 1, reg%size
      error = 1

      ! Each registry entry has exactly one associated typed pointer.
      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
      else if (associated(reg%params(i)%value%iarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%iarr, error=error)
      else if (associated(reg%params(i)%value%rarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%rarr, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.
    end do

  end subroutine


  subroutine Scan_Ini(fini, nprobes, probes_name, nmgl)
    use Finer, only: file_ini
    !use ARES_Lib_RotatingFrame, only: obj_rot
    use ARES_Config_Types_m, only: obj_sim_param
    use ARES_Global_m
    use IR_precision
    implicit none
    type(file_ini), intent(in) :: fini
    integer, intent(out) :: nprobes, nmgl
    character(len=16), allocatable, intent(out) :: probes_name(:)
    ! Local
    character(len=llen) :: wholestring
    integer :: error, i

    ! Read the MGL (HYDRA coupling reads it from a dedicated section)
    nmgl = 1
    if (.not. obj_sim_param%HYDRA_MG) then
      call fini%get(section_name=trim(codename)//'-Multigrid', option_name='levels', val=nmgl, error=error)
    else
      call fini%get(section_name='HYDRA-Multigrid', option_name='levels', val=nmgl, error=error)
    endif

    ! Count the probes
    nprobes = 0
    do 
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,nprobes+1)), val=wholestring, error=error)
      if (error/=0) exit
      nprobes = nprobes+1
    enddo

    allocate(probes_name(nprobes))
    do i = 1, nprobes
      call fini%get(section_name=trim(codename)//'-Probes', option_name='probe'//trim(str(.true.,i)), val=probes_name(i), error=error)
    enddo
    
    ! Count the stationary faces for the rotating frame
    !obj_rot%n_stationary = 0
    !i = 0
    !do
    !  call fini%get(section_name=trim(codename)//'-Rotating-Frame',option_name='stationary-face-'//trim(str(.true.,i+1)), val=wholestring, error=error)
    !  if (error /= 0) exit
    !  obj_rot%n_stationary = obj_rot%n_stationary + 1
    !enddo
    !i = i + 1 
    !if (obj_rot%n_stationary > 0) allocate(obj_rot%stationary_face_str(obj_rot%n_stationary))
    
  end subroutine Scan_Ini


  !> Keys of the input layout ARES had before it took the one of MOSE. The registry no longer
  !  reads them, so they would be ignored without a word (a turbulent case would run laminar):
  !  Check_Input stops on the message returned here, empty when there is none.
  function Obsolete_Keys(fini) result(out)
    use Finer, only: file_ini
    use ARES_Global_m
    implicit none
    type(file_ini), intent(in) :: fini
    character(len=hlen) :: out
    ! Local
    character(len=:), allocatable :: s, item

    out = ''
    s = trim(codename)
    item = new_line('a')//'      - '
    if (fini%index(section_name=s//'-Parameters', option_name='simulation-type') > 0) &
      out = trim(out)//item//'simulation-type of ['//s//'-Parameters] is now equations = euler or navier-stokes in ['// &
            s//'-Physics] (laminar: navier-stokes, turbulent: navier-stokes with turbulence = <model>)'
    if (fini%has_section(section_name=s//'-RANS')) &
      out = trim(out)//item//'['//s//'-RANS] is gone: turbulence-model is now turbulence in ['//s//'-Physics], Sc is schmidt in ['// &
            s//'-Physics], the other keys are in ['//s//'-Turbulence]'
    if (fini%index(section_name=s//'-IO', option_name='ini-format') > 0) &
      out = trim(out)//item//'ini-format of ['//s//'-IO] is now ic-format'
    if (len_trim(out) > 0) &
      out = '[ERROR] input.ini uses keys of the earlier ARES layout, ARES now reads the same keys as MOSE:'//trim(out)

  end function Obsolete_Keys

end module ARES_Backend_INI