module ARES_Mod_Space
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Space_Scheme

contains

  subroutine Setup_Space_Scheme()
    use ARES_Config_Types_m, only: obj_space_scheme
    implicit none
    ! The limiter in use: flux-limiter keeps the value of input.ini, which is read again during the run and checked against its
    ! rule
    character(len=len(obj_space_scheme%flux_limiter)) :: limiter

    obj_space_scheme%description = 'none'
    obj_space_scheme%warning_message = 'none' 
    obj_space_scheme%error_message = 'none'

    limiter = obj_space_scheme%flux_limiter
    ! MUSCL
    if (index(obj_space_scheme%space_reconstruction,'MUSCL')>0) then
      obj_space_scheme%description = 'MUSCL'
      ! flux-limiter omitted: the registry default is none
      if (len_trim(limiter)==0 .or. limiter=='none') then
        limiter='vanleer'
        obj_space_scheme%warning_message = '[WARNING] You have specified MUSCL reconstruction but not a flux limiter. Van Leer by default.'
      endif
    ! First order
    else
      obj_space_scheme%description = 'First-order'
      if (len_trim(limiter)>0 .and. limiter /= 'none') then
        obj_space_scheme%warning_message = '[WARNING] Flux limiter specified but MUSCL reconstruction is not enabled. Ignoring flux limiter.'
      end if
      ! First order means zero slopes: no limiter, so that Assign_Limiter selects rlimiter_IORD
      limiter = ''
    end if

    call Assign_Limiter(limiter, obj_space_scheme%limiter_name)

    ! Chimera faces (102), and with mean those between blocks, take the inverse-mean viscous metric; with stencil the rest take
    ! that of ARES_Lib_Stencil_Metrics, but 300 faces across one cell take the cell metric
    obj_space_scheme%stencil_metric = ( trim(obj_space_scheme%diffusive_metric) == 'stencil' )
    if ( obj_space_scheme%stencil_metric ) then
      obj_space_scheme%metric_mode = 2
    else
      obj_space_scheme%metric_mode = 0
    end if

  end subroutine Setup_Space_Scheme


  !> Points rlimiter at the limiter that key names (zero slopes for any other key) and gives its name for the set-up summary
  subroutine Assign_Limiter(key, name)
    use ARES_Lib_Limiters
    implicit none
    character(len=*), intent(in)  :: key
    character(len=*), intent(out) :: name

    name = ''
    select case (key)
      case default
        rlimiter => rlimiter_IORD
      case ('minmod')
        rlimiter => rlimiter_MINMOD
        name = 'Min-Mod'
      case ('vanleer')
        rlimiter => rlimiter_VANLEER
        name = 'Van Leer'
      case ('vanalbada')
        rlimiter => rlimiter_VANALBADA
        name = 'Van Albada'
      case ('MC', 'mc')
        rlimiter => rlimiter_MC
        name = 'MC'
      case ('LIMO3')
        rlimiter => rlimiter_LIMO3
        name = 'LIMO3'
      case ('superbee')
        rlimiter => rlimiter_SB
        name = 'Super-Bee'
    end select

  end subroutine Assign_Limiter

end module ARES_Mod_Space