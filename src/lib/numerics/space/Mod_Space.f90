module ARES_Mod_Space
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Space_Scheme

contains

  subroutine Setup_Space_Scheme()
    use ARES_Config_Types_m, only: obj_space_scheme
    implicit none

    obj_space_scheme%description = 'none'
    obj_space_scheme%warning_message = 'none' 
    obj_space_scheme%error_message = 'none'

    ! MUSCL
    if (index(obj_space_scheme%space_reconstruction,'MUSCL')>0) then
      obj_space_scheme%description = 'MUSCL'
      if (obj_space_scheme%flux_limiter=='none') then
        obj_space_scheme%flux_limiter='vanleer'
        obj_space_scheme%warning_message = '[WARNING] You have specified MUSCL reconstruction but not a flux limiter. Van Leer by default.'
      endif
    ! First order
    else
      obj_space_scheme%description = 'First-order'
      if (obj_space_scheme%flux_limiter /= 'none') then
        obj_space_scheme%warning_message = '[WARNING] Flux limiter specified but MUSCL reconstruction is not enabled. Ignoring flux limiter.'
      end if
    end if

    call Assign_Limiter(obj_space_scheme%flux_limiter)

  end subroutine Setup_Space_Scheme


  subroutine Assign_Limiter(marker)
    use ARES_Lib_Limiters
    implicit none
    character(len=*), intent(inout) :: marker

    select case (marker)
      case default
        rlimiter => rlimiter_IORD
      case ('minmod')
        rlimiter => rlimiter_MINMOD
        marker = 'Min-Mod'
      case ('vanleer')
        rlimiter => rlimiter_VANLEER
        marker = 'Van Leer'
      case ('vanalbada')
        rlimiter => rlimiter_VANALBADA
        marker = 'Van Albada'
      case ('MC')
        rlimiter => rlimiter_MC
        marker = 'MC'
      case ('LIMO3')
        rlimiter => rlimiter_LIMO3
        marker = 'LIMO3'
      case ('superbee')
        rlimiter => rlimiter_SB
        marker = 'Super-Bee'
    end select

  end subroutine Assign_Limiter

end module ARES_Mod_Space