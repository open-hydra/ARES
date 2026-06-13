module ARES_Mod_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Integration_Variables

  !> Concrete procedure pointing to one of the function realizations
  procedure(newstate_if), pointer, public :: RK_Newstate

  !> Abstract interface relative to the "Newstate" procedure
  abstract interface
    subroutine Newstate_if ( domain, irk )
      use ARES_Advanced_Types_m
      use ARES_Config_Types_m, only: obj_time_scheme, obj_irs
      use ARES_Global_m
      use FLINT_Lib_Thermodynamic
      use ARES_Lib_RK
      use ARES_Lib_RANS
      use ARES_Lib_IRS
      implicit none
      type(ARES_domain_type), intent(inout) :: domain
      integer, intent(in)                   :: irk
    end subroutine Newstate_if
  end interface

contains

  subroutine Assign_Integration_Variables ()
    use ARES_Config_Types_m, only: obj_time_scheme
    use ARES_Lib_Newstate
    implicit none

    select case (obj_time_scheme%integration_variables)
    case ('prim')
      RK_Newstate => Newstate_Primitive
      obj_time_scheme%integration_variables = 'Primitive'
    case ('prec')
      RK_Newstate => Newstate_Preconditioned
      obj_time_scheme%integration_variables = 'Preconditioned'
    end select

  end subroutine Assign_Integration_Variables

end module ARES_Mod_Newstate