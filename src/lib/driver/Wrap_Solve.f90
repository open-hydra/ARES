module ARES_Wrap_Solve
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: ARES_solve

contains

  subroutine ARES_solve ( simulation, External_Function )
    use ARES_Advanced_Types_m, only: ARES_simulation_type
    use ARES_Config_Types_m,   only: obj_time_scheme, obj_multigrid
    use ARES_Mod_Multigrid,    only: Prolongation
    use ARES_Mod_Explicit,     only: Explicit_Step
    implicit none
    type(ARES_simulation_type), intent(inout) :: simulation
    external :: External_Function

    select case (trim(obj_time_scheme%solver_type))

      case ('euler', 'RK2', 'RK3')

        ! 1. Calculate solution on local grid
        call Explicit_Step( simulation%domain, External_Function )

        ! 2. Interpolate solution on finer grid
        if (obj_multigrid%change_MG) then
          call Prolongation ( Fine=simulation%domain(obj_multigrid%MG_level-1), Coarse=simulation%domain(obj_multigrid%MG_level) )
        endif

    end select
    
  end subroutine ARES_solve

end module ARES_Wrap_Solve