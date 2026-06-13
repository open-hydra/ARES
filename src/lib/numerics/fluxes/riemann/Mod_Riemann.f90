module ARES_Mod_Riemann
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Riemann_Solver

  !> Concrete riemann solver procedure pointing to one of the subroutine realizations
  procedure(riemann_if), pointer, public :: Riemann

  !> Abstract interface relative to the riemann solver procedure
  abstract interface
    subroutine riemann_if(pl,ul,vl,wl,hl,pr,ur,vr,wr,hr,nx,ny,nz,F_r,F_u,F_v,F_w,F_E )
      use iso_fortran_env, only: I4 => int32, R8 => real64
      use ARES_Global_m
      use FLINT_Lib_Thermodynamic
      implicit none
      integer s
      real(kind=8), intent(in)  :: hl,ul,vl,wl,pl 
      real(kind=8), intent(in)  :: hr,ur,vr,wr,pr 
      real(kind=8), intent(out) :: F_r, F_u, F_v, F_w, F_e
      real(kind=8), intent(in)  :: nx, ny, nz

    end subroutine riemann_if
  end interface

contains

  subroutine Assign_Riemann_Solver()
    use ARES_Config_Types_m, only: obj_riemann
    use ARES_Lib_Riemann_HLL
    use ARES_Lib_Riemann_LF
    implicit none

    nullify(Riemann)

    select case (obj_riemann%description)

    !! Lax-Friedrichs-type solvers
    case ('Rusanov')
      Riemann => riemann_LLF
      obj_riemann%description = 'Local Lax-Friedrichs (Rusanov)'

    case ('PLLF')
      Riemann => riemann_PLLF
      obj_riemann%description = 'Preconditioned Local Lax-Friedrichs'

    !! HLL-type solvers
    case ('HLLE')
      Riemann => riemann_HLLE
      obj_riemann%description = 'HLLE'

    case ('HLLC Prec')
      Riemann => riemann_HLLCprec
      obj_riemann%description = 'Preconditioned HLLC'

    case ('HLLC Rotated')
      Riemann => riemann_HLLCHLLE
      obj_riemann%description = 'Rotated HLLC Batten / HLLE'

    case ('HLLC')
      Riemann => riemann_HLLC
      obj_riemann%description = 'HLLC Batten'

    !! Default: unknown solver
    case default
      obj_riemann%error_message = '[ERROR] Unknown Riemann solver: '//trim(obj_riemann%description)//'. Valid options: Rusanov, PLLF, HLLE, HLLC, HLLC Prec, HLLC Rotated.'
      Riemann => riemann_HLLC
      obj_riemann%description = 'HLLC Batten (fallback)'


    end select

  end subroutine Assign_Riemann_Solver

end module ARES_Mod_Riemann