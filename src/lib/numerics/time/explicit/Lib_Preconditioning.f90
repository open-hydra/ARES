module ARES_Lib_Preconditioning
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Parameters_m

  implicit none
  
  private
  public :: comp_Ur

  real(R8), parameter, private :: eps_min = 0.10d0 ! Default minimum eps

contains


  function comp_Ur ( vel, sound ) result ( Ur )
    use ARES_Global_m
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_RANS

    implicit none
    real(R8), intent(in) :: vel, sound
    real(R8) :: Ur
    ! Local
    real(R8) :: eps, Mach, Uref_, epsm
    
    ! Reference velocity
    Uref_ = Uref
    if ( Uref_ < 0.0d0 ) Uref_ = sound
    ! Reference local mach
    epsm = emin
    if ( epsm < 0.0d0 ) epsm = eps_min

    Mach = vel/Uref_
    eps = min(1.d0, max( epsm, Mach**2 ) )
    if ( vel < eps*Uref_ ) then
      Ur = Uref_*eps 
    elseif (vel > Uref_) then
      Ur = Uref_
    else
      Ur = vel
    end if
    
    end function comp_Ur


end module ARES_Lib_Preconditioning
