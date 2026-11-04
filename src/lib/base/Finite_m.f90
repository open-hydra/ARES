!> The test of a real value for NaN and infinities by its bits.
module ARES_Finite_m
  use iso_fortran_env, only: R8 => real64

  implicit none
  private
  public :: Not_Finite

contains

  !> .true. for NaN and infinities, the values with all the exponent bits set. The test reads the bits:
  !> -ffast-math (-ffinite-math-only) lets gfortran fold isnan, ieee_is_nan and ieee_is_finite to constants. Reading
  !> the bits is an operation on integers, which raises no floating-point exception, and it is no comparison, which a
  !> NaN passes as false: the state update (Lib_Newstate) and the readers of the input values (the initial condition in
  !> Wrap_Setup, the time tables in Series_Data_m, the lines of bc.txt in IO_BC, the reals of input.ini in the registry)
  !> check their values with it, also in a build that traps the invalid operations, such as the DEBUG build.
  elemental logical function Not_Finite ( x )
    use iso_fortran_env, only: int64
    implicit none
    real(R8), intent(in) :: x
    integer(int64), parameter :: expo = int(z'7FF', int64)

    Not_Finite = iand ( shiftr ( transfer ( x, 0_int64 ), 52 ), expo ) == expo

  end function Not_Finite

end module ARES_Finite_m
