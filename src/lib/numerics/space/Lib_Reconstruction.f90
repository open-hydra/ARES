module ARES_Lib_Reconstruction
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Reconstruction

contains

  elemental subroutine Reconstruction ( prev, local, next, next2, dl0, dl1, dl2, dll, dlr, priml, primr )
    use ARES_Lib_Limiters, only : rlimiter
    implicit none
    real(R8), intent(in)    :: prev, local, next, next2
    real(R8), intent(inout) :: priml, primr
    real(R8), intent(in) :: dl0, dl1, dl2, dll, dlr
    ! Local
    real(R8) :: slope0, slope1, slope2, slopel, sloper

    ! Piecewise Linear Reconstruction
    slope0  = ( local - prev ) / dl0
    slope1  = ( next - local ) / dl1
    slope2  = ( next2 - next ) / dl2
    slopel = rlimiter ( slope1, slope0 )
    sloper = rlimiter ( slope2, slope1 )
    priml = local + slopel * dll
    primr =  next - sloper * dlr

  end subroutine Reconstruction

end module ARES_Lib_Reconstruction