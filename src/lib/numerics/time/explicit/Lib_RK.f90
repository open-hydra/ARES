module ARES_Lib_RK
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: RK_stage

  real(R8), parameter, private :: RKcoeff(3,3) = reshape([  1._R8,   0._R8,      0._R8,  &
                                                            1._R8, 0.50_R8,      0._R8,  &
                                                            1._R8, 0.25_R8,2._R8/3._R8 ],&
                                                            shape(RKcoeff), order=[2,1])

contains

  pure function RK_stage ( irk, n_rk, state, oldstate, residual ) result ( newstate )
    use ARES_Global_m, only: nprim
    implicit none
    integer, intent(in) :: irk
    integer, intent(in) :: n_rk
    real(R8), dimension(nprim), intent(in) :: state, oldstate, residual
    ! Local
    real(R8), dimension(nprim) :: newstate, residual_

    residual_ = ( state - oldstate + residual ) * rkcoeff(n_rk, irk)
    newstate = oldstate + residual_
    
  end function RK_stage
      
end module ARES_Lib_RK
