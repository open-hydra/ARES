module ARES_Lib_Limiters
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

  !> Concrete "rlimiter" procedure pointing to one of the function realizations
  procedure(rlimiter_if), pointer :: rlimiter

  !> Abstract interface relative to the "rlimiter" procedure
  abstract interface
    pure function rlimiter_if(x,y) result(result)
      use iso_fortran_env, only: I4 => int32, R8 => real64
      implicit none
      real(R8), intent(in) :: x, y
      real(R8) :: ru
      real(R8) :: phi
      real(R8) :: result
    end function rlimiter_if
  end interface

contains

  !> First order scheme
  pure function rlimiter_IORD(x,y) result(result)
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: result
    
    result = 0.0d0

  end function rlimiter_IORD


  !> Min-Mod flux limiter
  pure function rlimiter_MINMOD(x,y) result(result)
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: result
    
    result = 0.d0
    if (x*y<0.d0) return
    result = x
    if (abs(x)>abs(y)) result = y

  end function rlimiter_MINMOD


  !> Van Leer flux limiter
  pure function rlimiter_VANLEER(x,y) result(result)
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x + 1d-30)
    phi = (ru+abs(ru))/(1d0+abs(ru))
    result = phi*x

  end function rlimiter_VANLEER


  !> Van Albada flux limiter
  pure function rlimiter_VANALBADA(x,y) result(result)
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x+1d-30)
    phi = (ru*ru+ru)/(1d0+ru*ru)
    result = phi*x

  end function rlimiter_VANALBADA


  !> MC flux limiter
  pure function rlimiter_MC(x,y) result(result)
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,minval([2d0*ru,0.5d0*(1d0+ru),2d0])])
    result = phi*x

  end function rlimiter_MC


  !> Super-bee flux limiter
  pure function rlimiter_SB(x,y) result(result)
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval([0d0,min(2d0*ru,1d0),min(ru,2d0)])
    result = phi*x

  end function rlimiter_SB

  
  !> 3rd order limiter
  pure function rlimiter_LIMO3(x,y) result(result)
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    ru = y/(x+1d-30)
    phi = maxval((/0d0,minval((/((2.d0+ru)/3.d0),maxval((/-0.5d0*ru,0.d0,minval((/2.d0*ru,((2.d0+ru)/3.d0),1.6d0/))/))/))/))
    result = phi*x

  end function rlimiter_LIMO3

end module ARES_Lib_Limiters