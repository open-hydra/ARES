!>@brief: Fluid mechanics relations.
module ARES_Lib_Fluid
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  real(kind=8), parameter :: negativetwothird = -2.d0/3.d0 

  ! Concrete procedure pointing to one of the function realizations
  procedure (Stress_Vector_if), pointer, public :: Stress_Vector

  ! Abstract interface relative to the "Stress_Vector" procedure
  abstract interface
    pure function Stress_Vector_if ( Gradvel, Normal, mi_l, mi_t, var ) result ( Stress )
      implicit none
      real(kind=8), intent(in) :: Gradvel(3,3), Normal(3), mi_l, mi_t, var(:)
      real(kind=8) :: Stress(3)
    end function Stress_Vector_if
  end interface

contains

  pure function Stress_Vector_Std ( Gradvel, Normal, mi_l, mi_t, var ) result ( Stress )
    implicit none
    real(kind=8), intent(in) :: Gradvel(3,3), Normal(3), mi_t, mi_l, var(:)
    real(kind=8) :: Stress(3)
    ! Local
    real(kind=8) :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
    real(kind=8) :: mi, lam, Divervel, lamDivervel

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1)
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2)
    Tau_zz = lamDivervel + 2d0 * mi * Gradvel(3,3)
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )
    Tau_xz = mi * ( Gradvel(1,3) + Gradvel(3,1) )
    Tau_yz = mi * ( Gradvel(2,3) + Gradvel(3,2) )

    Stress(1) = Dot_Product ( [Tau_xx, Tau_xy, Tau_xz], Normal )
    Stress(2) = Dot_Product ( [Tau_xy, Tau_yy, Tau_yz], Normal )
    Stress(3) = Dot_Product ( [Tau_xz, Tau_yz, Tau_zz], Normal )
  end function Stress_Vector_Std


  pure function Stress_Vector_2eq ( Gradvel, Normal, mi_l, mi_t, rhok ) result ( Stress )
    implicit none
    real(kind=8), intent(in) :: Gradvel(3,3), Normal(3), mi_t, mi_l, rhok(:)
    real(kind=8) :: Stress(3)
    ! Local
    real(kind=8) :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
    real(kind=8) :: mi, lam, Divervel, lamDivervel, rhokdelta

    mi = mi_l + mi_t
    lam = negativetwothird * mi
    Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)
    lamDivervel = lam*Divervel
    rhokdelta = negativetwothird * rhok(1)

    Tau_xx = lamDivervel + 2d0 * mi * Gradvel(1,1) + rhokdelta
    Tau_yy = lamDivervel + 2d0 * mi * Gradvel(2,2) + rhokdelta
    Tau_zz = lamDivervel + 2d0 * mi * Gradvel(3,3) + rhokdelta
    Tau_xy = mi * ( Gradvel(1,2) + Gradvel(2,1) )
    Tau_xz = mi * ( Gradvel(1,3) + Gradvel(3,1) )
    Tau_yz = mi * ( Gradvel(2,3) + Gradvel(3,2) )

    Stress(1) = Dot_Product ( [Tau_xx, Tau_xy, Tau_xz], Normal )
    Stress(2) = Dot_Product ( [Tau_xy, Tau_yy, Tau_yz], Normal )
    Stress(3) = Dot_Product ( [Tau_xz, Tau_yz, Tau_zz], Normal )
  end function Stress_Vector_2eq


  pure function Stress_Vector_RSM ( Gradvel, Normal, mi_l, mi_t, tij ) result ( Stress )
    implicit none
    real(kind=8), intent(in) :: Gradvel(3,3), Normal(3), mi_t, mi_l, tij(:)
    real(kind=8) :: Stress(3)
    ! Local
    real(kind=8) :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
    real(kind=8) :: lam, Divervel, lamDivervel

    lam = negativetwothird * mi_l
    Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)
    lamDivervel = lam*Divervel

    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) - tij(1)
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) - tij(2)
    Tau_zz = lamDivervel + 2d0 * mi_l * Gradvel(3,3) - tij(3)
    Tau_xy = mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) - tij(4)
    Tau_xz = mi_l * ( Gradvel(1,3) + Gradvel(3,1) ) - tij(5)
    Tau_yz = mi_l * ( Gradvel(2,3) + Gradvel(3,2) ) - tij(6)

    Stress(1) = Dot_Product ( [Tau_xx, Tau_xy, Tau_xz], Normal )
    Stress(2) = Dot_Product ( [Tau_xy, Tau_yy, Tau_yz], Normal )
    Stress(3) = Dot_Product ( [Tau_xz, Tau_yz, Tau_zz], Normal )
  end function Stress_Vector_RSM


  pure function Vorticity_Vector ( Gradvel ) result ( omega )
    implicit none
    real(kind=8), intent(in), dimension(3,3) :: Gradvel
    real(kind=8), dimension(3) :: omega

    omega(1) = Gradvel(3,2) - Gradvel(2,3)
    omega(2) = Gradvel(1,3) - Gradvel(3,1)
    omega(3) = Gradvel(2,1) - Gradvel(1,2)
  end function


  pure function Vorticity_Tensor ( omega ) result ( W )
    implicit none
    real(kind=8), intent(in), dimension(3) :: omega
    real(kind=8), dimension(3,3) :: W

    W(1,1) = 0.d0 ; W(2,2) = 0.d0 ; W(3,3) = 0.d0

    W(1,2) = -0.5d0*omega(3)
    W(2,1) =  0.5d0*omega(3)

    W(1,3) =  0.5d0*omega(2)
    W(3,1) = -0.5d0*omega(2)

    W(2,3) =  0.5d0*omega(1)
    W(3,2) = -0.5d0*omega(1)
  end function


  pure function Strain_Tensor ( Gradvel ) result ( E )
    implicit none
    real(kind=8), intent(in), dimension(3,3) :: Gradvel
    real(kind=8), dimension(3,3) :: E

    E(1,1) = Gradvel(1,1)
    E(2,2) = Gradvel(2,2)
    E(3,3) = Gradvel(3,3)

    E(1,2) = 0.5d0*(Gradvel(1,2) + Gradvel(2,1))
    E(1,3) = 0.5d0*(Gradvel(1,3) + Gradvel(3,1))
    E(2,3) = 0.5d0*(Gradvel(2,3) + Gradvel(3,2))

    E(2,1) = E(1,2)
    E(3,1) = E(1,3)
    E(3,2) = E(2,3)

  end function
  
end module ARES_Lib_Fluid