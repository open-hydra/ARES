!> @brief Module for Quadratic Constitutive Relations.
module ARES_Lib_QCR2000
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Lib_Fluid, only: Vorticity_Vector, Vorticity_Tensor

  implicit none
  real(R8), parameter, private :: C_cr1 = 0.3

contains

  pure function Stress_Vector_QCR2000 ( Gradvel, Normal, mi_l, mi_t, rhok ) result ( Stress )
    implicit none
    real(kind=8), intent(in) :: Gradvel(3,3), Normal(3), mi_l, mi_t, rhok(:)
    real(kind=8) :: Stress(3)
    ! Local
    real(kind=8) :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz, lam, Divervel, lamDivervel
    real(kind=8) :: Tau_xx_QCR, Tau_xy_QCR, Tau_xz_QCR, Tau_yy_QCR, Tau_yz_QCR, Tau_zz_QCR
    real(kind=8) :: omega(3), W_ij(3,3), O_ij(3,3), Den

    Divervel = Gradvel(1,1) + Gradvel(2,2) + Gradvel(3,3)

    ! Turbulent stress tensor
    Tau_xx = 2d0 * mi_t * ( Gradvel(1,1) - 1d0/3d0 * Divervel )
    Tau_yy = 2d0 * mi_t * ( Gradvel(2,2) - 1d0/3d0 * Divervel )
    Tau_zz = 2d0 * mi_t * ( Gradvel(3,3) - 1d0/3d0 * Divervel )
    Tau_xy = mi_t * ( Gradvel(1,2) + Gradvel(2,1) )
    Tau_xz = mi_t * ( Gradvel(1,3) + Gradvel(3,1) )
    Tau_yz = mi_t * ( Gradvel(2,3) + Gradvel(3,2) )

    ! QCR correction
    omega = Vorticity_Vector ( Gradvel )
    W_ij  = Vorticity_Tensor ( omega )
    Den = Sqrt ( Sum ( Gradvel**2 ) ) + 1d-40 ! Just to avoid division by zero
    O_ij = 2d0 * W_ij / Den
    
    ! T(i,j)_QCR = T(i,j) - C_cr1 * [ O(i,k) * T(j,k) + O(j,k) * T(i,k) ] with sum on k
    Tau_xx_QCR = Tau_xx - 2d0 * ( O_ij(1,2) * Tau_xy + O_ij(1,3) * Tau_xz ) * C_cr1
    Tau_yy_QCR = Tau_yy - 2d0 * ( O_ij(2,1) * Tau_xy + O_ij(2,3) * Tau_yz ) * C_cr1
    Tau_zz_QCR = Tau_zz - 2d0 * ( O_ij(3,1) * Tau_xz + O_ij(3,2) * Tau_yz ) * C_cr1

    Tau_xy_QCR = Tau_xy - ( O_ij(1,2) * Tau_yy + O_ij(1,3) * Tau_yz + &
                            O_ij(2,1) * Tau_xx + O_ij(2,3) * Tau_xz ) * C_cr1
    Tau_xz_QCR = Tau_xz - ( O_ij(1,2) * Tau_yz + O_ij(1,3) * Tau_zz + &
                            O_ij(3,1) * Tau_xx + O_ij(3,2) * Tau_xy ) * C_cr1
    Tau_yz_QCR = Tau_yz - ( O_ij(2,1) * Tau_xz + O_ij(2,3) * Tau_zz + &
                            O_ij(3,1) * Tau_xy + O_ij(3,2) * Tau_yy ) * C_cr1

    ! Laminar stress tensor
    lam = -2d0/3d0 * mi_l
    lamDivervel = lam*Divervel
 
    ! Laminar + turbulent stress tensor
    Tau_xx = lamDivervel + 2d0 * mi_l * Gradvel(1,1) + Tau_xx_QCR
    Tau_yy = lamDivervel + 2d0 * mi_l * Gradvel(2,2) + Tau_yy_QCR
    Tau_zz = lamDivervel + 2d0 * mi_l * Gradvel(3,3) + Tau_zz_QCR
    Tau_xy = mi_l * ( Gradvel(1,2) + Gradvel(2,1) ) + Tau_xy_QCR
    Tau_xz = mi_l * ( Gradvel(1,3) + Gradvel(3,1) ) + Tau_xz_QCR
    Tau_yz = mi_l * ( Gradvel(2,3) + Gradvel(3,2) ) + Tau_yz_QCR            
 
    Stress(1) = Dot_Product ( [Tau_xx, Tau_xy, Tau_xz], Normal )
    Stress(2) = Dot_Product ( [Tau_xy, Tau_yy, Tau_yz], Normal )
    Stress(3) = Dot_Product ( [Tau_xz, Tau_yz, Tau_zz], Normal )

  end function Stress_Vector_QCR2000
  
end module ARES_Lib_QCR2000
