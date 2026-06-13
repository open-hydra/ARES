!> @brief Module for Spalart-Allmaras Prt correction for rough walls, see https://doi.org/10.1016/j.ast.2022.107672
module ARES_Lib_Prt_Correction
    use iso_fortran_env, only: I4 => int32, R8 => real64

    implicit none

    private
    public :: delta_Prt

    real(R8), parameter, private :: &
        a1 = -0.0002346d0, &
        a2 =  0.002102d0,  &
        a3 =  0.003542d0,  &
        b1 = -0.002303d0,  & 
        b2 =  0.05588d0,   & 
        b3 = -0.003043d0,  & 
        vk =  0.41d0,      &  ! von Kármán constant
        nk = dexp(1.3325d0)   ! Nikuradse constant 
    contains 

    pure function delta_Prt ( nil, cp, cond, nit, hs, yn ) result( dPrt )
        implicit none
 
        real(R8), intent(in) :: nil     ! dinamic viscosity
        real(R8), intent(in) :: cp      ! specific heat @ constant pressure
        real(R8), intent(in) :: cond    ! thermal conductivity
        real(R8), intent(in) :: nit     ! SA variable
        real(R8), intent(in) :: hs      ! roughness
        real(R8), intent(in) :: yn      ! Wall distance
        real(R8)             :: dPrt    ! delta Prt considering rughness
        ! Local
        real(R8)             :: Pr      ! Prandtl number
        real(R8)             :: dUplus  ! delta Dimensionless velocity considering roughness
        real(R8)             :: a, b    ! Correlation coefficient
        
        ! Local Prandtl number and related coefficients
        Pr  = nil * cp / cond
        a   = a1*Pr*Pr + a2*Pr + a3
        b   = b1*Pr*Pr + b2*Pr + b3 

        ! Roughness correction
        dUplus  = 1.0d0/vk * dlog( 1.0d0 + ( nit * hs ) / ( vk * nil * nk * ( yn + 0.03d0 * hs + 1.0d-20 ) ) )
        dPrt    = ( a*dUplus*dUplus + b*dUplus ) * dexp ( -yn / ( hs + 1.0d-20 ) ) 

    end function delta_Prt

end module ARES_Lib_Prt_Correction
