module ARES_Lib_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  
  !> Concrete procedures pointing to one of the subroutine realizations.
  procedure(Eddy_Viscosity_if),             pointer, public :: Eddy_Viscosity
  procedure(RANS_Diffusive_Flux_if),        pointer, public :: RANS_Diffusive_Flux
  procedure(RANS_Source_Terms_if),          pointer, public :: RANS_Source_Terms
  procedure(RANS_Set_Wall_Values_if),       pointer, public :: RANS_Set_Wall_Values
  procedure(RANS_Extrapolate_Wall_if),      pointer, public :: RANS_Extrapolate_Wall
  procedure(RANS_Enforce_Realizability_if), pointer, public :: RANS_Enforce_Realizability
  procedure(RANS_Set_Blowing_Wall_if),      pointer, public :: RANS_Set_Blowing_Wall

  !> Abstract interface relative to the eddy viscosity computation procedure.
  abstract interface
    subroutine Eddy_Viscosity_if ( mut, rans_variables, mul, rho, vel_gradient, walldist, ks )
      use iso_fortran_env, only: I4 => int32, R8 => real64
      use ARES_Global_m
      implicit none
      real(R8), intent(in), dimension(nRANS) :: rans_variables  ! : vector of RANS eqs variables in form (rho*R)
      real(R8), intent(in)                   :: mul             ! : Laminar viscosity
      real(R8), intent(in)                   :: rho             ! : Density
      real(R8), intent(in), dimension(3,3)   :: vel_gradient    ! : Velocity gradient
      real(R8), intent(in)                   :: walldist        ! : Distance nearest wall
      real(R8), intent(in)                   :: ks              ! : Roughness nearest wall
      real(R8), intent(out)                  :: mut             ! : Eddy viscosity
    end subroutine Eddy_Viscosity_if
  end interface

  !> Abstract interface relative to the computation of conservative diffusive fluxes for RANS equations.
  abstract interface
    subroutine RANS_Diffusive_Flux_if(flux, rans_variables, mul, vel_gradient, rans_gradient, rho, area, normal, dist)
      use iso_fortran_env, only: I4 => int32, R8 => real64
      use ARES_Global_m
      implicit none
      real(R8), intent(in), dimension(nrans)   :: rans_variables         ! : vector of RANS eqs variables in form (rho*R)
      real(R8), intent(in)                     :: mul                    ! : laminar viscosity
      real(R8), intent(in), dimension(3,3)     :: vel_gradient           ! : Velocity gradient
      real(R8), intent(in), dimension(nrans,3) :: rans_gradient          ! : Gradient in form dRi/dxj
      real(R8), intent(in)                     :: rho                    ! : Density
      real(R8), intent(in)                     :: area, normal(3), dist  ! : Metrics (face area, face normal, distance from wall)
      real(R8), intent(out), dimension(nrans)  :: flux                   ! : Diffusive flux
    end subroutine RANS_Diffusive_Flux_if
  end interface

  !> Abstract interface relative to the computation of RHS terms for RANS equations.
  abstract interface
    subroutine RANS_Source_Terms_if(domain)
      use ARES_Advanced_Types_m
      implicit none
      type(ARES_domain_type), intent(inout) :: domain
    end subroutine RANS_Source_Terms_if
  end interface

  !> Abstract interface relative to the computation of wall values for RANS variables.
  abstract interface
    subroutine RANS_Set_Wall_Values_if(rans_cell, mil_wall, rans_wall, metric, dist, ks)
      use ARES_Global_m
      implicit none
      real(kind=8), intent(out), dimension(nRANS)  :: rans_wall        ! : ... at wall
      real(kind=8), intent(in),  dimension(nRANS)  :: rans_cell        ! : ... at cell center
      real(kind=8), intent(in)                     :: mil_wall         ! : laminar viscosity at wall
      real(kind=8), intent(in)                     :: metric           ! MetricTensor(dir) dot_product FaceNormal
      real(kind=8), intent(in)                     :: dist             ! : cell center-wall distance
      real(kind=8), intent(in)                     :: ks               ! roughness nearest wall
    end subroutine RANS_Set_Wall_Values_if
  end interface

  abstract interface
    subroutine RANS_Extrapolate_Wall_if ( prim, wall_rans, rho, wall_rho, ghost_rans )
      use iso_fortran_env, only: I4 => int32, R8 => real64
      use ARES_Global_m
      implicit none
      real(R8), intent(in),  dimension(nRANS) :: prim, wall_rans  ! : boundary cell and wall variables
      real(R8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
      real(R8), intent(out), dimension(nRANS) :: ghost_rans ! : ghost cell RANS variables
    end subroutine RANS_Extrapolate_Wall_if
  end interface

  abstract interface
    subroutine RANS_Enforce_Realizability_if ( var )
      use iso_fortran_env, only: I4 => int32, R8 => real64
      implicit none
      real(R8), intent(inout),  dimension(:) :: var
    end subroutine RANS_Enforce_Realizability_if
  end interface

    !> Abstract interface relative to the computation of blowing wall values for RANS variables.
  abstract interface
    subroutine RANS_Set_Blowing_Wall_if(rho, mil, tau, rans_variables, mdot, dist)
      use iso_fortran_env, only: I4 => int32, R8 => real64
      use ARES_Global_m
      implicit none
      real(R8), intent(out), dimension(nRANS)  :: rans_variables
      real(R8), intent(in)                     :: rho, mil, mdot, dist
      real(R8), intent(in), dimension(3)       :: tau
    end subroutine RANS_Set_Blowing_Wall_if
  end interface

contains

  ! Setting the value of turbulent stresses at symmetry interfaces and in ghost cells
  ! NOTE: symmetry interface MUST be either an x-plane, y-plane or z-plane
  subroutine RSM_Symmetry ( Rij, Rij_w, Rij_g, n )
    use ARES_Global_m
    implicit none
    real(R8), intent(in), dimension(:) :: Rij
    real(R8), intent(out), dimension(:) :: Rij_w, Rij_g
    real(R8), intent(in), dimension(3) :: n

    Rij_w(1:3) = Rij(1:3) ! zero-gradient

    if ( n(1) > n(2)+n(3)  ) then ! x-plane
      Rij_w(4) = 0d0
      Rij_w(5) = 0d0
      Rij_w(6) = Rij(6)
      Rij_g(4) = -Rij(4)
      Rij_g(5) = -Rij(5)
      Rij_g(6) = Rij(6)
    elseif ( n(2) > n(1)+n(3) ) then ! y-plane
      Rij_w(4) = 0d0
      Rij_w(5) = Rij(5)
      Rij_w(6) = 0d0
      Rij_g(4) = -Rij(4)
      Rij_g(5) = Rij(5)
      Rij_g(6) = -Rij(6)
    elseif ( n(3) > n(1)+n(2) ) then ! y-plane
      Rij_w(4) = Rij(4)
      Rij_w(5) = 0d0
      Rij_w(6) = 0d0
      Rij_g(4) = Rij(4)
      Rij_g(5) = -Rij(5)
      Rij_g(6) = -Rij(6)
    else
      error stop ( 'Error: with RSM symmetry must be an x-, y- or z-plane.' )
    end if

  end subroutine RSM_Symmetry

end module ARES_Lib_RANS