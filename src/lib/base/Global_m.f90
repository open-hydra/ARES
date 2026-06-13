module ARES_Global_m
  use ARES_Parameters_m

  implicit none

  character(clen) :: ARES_phase_prefix = ''

  integer :: np=1        ! Number of species - Pressure index
  integer :: nu, nv, nw  ! Velocity/momentum components indexes in prim/residuals vector
  integer :: nh          ! Entalpy/energy index in prim/residuals vector
  integer :: nt, nrans   ! First RANS variable index in prim/residuals vector and number of RANS variables
  integer :: nprim       ! Number of primitive variables
  integer :: nres        ! Number of variables whose residuals are printed
  integer :: ndir        ! Number of dimensions of computational frame
  integer :: gc=2        ! ghost cells
  real(8) :: Uref, emin  ! Preconditioning parameters
  integer :: model       ! Model index (e.g. 0: Euler, 1: laminar NS, 2:RANS, etc)
  ! EOS variables - thermo and transport properties
  !real(8) :: pmin, pmax, deltap, hmin, hmax, deltah
  !real(8) :: pmin2,pmax2,deltap2,Tminn, Tmaxx, deltaT
  !integer :: intp, inth, intp2, intT

end module ARES_Global_m