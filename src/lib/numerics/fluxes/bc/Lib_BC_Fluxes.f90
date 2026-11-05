module ARES_Lib_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  public

contains

  pure subroutine Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    implicit none
    integer, intent(in) :: fm
    integer, intent(out) :: modfm, modfm1, modfm2, modfm3

    modFm  = mod(Fm,2)
    modfm1 = 1-modFm  
    modfm2 = 1-2*modFm
    modfm3 = 2*( 2*modFm-1 )

  end subroutine Compute_Modfm


  pure subroutine Face_Index ( face, dir, cell_i, cell_j, cell_k, face_i, face_j, face_k )
    implicit none
    integer, intent(in) :: face, cell_i, cell_j, cell_k
    integer, intent(out) :: dir, face_i, face_j, face_k

    select case ( face )
      case(1:2)
        dir = 1
        face_i = cell_i - mod ( face, 2 )
        face_j = cell_j
        face_k = cell_k
      case(3:4)
        dir = 2
        face_i = cell_i
        face_j = cell_j - mod ( face, 2 )
        face_k = cell_k
      case(5:6)
        dir = 3
        face_i = cell_i
        face_j = cell_j
        face_k = cell_k - mod ( face, 2 )
    end select

  end subroutine Face_Index


  !> The properties of a wall face for the wall output, the same on every face: the stress is the traction of the gas on
  !> the wall, tau . n with n the normal into the gas; qw is the heat flux into the gas (out of the solid). The callers
  !> give stress and qw along the face normal of the solver, toward growing index: stress = tau . n, qw = k grad T . n
  !> (the heat flux into the wall times into), where into = +1 when that normal points into the gas (faces 1, 3, 5) and
  !> -1 when it points into the wall (faces 2, 4, 6). A zero is written as +0 (Positive_Zero).
  subroutine Compute_Wall_Properties(stress,pw,into,Tw,rhow,mu,qw,y,exit_array)
    use iso_fortran_env, only: int64
    implicit none
    real(R8), dimension(:), intent(in)    :: stress
    real(R8), intent(in)                  :: pw
    integer, intent(in)                   :: into
    real(R8), dimension(8), intent(inout) :: exit_array
    real(R8), optional :: Tw, qw, rhow, mu, y
    ! Local
    real(R8) :: tauw, ut

    ! y+
    tauw = sqrt( sum ( Stress**2 ) )
    ut = sqrt ( tauw/rhow )
    exit_array(1) = rhow*y*ut/mu
    ! Wall stress: the traction of the gas on the wall
    exit_array(2:4) = Positive_Zero ( into * stress )
    ! Wall pressure
    exit_array(5) = pw
    ! Wall temperature
    if (present(Tw)) exit_array(6) = Tw
    ! Heat flux into the gas
    if (present(qw)) exit_array(7) = Positive_Zero ( - into * qw )

  contains

    !> x, with -0 written as +0. The wall output of a run on several ranks is gathered with a sum, which turns -0 into +0,
    !> so a -0 left here would make it differ, bit for bit, from that of a run on one rank. The test is on the bits: under a
    !> fast floating-point model (as -ffast-math allows) the compilers may drop an arithmetic one, such as x + 0.
    elemental function Positive_Zero ( x ) result ( z )
      real(R8), intent(in) :: x
      real(R8) :: z
      integer(int64) :: b
      b = transfer ( x, b )
      if ( ibclr ( b, 63 ) == 0_int64 ) b = 0_int64
      z = transfer ( b, z )
    end function Positive_Zero

  end subroutine Compute_Wall_Properties

end module ARES_Lib_BC_Fluxes