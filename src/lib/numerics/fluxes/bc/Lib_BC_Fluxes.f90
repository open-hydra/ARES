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


  subroutine Compute_Wall_Properties(stress,pw,Tw,rhow,mu,qw,y,exit_array)
    implicit none
    real(R8), dimension(:), intent(in)    :: stress
    real(R8), intent(in)                  :: pw
    real(R8), dimension(8), intent(inout) :: exit_array
    real(R8), optional :: Tw, qw, rhow, mu, y
    ! Local
    real(R8) :: tauw, ut

    ! y+
    tauw = sqrt( sum ( Stress**2 ) )
    ut = sqrt ( tauw/rhow )
    exit_array(1) = rhow*y*ut/mu
    ! Wall stress
    exit_array(2:4) = stress
    ! Wall pressure
    exit_array(5) = pw
    ! Wall temperature
    if (present(Tw)) exit_array(6) = Tw
    ! Heat flux
    if (present(qw)) exit_array(7) = qw

  end subroutine Compute_Wall_Properties

end module ARES_Lib_BC_Fluxes