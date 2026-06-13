module ARES_Lib_Multigrid
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Check_Multigrid, Coarse_Grid, Coarse_IOfield
  public :: fine2coarse_prim, coarse2fine_prim

contains

  subroutine Check_Multigrid ( domain )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_multigrid
    use ir_precision, only: str
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, d, rap, check

    rap = 2**(obj_multigrid%MGL-1)
    do b = 1, domain % nb
      check = 0
      do d = 1, 3
        if ( Mod ( domain % Blk(b) % Dim(d), rap ) == 0 ) check = check + 1
        if (d == 3) then
          if ( domain % Blk(b) % Dim(d) == 1 ) check = check + 1
        endif
      enddo
      if ( check < 3 ) then
        write(*,'(A90)') ' [ERROR] in Check_Multigrid, block: '//trim(str(.true.,b))
        stop
      endif
    enddo

  end subroutine Check_Multigrid


  subroutine fine2coarse_prim ( fPrim, cPrim, fVol, cVol, fDim, cDim )
    use ARES_Global_m
    implicit none
    integer,  intent(in),  dimension(3) :: fDim, cDim
    real(R8), intent(in),  dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(out), dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    real(R8), intent(in),  dimension(1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fVol
    real(R8), intent(in),  dimension(1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cVol
    ! Local
    integer  :: i, j, k, i2, j2, k2, i2d, j2d, k2d

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)
    do i = 1, cDim(1)

      i2 = 2*i
      j2 = 2*j
      k2 = 2*k
      i2d = i2-1
      j2d = j2-1
      k2d = k2-1

      if ( fDim(3) == 1 ) then ! 2D
        cPrim(:,i,j,k) = fPrim(:,i2d,j2d,1) * fVol(i2d,j2d,1)  &
                        + fPrim(:,i2, j2d,1) * fVol(i2, j2d,1)  &
                        + fPrim(:,i2d,j2, 1) * fVol(i2d,j2, 1)  &
                        + fPrim(:,i2, j2, 1) * fVol(i2, j2, 1)
      else ! 3D
        cPrim(:,i,j,k) = fPrim(:,i2d,j2d,k2d) * fVol(i2d,j2d,k2d)  &
                        + fPrim(:,i2, j2d,k2d) * fVol(i2, j2d,k2d)  &
                        + fPrim(:,i2d,j2, k2d) * fVol(i2d,j2, k2d)  &
                        + fPrim(:,i2, j2, k2d) * fVol(i2, j2, k2d)  &
                        + fPrim(:,i2d,j2d,k2 ) * fVol(i2d,j2d,k2 )  &
                        + fPrim(:,i2, j2d,k2 ) * fVol(i2, j2d,k2 )  &
                        + fPrim(:,i2d,j2, k2 ) * fVol(i2d,j2, k2 )  &
                        + fPrim(:,i2, j2, k2 ) * fVol(i2, j2, k2 )
      end if ! 2D

      cPrim(:,i,j,k) = cPrim(:,i,j,k) / cVol(i,j,k)

    enddo; enddo; enddo

  end subroutine fine2coarse_prim


  subroutine coarse2fine_prim ( fPrim, cPrim, fDim, cDim )
    use ARES_Global_m
    implicit none
    integer, intent(in), dimension(3) :: fDim, cDim
    real(R8), intent(out), dimension(nprim,1-gc:fDim(1)+gc,1-gc:fDim(2)+gc,1-gc:fDim(3)+gc) :: fPrim
    real(R8), intent(in),  dimension(nprim,1-gc:cDim(1)+gc,1-gc:cDim(2)+gc,1-gc:cDim(3)+gc) :: cPrim
    ! Local
    integer :: i, j, k, rap, i2, j2, k2, i2d, j2d, k2d, im, jm, km, ip, jp, kp
    integer :: ii, jj, kk, counter, mask(3), id(6)
    real(R8) :: a1, a2, a3, a4, coeffs(8), interp(nprim)

    a1 = 27d0/64d0
    a2 = 9d0/64d0
    a3 = 3d0/64d0
    a4 = 1d0/64d0
    coeffs(1:8) = [ a1, a2, a2, a2, a3, a3, a3, a4 ]
    rap = 2
    
    if ( fDim(3) == 1 ) then ! 2D
      a1 = 9.d0/16.d0
      a2 = 3.d0/16.d0
      a3 = 1.d0/16.d0
      a4 = 0d0
      coeffs(1:8) = [ a1, a2, a2, a4, a3, a4, a4, a4 ]
    end if

    !$omp do collapse (2)
    do k = 1, cDim(3)
    do j = 1, cDim(2)            
    do i = 1, cDim(1)
            
      i2 = rap*i
      j2 = rap*j
      k2 = rap*k
      i2d = i2-(rap-1)
      j2d = j2-(rap-1)
      k2d = k2-(rap-1)
          
      im = Max (1,i-1)
      jm = Max (1,j-1)
      km = Max (1,k-1)

      ip = Min (cDim(1), i+1)
      jp = Min (cDim(2), j+1)
      kp = Min (cDim(3), k+1)

      if ( fDim(3) == 1 ) then
        k2 = 1
        k2d = 1 
      end if ! 2D

      id(1:6) = [ im, jm, km, ip, jp, kp ]
          
      counter = 1
          
      do kk = k2d, k2          
      do jj = j2d, j2            
      do ii = i2d, i2

        if ( fDim(3) == 1 ) then ! 2D
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
        else
          if (counter==1) mask(1:3) = [1,2,3]
          if (counter==2) mask(1:3) = [4,2,3]
          if (counter==3) mask(1:3) = [1,5,3]
          if (counter==4) mask(1:3) = [4,5,3]
          if (counter==5) mask(1:3) = [1,2,6]
          if (counter==6) mask(1:3) = [4,2,6]
          if (counter==7) mask(1:3) = [1,5,6]
          if (counter==8) mask(1:3) = [4,5,6]
        endif  

        interp = coeffs(1) * cPrim(:,i,j,k) &
               + coeffs(2) * cPrim(:,id(mask(1)),j,k) &
               + coeffs(3) * cPrim(:,i,id(mask(2)),k)    &
               + coeffs(4) * cPrim(:,i,j,id(mask(3)))    &
               + coeffs(5) * cPrim(:,id(mask(1)),id(mask(2)),k)   &
               + coeffs(6) * cPrim(:,id(mask(1)),j,id(mask(3)))   &
               + coeffs(7) * cPrim(:,i,id(mask(2)),id(mask(3)))   &
               + coeffs(8) * cPrim(:,id(mask(1)),id(mask(2)),id(mask(3)))

        fPrim(:,ii,jj,kk) = interp

        counter = counter + 1          
        
      enddo; enddo; enddo
    
    enddo; enddo; enddo
      
  end subroutine coarse2fine_prim


  subroutine Coarse_Grid ( Fine, Coarse )
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b

    do b = 1, Coarse % nb ! Loop over blocks
      call Coarse_Grid_Blk ( Fine % Blk(b) % node, Coarse % Blk(b) % node, &
                             Fine % Blk(b) % dim,  Coarse % Blk(b) % dim )
    enddo

  contains

    subroutine Coarse_Grid_Blk ( fNode, cNode, fDim, cDim )
      use ARES_Base_Types_m
      implicit none
      integer, intent(in), dimension(3) :: fDim, cDim
      type(ARES_vector_3D_type), intent(in),  dimension(0:fDim(1),0:fDim(2),0:fDim(3)) :: fNode
      type(ARES_vector_3D_type), intent(out), dimension(0:cDim(1),0:cDim(2),0:cDim(3)) :: cNode
      ! Local
      integer :: i, j, k, i2, j2, k2

      !$omp do collapse (2)
      do k = 0, fDim(3), 2-mod(fDim(3),2)
      do j = 0, fDim(2), 2
      do i = 0, fDim(1), 2
        
        i2 = i / 2
        j2 = j / 2
        k2 = k / 2

        if ( fDim(3) == 1 ) then ! 2D
          k2 = k
        end if

        cNode(i2,j2,k2) % c = fNode(i,j,k) % c
      
      enddo; enddo; enddo
    
    end subroutine Coarse_Grid_Blk

  end subroutine Coarse_Grid


  subroutine Coarse_IOfield ( fIOfield, cIOfield )
    use Lib_ORION_data
    use ARES_Global_m
    implicit none
    type(ORION_data), intent(in)    :: fIOfield
    type(ORION_data), intent(inout) :: cIOfield
    ! Local
    integer :: b, nb, Ni, Nj, Nk, Onvar

    Onvar = nprim + 3  ! T, gamma, R

    nb = size( fIOfield % block )
    allocate ( cIOfield % block( nb ) )
    
    do b = 1, nb
      Ni = fIOfield % block(b)%Ni /2
      Nj = fIOfield % block(b)%Nj /2
      Nk = fIOfield % block(b)%Nk /2
      Nk = Max ( 1, Nk ) ! 2D case
      cIOfield % block(b) % Ni = Ni
      cIOfield % block(b) % Nj = Nj
      cIOfield % block(b) % Nk = Nk
      allocate( cIOfield%block(b)%mesh(1:3,0:Ni,0:Nj,0:Nk) )
      allocate( cIOfield%block(b)%vars(1:Onvar,Ni,Nj,Nk) )
    enddo

  end subroutine Coarse_IOfield

end module ARES_Lib_Multigrid