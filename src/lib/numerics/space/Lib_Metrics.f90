module ARES_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Base_Types_m
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m

  implicit none
  real(R8), public :: delthe   ! grid axisymmetric angle

contains

  subroutine Check_Mesh_Type( domain )
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    real(R8) :: theta1, theta2, theta(2)

    ! Mesh definition (2D,2Daxi,3D)
    if (domain%Blk(1)%dim(3)>1) then
      ! 3D
      ndir = 3
      associate( node => domain%Blk(1)%node, jm => domain%Blk(1)%dim(2), km => domain%Blk(1)%dim(3) )
        theta1 = atan2( node(0,1,0)%c(3),node(0,1,0)%c(2) ) ! computes axisymmetric angle
        theta2 = atan2( node(0,1,km)%c(3),node(0,1,km)%c(2) )
        theta(1) = theta2-theta1
        theta1 = atan2( node(0,jm,0)%c(3), node(0,jm,0)%c(2) )
        theta2 = atan2( node(0,jm,km)%c(3), node(0,jm,km)%c(2) )
        theta(2) = theta2-theta1
      end associate
      if ( (theta(1)-theta(2)) < 1.d-5 ) then
        delthe = theta(1)
      else
        delthe = 0d0
      endif
    elseif (domain%Blk(1)%dim(3)==1 .and. domain%Blk(1)%dim(2)==1) then
      ! 1D
      ndir = 1
    else
      ! 2D
      ndir = 2
      associate( node => domain%Blk(1)%node, jm => domain%Blk(1)%dim(2) )
        theta1 = atan2( node(0,1,0)%c(3),node(0,1,0)%c(2) ) ! computes axisymmetric angle
        theta2 = atan2( node(0,1,1)%c(3),node(0,1,1)%c(2) )
        theta(1) = theta2-theta1
        theta1 = atan2( node(0,jm,0)%c(3), node(0,jm,0)%c(2) )
        theta2 = atan2( node(0,jm,1)%c(3), node(0,jm,1)%c(2) )
      end associate
      theta(2) = theta2-theta1
      if ( (theta(1)-theta(2)) < 1.d-5 ) then
        ! 2Daxi
        delthe = theta(1)
      else
        ! 2D
        delthe = 0.d0
      endif
    endif

  end subroutine Check_Mesh_Type


  subroutine Compute_Norm_Area( b )
    implicit none
    class(ARES_block_type) :: b  
    ! Local
    real(R8) :: Ai, snix, sniy, sniz, Aj, snjx, snjy, snjz, Ak, snkx, snky, snkz
    integer      :: i, j, k, im, jm, km
    real(R8) :: d1(3), d2(3), d3(3)
    real(R8) :: snixx, sniyy, snizz, snjxx, snjyy, snjzz, snkxx, snkyy, snkzz
    real(R8) :: scal, signi, signj, signk
  
    ! Peliminary operations
    im = b%dim(1)
    jm = b%dim(2)
    km = b%dim(3)
  
    ! compute sign of normal vectors to the intefaces
  
    ! i direction
    d1 = b%node(1,1,1)%c - b%node(1,0,0)%c
    d2 = b%node(1,0,1)%c - b%node(1,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snixx=d3(1)
    sniyy=d3(2)
    snizz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(1,0,1)%c + b%node(1,0,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(0,1,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snixx+d3(2)*sniyy+d3(3)*snizz
  
    signi=sign(1.d0,scal)

    ! j direction
    d1 = b%node(0,1,1)%c - b%node(1,1,0)%c
    d2 = b%node(1,1,1)%c - b%node(0,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snjxx=d3(1)
    snjyy=d3(2)
    snjzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,1,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(1,0,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snjxx+d3(2)*snjyy+d3(3)*snjzz
  
    signj=sign(1.d0,scal)

    ! k direction
    d1 = b%node(0,1,1)%c - b%node(1,0,1)%c
    d2 = b%node(0,0,1)%c - b%node(1,1,1)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snkxx=d3(1)
    snkyy=d3(2)
    snkzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(1,0,1)%c )
    d2 = 0.25d0*( b%node(1,1,0)%c + b%node(0,1,0)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snkxx+d3(2)*snkyy+d3(3)*snkzz
  
    signk=sign(1.d0,scal)
    
    ! compute metrics: n, A
    
    !$omp parallel private (d1,d2,d3,i,j,k,snix,sniy,sniz,Ai,Aj,snjx,snjy,snjz,Ak,snkx,snky,snkz)
    
    ! i direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 1, jm ; do i = 0, im
  
      d1 = b%node(i,j,k)%c - b%node(i,j-1,k-1)%c
      d2 = b%node(i,j-1,k)%c - b%node(i,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ai = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Ai == 0d0) then
        snix = 0d0
        sniy = 0d0
        sniz = 0d0
      else
        snix = d3(1)/Ai*signi
        sniy = d3(2)/Ai*signi
        sniz = d3(3)/Ai*signi
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(1)%f(i,j,k)%a = Ai
      b%dir(1)%f(i,j,k)%n = [ snix, sniy, sniz ]
      
    enddo; enddo; enddo
  
    ! j direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 0, jm ; do  i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j,k-1)%c
      d2 = b%node(i,j,k)%c - b%node(i-1,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Aj = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Aj == 0d0) then
        snjx = 0d0
        snjy = 0d0
        snjz = 0d0
      else
        snjx = d3(1)/Aj*signj
        snjy = d3(2)/Aj*signj
        snjz = d3(3)/Aj*signj
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(2)%f(i,j,k)%A = Aj
      b%dir(2)%f(i,j,k)%n = [ snjx, snjy, snjz ]
  
    enddo; enddo; enddo
  
    ! k direction
    !$omp do collapse(3)
    do k = 0, km ; do j = 1, jm ; do i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j-1,k)%c
      d2 = b%node(i-1,j-1,k)%c - b%node(i,j,k)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ak = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Ak == 0d0) then
        snkx = 0d0
        snky = 0d0
        snkz = 0d0
      else
        snkx = d3(1)/Ak*signk
        snky = d3(2)/Ak*signk
        snkz = d3(3)/Ak*signk
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(3)%f(i,j,k)%A = Ak
      b%dir(3)%f(i,j,k)%n = [ snkx, snky, snkz ]
  
    enddo; enddo; enddo
    !$omp end parallel
  
  end subroutine Compute_Norm_Area


  subroutine Compute_Metric_Tensor( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol )
    implicit none
    type(ARES_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(ARES_tensor_3D_type), intent(out) :: M
    type(ARES_vector_3D_type), intent(out) :: dl
    real(R8), intent(out)             :: vol
    ! Local
    integer  :: h
    real(R8) :: det, A(3,3), cofactor(3,3)
    real(R8) :: vx(8), vy(8), vz(8)

    ! Compute average cell dimension. A is M^-1, inverse of metric tensor.
    ! M1 = [ xcs, ycs, zcs ; xet, yet, zet ; xzi, yzi, zzi ].
    
    A(1,:) = N5 % c - N1 % c + &
             N6 % c - N2 % c + &
             N7 % c - N3 % c + &
             N8 % c - N4 % c
  
    A(2,:) = N3 % c - N1 % c + &
             N4 % c - N2 % c + &
             N7 % c - N5 % c + &
             N8 % c - N6 % c
             
    A(3,:) = N2 % c - N1 % c + &
             N4 % c - N3 % c + &
             N6 % c - N5 % c + &
             N8 % c - N7 % c

    A = 0.25d0 * A

    ! Determinant of M^-1.
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    if ( abs(det) == 0d0 ) then
      write(*,'(A90)') ' [ERROR] Metric tensor det=0. Should not happen, but going on with M==I'
      stop

      M % c = 0d0
      do h = 1, 3
        M % c(h,h) = 1d0
      enddo

    else
      
      cofactor(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))
      cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      cofactor(1,3) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))
      cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      cofactor(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))
      cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      cofactor(3,1) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))
      cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      cofactor(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))

      M % c = cofactor / det

    endif

    ! Cell length in i,j,k
    do h = 1, 3
      dl % c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
    enddo

    ! Cell volume computation
    vx(1) = N1%c(1)
    vy(1) = N1%c(2)
    vz(1) = N1%c(3)

    vx(2) = N5%c(1)
    vy(2) = N5%c(2)
    vz(2) = N5%c(3)

    vx(3) = N3%c(1)
    vy(3) = N3%c(2)
    vz(3) = N3%c(3)

    vx(4) = N7%c(1)
    vy(4) = N7%c(2)
    vz(4) = N7%c(3)
  
    vx(5) = N2%c(1)
    vy(5) = N2%c(2)
    vz(5) = N2%c(3)

    vx(6) = N6%c(1)
    vy(6) = N6%c(2)
    vz(6) = N6%c(3)

    vx(7) = N4%c(1)
    vy(7) = N4%c(2)
    vz(7) = N4%c(3)
  
    vx(8) = N8%c(1)
    vy(8) = N8%c(2)
    vz(8) = N8%c(3)
  
    vol = tvol(vx,vy,vz,1,2,3,5) + tvol(vx,vy,vz,2,4,3,8) &
        + tvol(vx,vy,vz,5,8,6,2) + tvol(vx,vy,vz,5,7,8,3) &
        + tvol(vx,vy,vz,5,8,2,3)
  
    if( vol <= 0d0 ) then
      write(*,'(A90)') ' [ERROR] Negative volume'
      stop
    endif

    contains
  
      pure function tvol(vx, vy, vz, i1, i2, i3, i4) result(volume)
        implicit none
        real(R8), intent(in) :: vx(8), vy(8), vz(8)
        integer, intent(in)  :: i1, i2, i3, i4
        real(R8)             :: volume
  
        volume = abs(((vx(i2)-vx(i1))* &
          ((vy(i3)-vy(i1))*(vz(i4)-vz(i1))-(vy(i4)-vy(i1))*(vz(i3)-vz(i1)))+ &
                                  (vy(i2)-vy(i1))* &
          ((vx(i4)-vx(i1))*(vz(i3)-vz(i1))-(vx(i3)-vx(i1))*(vz(i4)-vz(i1)))+ &
                                  (vz(i2)-vz(i1))* &
          ((vx(i3)-vx(i1))*(vy(i4)-vy(i1))-(vx(i4)-vx(i1))*(vy(i3)-vy(i1)))) &
          /6.d0)
  
      end function tvol
  
  end subroutine Compute_Metric_Tensor

  
  subroutine Compute_Yn ( blk, bc, nb, nbound )
    use ARES_Config_Types_m, only: obj_rans
    implicit none
    integer, intent(in)                                  :: nb, nbound
    type(ARES_block_type), dimension(nb), intent(inout)  :: blk
    type(ARES_bc_type), dimension(nbound), intent(in)    :: bc
    ! Local
    real(R8), dimension(nbound) :: facex, facey, facez, facehs
    real(R8) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
    integer :: Fm, Bm, Im, Jm, Km, l, b, nwall
    integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

    nwall=0
    
    ! bc vector processing.
    do l = 1, nbound

      if ( Is_Wall(bc(l) % type) ) then
        
        Bm = bc(l) % b
        Im = bc(l) % i
        Jm = bc(l) % j
        Km = bc(l) % k
        Fm = bc(l) % f
        nwall = nwall+1

        select case(Fm)
          case(1:2)
            i1 = Im - mod(Fm,2)
            j1 = Jm - 1
            k1 = Km - 1
            i2 = Im - mod(Fm,2)
            j2 = Jm
            k2 = Km - 1
            i3 = Im - mod(Fm,2)
            j3 = Jm - 1
            k3 = Km
            i4 = Im - mod(Fm,2)
            j4 = Jm
            k4 = Km
          case(3:4)
            i1 = Im - 1
            j1 = Jm - mod(Fm,2)
            k1 = Km - 1
            i2 = Im
            j2 = Jm - mod(Fm,2)
            k2 = Km - 1
            i3 = Im - 1
            j3 = Jm - mod(Fm,2)
            k3 = Km
            i4 = Im
            j4 = Jm - mod(Fm,2)
            k4 = Km
          case(5:6)
            i1 = Im - 1
            j1 = Jm - 1
            k1 = Km - mod(Fm,2)
            i2 = Im
            j2 = Jm - 1
            k2 = Km - mod(Fm,2)
            i3 = Im - 1
            j3 = Jm
            k3 = Km - mod(Fm,2)
            i4 = Im
            j4 = Jm
            k4 = Km - mod(Fm,2)
        end select

        ! Boundary face nodes.
        x1 = blk(Bm) % node(i1,j1,k1) % c(1)
        x2 = blk(Bm) % node(i2,j2,k2) % c(1)
        x3 = blk(Bm) % node(i3,j3,k3) % c(1)
        x4 = blk(Bm) % node(i4,j4,k4) % c(1)

        y1 = blk(Bm) % node(i1,j1,k1) % c(2)
        y2 = blk(Bm) % node(i2,j2,k2) % c(2)
        y3 = blk(Bm) % node(i3,j3,k3) % c(2)
        y4 = blk(Bm) % node(i4,j4,k4) % c(2)

        z1 = blk(Bm) % node(i1,j1,k1) % c(3)
        z2 = blk(Bm) % node(i2,j2,k2) % c(3)
        z3 = blk(Bm) % node(i3,j3,k3) % c(3)
        z4 = blk(Bm) % node(i4,j4,k4) % c(3)
        
        ! Face center coordinates.
        facex(nwall)  = 0.25d0*( x1 + x2 + x3 + x4 )
        facey(nwall)  = 0.25d0*( y1 + y2 + y3 + y4 )
        facez(nwall)  = 0.25d0*( z1 + z2 + z3 + z4 )
        facehs(nwall) = bc(l) % k_rough
      
      endif
    
    enddo

    ! After bc processing nwall is the viscous face counter.
    ! Compute distance in block b using viscous face coordinates facex, facey, facez.
    do b = 1, nb
      call Wall_Distance_Blk ( blk(b), facex, facey, facez, facehs, nwall )
    enddo

    contains
      
      subroutine Wall_Distance_Blk ( blk, fx, fy, fz, fhs, nwall )
        implicit none
        type(ARES_block_type), intent(inout)       :: blk
        integer, intent(in)                        :: nwall
        real(R8), intent(in), dimension(nwall) :: fx, fy, fz, fhs
        ! Local
        integer :: l, i, j, k
        real(R8) :: center(3), dummy

        !$omp parallel
        !$omp do collapse(3) private(i, j, k, l, center, dummy)
        do k = 1, blk % dim(3)
        do j = 1, blk % dim(2)
        do i = 1, blk % dim(1)

          ! Cell center coordinates.
          center = 1d0/8d0 * ( blk % node(i  ,j  ,k  ) % c + blk % node(i  ,j-1,k  ) % c + &
                              blk % node(i  ,j-1,k-1) % c + blk % node(i  ,j  ,k-1) % c + &
                              blk % node(i-1,j  ,k  ) % c + blk % node(i-1,j-1,k  ) % c + &
                              blk % node(i-1,j-1,k-1) % c + blk % node(i-1,j  ,k-1) % c )

          ! Minimum distance initialization.
          blk % yn (i,j,k) = 1d8

          ! Roughness initialization
          if ( obj_rans%rough ) then
            blk%ks(i,j,k) = 0d0
          end if

          ! Wall face vector processing.
          do l = 1, nwall
            dummy = sqrt( (center(1) - fx(l))**2 + (center(2) - fy(l))**2 + (center(3) - fz(l))**2 )
            if (obj_rans%rough .and. ( dummy < blk%yn(i,j,k)) ) then
              blk%ks(i,j,k) = fhs(l)
            end if
            blk % yn(i,j,k) = min ( dummy, blk % yn(i,j,k) )
          enddo

        enddo; enddo; enddo
        !$omp end parallel

      end subroutine Wall_Distance_Blk

  end subroutine Compute_Yn
  

  subroutine Yn_Connection ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, blks )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm, Is, Js, Ks
    type(ARES_block_type), intent(inout) :: blks
    type(ARES_block_type), intent(inout) :: blkm
    ! Local
    integer :: Ig, Jg, Kg

    ! Ghost cell coordinates
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    blkm % yn(Ig,Jg,Kg) = blks % yn(Is,Js,Ks)
    blkm % ks(Ig,Jg,Kg) = blks % ks(Is,Js,Ks)

  end subroutine Yn_Connection


  subroutine BC_Extrapolate_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blk
    type(ARES_tensor_3D_type), intent(out)    :: Mg(2)
    type(ARES_vector_3D_type), intent(out)    :: dlg(2)
    real(R8), intent(out)                :: volg(2)
    ! Local
    integer :: g
    type(ARES_vector_3D_type), dimension(-gc:gc):: N1, N2, N3, N4, N5, N6, N7, N8

    select case (Fm)
      case(1)
        do g = 0, gc
          N1(g)%c = blk % node (Im+g-1,Jm-1,Km-1) % c
          N2(g)%c = blk % node (Im+g-1,Jm-1,Km  ) % c
          N3(g)%c = blk % node (Im+g-1,Jm  ,Km-1) % c
          N4(g)%c = blk % node (Im+g-1,Jm  ,Km  ) % c
        enddo
        do g = -1, -gc, -1
          N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
          N2(g)%c = 3d0*N2(g+1)%c - 3d0*N2(g+2)%c + N2(g+3)%c
          N3(g)%c = 3d0*N3(g+1)%c - 3d0*N3(g+2)%c + N3(g+3)%c
          N4(g)%c = 3d0*N4(g+1)%c - 3d0*N4(g+2)%c + N4(g+3)%c
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g+1), N2(g+1), N3(g+1), N4(g+1), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km))
          Mg(-g)   = blk % M(Im+g, Jm, Km)
          dlg(-g)  = blk % dl(Im+g, Jm, Km)
          volg(-g) = blk % vol(Im+g, Jm, Km)
        enddo
        
      case(2)
        do g = 0, -gc, -1
          N5(g)%c = blk % node (Im+g,Jm-1,Km-1) % c
          N6(g)%c = blk % node (Im+g,Jm-1,Km  ) % c
          N7(g)%c = blk % node (Im+g,Jm  ,Km-1) % c
          N8(g)%c = blk % node (Im+g,Jm  ,Km  ) % c
        enddo
        do g = 1, gc
          N5(g)%c = 3d0*N5(g-1)%c - 3d0*N5(g-2)%c + N5(g-3)%c
          N6(g)%c = 3d0*N6(g-1)%c - 3d0*N6(g-2)%c + N6(g-3)%c
          N7(g)%c = 3d0*N7(g-1)%c - 3d0*N7(g-2)%c + N7(g-3)%c
          N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km))
          Mg(g)   = blk % M(Im+g, Jm, Km)
          dlg(g)  = blk % dl(Im+g, Jm, Km)
          volg(g) = blk % vol(Im+g, Jm, Km)
        enddo

      case(3)
        if ( blk % dim(2) > 2 ) then
          do g = 0, gc
            N1(g)%c = blk % node (Im-1,Jm+g-1,Km-1) % c
            N2(g)%c = blk % node (Im-1,Jm+g-1,Km  ) % c
            N5(g)%c = blk % node (Im  ,Jm+g-1,Km-1) % c
            N6(g)%c = blk % node (Im  ,Jm+g-1,Km  ) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
            N2(g)%c = 3d0*N2(g+1)%c - 3d0*N2(g+2)%c + N2(g+3)%c
            N5(g)%c = 3d0*N5(g+1)%c - 3d0*N5(g+2)%c + N5(g+3)%c
            N6(g)%c = 3d0*N6(g+1)%c - 3d0*N6(g+2)%c + N6(g+3)%c
            call Compute_Metric_Tensor ( N1(g), N2(g), N1(g+1), N2(g+1), N5(g), N6(g), N5(g+1), N6(g+1), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(-g)   = blk % M(Im, Jm+g, Km)
            dlg(-g)  = blk % dl(Im, Jm+g, Km)
            volg(-g) = blk % vol(Im, Jm+g, Km)
          enddo
        else
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm+g-1,Km-1) % c
            N2(g)%c = blk % node (Im-1,Jm+g-1,Km  ) % c
            N5(g)%c = blk % node (Im  ,Jm+g-1,Km-1) % c
            N6(g)%c = blk % node (Im  ,Jm+g-1,Km  ) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
            N2(g)%c = 2d0*N2(g+1)%c - N2(g+2)%c
            N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
            N6(g)%c = 2d0*N6(g+1)%c - N6(g+2)%c
            call Compute_Metric_Tensor ( N1(g), N2(g), N1(g+1), N2(g+1), N5(g), N6(g), N5(g+1), N6(g+1), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(-g)   = blk % M(Im, Jm+g, Km)
            dlg(-g)  = blk % dl(Im, Jm+g, Km)
            volg(-g) = blk % vol(Im, Jm+g, Km)
          enddo
        endif

      case(4)
        if ( blk % dim(2) > 2 ) then
          do g = 0, -gc, -1
            N3(g)%c = blk % node (Im-1,Jm+g,Km-1) % c
            N4(g)%c = blk % node (Im-1,Jm+g,Km  ) % c
            N7(g)%c = blk % node (Im  ,Jm+g,Km-1) % c
            N8(g)%c = blk % node (Im  ,Jm+g,Km  ) % c
          enddo
          do g = 1, gc
            N3(g)%c = 3d0*N3(g-1)%c - 3d0*N3(g-2)%c + N3(g-3)%c
            N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
            N7(g)%c = 3d0*N7(g-1)%c - 3d0*N7(g-2)%c + N7(g-3)%c
            N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
            call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(g)   = blk % M(Im, Jm+g, Km)
            dlg(g)  = blk % dl(Im, Jm+g, Km)
            volg(g) = blk % vol(Im, Jm+g, Km)
          enddo
        else
          do g = 0, -1, -1
            N3(g)%c = blk % node (Im-1,Jm+g,Km-1) % c
            N4(g)%c = blk % node (Im-1,Jm+g,Km  ) % c
            N7(g)%c = blk % node (Im  ,Jm+g,Km-1) % c
            N8(g)%c = blk % node (Im  ,Jm+g,Km  ) % c
          enddo
          do g = 1, gc
            N3(g)%c = 2d0*N3(g-1)%c - N3(g-2)%c
            N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
            N7(g)%c = 2d0*N7(g-1)%c - N7(g-2)%c
            N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
            call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                         blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km))
            Mg(g)   = blk % M(Im, Jm+g, Km)
            dlg(g)  = blk % dl(Im, Jm+g, Km)
            volg(g) = blk % vol(Im, Jm+g, Km)
          enddo
        endif

      case(5)
        if ( ndir==2 .and. delthe==0d0 .or. ndir==1 ) then
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
            N3(g)%c = 2d0*N3(g+1)%c - N3(g+2)%c
            N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
            N7(g)%c = 2d0*N7(g+1)%c - N7(g+2)%c
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        elseif (ndir==2 .and. delthe/=0d0) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c(1) = N1(0) % c(1) 
            N1(g)%c(2) = N1(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N1(g)%c(3) = N1(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N3(g)%c(1) = N3(0) % c(1) 
            N3(g)%c(2) = N3(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N3(g)%c(3) = N3(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N5(g)%c(1) = N5(0) % c(1) 
            N5(g)%c(2) = N5(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N5(g)%c(3) = N5(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N7(g)%c(1) = N7(0) % c(1) 
            N7(g)%c(2) = N7(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N7(g)%c(3) = N7(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo         
        elseif (ndir==3) then
          do g = 0, gc
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 3d0*N1(g+1)%c - 3d0*N1(g+2)%c + N1(g+3)%c
            N3(g)%c = 3d0*N3(g+1)%c - 3d0*N3(g+2)%c + N3(g+3)%c
            N5(g)%c = 3d0*N5(g+1)%c - 3d0*N5(g+2)%c + N5(g+3)%c
            N7(g)%c = 3d0*N7(g+1)%c - 3d0*N7(g+2)%c + N7(g+3)%c
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      case(6)
        if ( ndir==2 .and. delthe==0d0 .or. ndir==1 ) then
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c = 2d0*N2(g-1)%c - N2(g-2)%c
            N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
            N6(g)%c = 2d0*N6(g-1)%c - N6(g-2)%c
            N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        elseif (ndir==2 .and. delthe/=0d0) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c(1) = N2(0) % c(1) 
            N2(g)%c(2) = N2(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N2(g)%c(3) = N2(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N4(g)%c(1) = N4(0) % c(1) 
            N4(g)%c(2) = N4(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N4(g)%c(3) = N4(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N6(g)%c(1) = N6(0) % c(1) 
            N6(g)%c(2) = N6(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N6(g)%c(3) = N6(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N8(g)%c(1) = N8(0) % c(1) 
            N8(g)%c(2) = N8(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N8(g)%c(3) = N8(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo         
        elseif (ndir==3) then
          do g = 0, -gc, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c = 3d0*N2(g-1)%c - 3d0*N2(g-2)%c + N2(g-3)%c
            N4(g)%c = 3d0*N4(g-1)%c - 3d0*N4(g-2)%c + N4(g-3)%c
            N6(g)%c = 3d0*N6(g-1)%c - 3d0*N6(g-2)%c + N6(g-3)%c
            N8(g)%c = 3d0*N8(g-1)%c - 3d0*N8(g-2)%c + N8(g-3)%c
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g) )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      end select

  end subroutine BC_Extrapolate_Metrics


  subroutine BC_Symmetry_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout)   :: blk
    type(ARES_tensor_3D_type), intent(out) :: Mg(2)
    type(ARES_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    ! Local
    integer :: g

    do g = 1, gc
      select case (Fm)
        case(1)
          blk % M(Im-g,Jm,Km)   = blk % M(Im+g-1,Jm,Km)
          blk % dl(Im-g,Jm,Km)  = blk % dl(Im+g-1,Jm,Km)
          blk % vol(Im-g,Jm,Km) = blk % vol(Im+g-1,Jm,Km)
          Mg(g)   = blk % M(Im-g,Jm,Km)
          dlg(g)  = blk % dl(Im-g,Jm,Km)
          volg(g) = blk % vol(Im-g,Jm,Km)
        case(2)
          blk % M(Im+g,Jm,Km)   = blk % M(Im-g+1,Jm,Km)
          blk % dl(Im+g,Jm,Km)  = blk % dl(Im-g+1,Jm,Km)
          blk % vol(Im+g,Jm,Km) = blk % vol(Im-g+1,Jm,Km)
          Mg(g)   = blk % M(Im+g,Jm,Km)
          dlg(g)  = blk % dl(Im+g,Jm,Km)
          volg(g) = blk % vol(Im+g,Jm,Km)
        case(3)
          blk % M(Im,Jm-g,Km)   = blk % M(Im,Jm+g-1,Km)
          blk % dl(Im,Jm-g,Km)  = blk % dl(Im,Jm+g-1,Km)
          blk % vol(Im,Jm-g,Km) = blk % vol(Im,Jm+g-1,Km)
          Mg(g)   = blk % M(Im,Jm-g,Km)
          dlg(g)  = blk % dl(Im,Jm-g,Km)
          volg(g) = blk % vol(Im,Jm-g,Km)
        case(4)
          blk % M(Im,Jm+g,Km)   = blk % M(Im,Jm-g+1,Km)
          blk % dl(Im,Jm+g,Km)  = blk % dl(Im,Jm-g+1,Km)
          blk % vol(Im,Jm+g,Km) = blk % vol(Im,Jm-g+1,Km)
          Mg(g)   = blk % M(Im,Jm+g,Km)
          dlg(g)  = blk % dl(Im,Jm+g,Km)
          volg(g) = blk % vol(Im,Jm+g,Km)
        case(5)
          blk % M(Im,Jm,Km-g)   = blk % M(Im,Jm,Km+g-1)
          blk % dl(Im,Jm,Km-g)  = blk % dl(Im,Jm,Km+g-1)
          blk % vol(Im,Jm,Km-g) = blk % vol(Im,Jm,Km+g-1)
          Mg(g)   = blk % M(Im,Jm,Km-g)
          dlg(g)  = blk % dl(Im,Jm,Km-g)
          volg(g) = blk % vol(Im,Jm,Km-g)
        case(6)
          blk % M(Im,Jm,Km+g)   = blk % M(Im,Jm,Km-g+1)
          blk % dl(Im,Jm,Km+g)  = blk % dl(Im,Jm,Km-g+1)
          blk % vol(Im,Jm,Km+g) = blk % vol(Im,Jm,Km-g+1)
          Mg(g)   = blk % M(Im,Jm,Km+g)
          dlg(g)  = blk % dl(Im,Jm,Km+g)
          volg(g) = blk % vol(Im,Jm,Km+g)
      end select
    enddo

  end subroutine BC_Symmetry_Metrics


  subroutine BC_Connect_Metrics ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg )
    implicit none
    integer, intent(in)                  :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(ARES_block_type), intent(in)    :: blks
    type(ARES_block_type), intent(inout) :: blkm
    type(ARES_tensor_3D_type), intent(out) :: Mg(2)
    type(ARES_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    ! Local
    integer      :: g, Is1, Js1, Ks1
    integer      :: II, JJ, KK, III, JJJ, KKK      
    integer      :: guidem(6,3), guides(6,3), guidem2(6,3), guides2(6,3)
    type(ARES_vector_3D_type), dimension(0:gc) :: N1, N2, N3, N4, N5, N6, N7, N8
    
    ! ---------------------------------------------------------------------------------------------
    ! Preliminary definitions
    guidem(1,1)=-2 ; guidem(1,2)= 0 ; guidem(1,3)= 0 ; guidem(2,1)=1 ; guidem(2,2)=0 ; guidem(2,3)=0
    guidem(3,1)= 0 ; guidem(3,2)=-2 ; guidem(3,3)= 0 ; guidem(4,1)=0 ; guidem(4,2)=1 ; guidem(4,3)=0
    guidem(5,1)= 0 ; guidem(5,2)=0  ; guidem(5,3)=-2 ; guidem(6,1)=0 ; guidem(6,2)=0 ; guidem(6,3)=1

    guidem2(1,1)=-3 ; guidem2(1,2)= 0 ; guidem2(1,3)= 0 ; guidem2(2,1)=2 ; guidem2(2,2)=0 ; guidem2(2,3)=0
    guidem2(3,1)= 0 ; guidem2(3,2)=-3 ; guidem2(3,3)= 0 ; guidem2(4,1)=0 ; guidem2(4,2)=2 ; guidem2(4,3)=0
    guidem2(5,1)= 0 ; guidem2(5,2)= 0 ; guidem2(5,3)=-3 ; guidem2(6,1)=0 ; guidem2(6,2)=0 ; guidem2(6,3)=2

    guides(1,1)= 0 ; guides(1,2)= 0 ; guides(1,3)= 0 ; guides(2,1)=-1 ; guides(2,2)= 0 ; guides(2,3)= 0
    guides(3,1)= 0 ; guides(3,2)= 0 ; guides(3,3)= 0 ; guides(4,1)= 0 ; guides(4,2)=-1 ; guides(4,3)= 0
    guides(5,1)= 0 ; guides(5,2)= 0 ; guides(5,3)= 0 ; guides(6,1)= 0 ; guides(6,2)= 0 ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0 ; guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0 ; guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1 ; guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2
    ! ---------------------------------------------------------------------------------------------

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
      case(2)
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(3)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
      case(4)
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(5)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
      case(6)
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
    end select

    do g = 1, gc

      ! connected cell node index
      Is1 = Is + (gc-g)*guides(Fs,1)+(g-1)*guides2(Fs,1)
      Js1 = Js + (gc-g)*guides(Fs,2)+(g-1)*guides2(Fs,2)
      Ks1 = Ks + (gc-g)*guides(Fs,3)+(g-1)*guides2(Fs,3)

      ! compute the indexes in the connected block corresponding to the ghost nodes to update in the local block 
      select case(Fs)
        case(1:2) ! i-faces
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4) ! j-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6) ! k-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      if (Fm == 1) then
          N4(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 2) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 3) then
          N6(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 4) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 5) then
          N7(g)%c = blks%node(ii,jj,kk)%c
      elseif (Fm == 6) then
          N8(g)%c = blks%node(ii,jj,kk)%c
      endif

      ! update the remaining 3 indexes (only necessary in the corners, but for good measure)
      ! node i1-1,i2 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s
          KKK=KK-d12s
        case(3:4)
          III=II-d11s
          JJJ=JJ
          KKK=KK-d12s
        case(5:6)
          III=II-d11s
          JJJ=JJ-d12s
          KKK=KK
      end select
      if (Fm == 1) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d21s
          KKK=KK-d22s
        case(3:4)
          III=II-d21s
          JJJ=JJ
          KKK=KK-d22s
        case(5:6)
          III=II-d21s
          JJJ=JJ-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! node i1-1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s-d21s
          KKK=KK-d12s-d22s
        case(3:4)
          III=II-d11s-d21s
          JJJ=JJ
          KKK=KK-d12s-d22s
        case(5:6)
          III=II-d11s-d21s
          JJJ=JJ-d12s-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 2) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 3) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 4) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 5) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c
      elseif (Fm == 6) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c
      endif

      ! compute metric variables
      select case(Fm)
        case(1)
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g-1), N2(g-1), N3(g-1), N4(g-1), &
                                       blkm%M(Im-g, Jm, Km), blkm%dl(Im-g, Jm, Km), blkm%vol(Im-g, Jm, Km))
          Mg(g)   = blkm % M(Im-g, Jm, Km)
          dlg(g)  = blkm % dl(Im-g, Jm, Km)
          volg(g) = blkm % vol(Im-g, Jm, Km)
        case(2)
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blkm%M(Im+g, Jm, Km), blkm%dl(Im+g, Jm, Km), blkm%vol(Im+g, Jm, Km))
          Mg(g)   = blkm % M(Im+g, Jm, Km)
          dlg(g)  = blkm % dl(Im+g, Jm, Km)
          volg(g) = blkm % vol(Im+g, Jm, Km)
        case(3)
          call Compute_Metric_Tensor ( N1(g), N2(g), N1(g-1), N2(g-1), N5(g), N6(g), N5(g-1), N6(g-1), &
                                       blkm%M(Im, Jm-g, Km), blkm%dl(Im, Jm-g, Km), blkm%vol(Im, Jm-g, Km))
          Mg(g)   = blkm % M(Im, Jm-g, Km)
          dlg(g)  = blkm % dl(Im, Jm-g, Km)
          volg(g) = blkm % vol(Im, Jm-g, Km)
        case(4)
          call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                       blkm%M(Im, Jm+g, Km), blkm%dl(Im, Jm+g, Km), blkm%vol(Im, Jm+g, Km))
          Mg(g)   = blkm % M(Im, Jm+g, Km)
          dlg(g)  = blkm % dl(Im, Jm+g, Km)
          volg(g) = blkm % vol(Im, Jm+g, Km)
        case(5)
          call Compute_Metric_Tensor ( N1(g), N1(g-1), N3(g), N3(g-1), N5(g), N5(g-1), N7(g), N7(g-1), &
                                       blkm%M(Im, Jm, Km-g), blkm%dl(Im, Jm, Km-g), blkm%vol(Im, Jm, Km-g))
          Mg(g)   = blkm % M(Im, Jm, Km-g)
          dlg(g)  = blkm % dl(Im, Jm, Km-g)
          volg(g) = blkm % vol(Im, Jm, Km-g)
        case(6)
          call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                       blkm%M(Im, Jm, Km+g), blkm%dl(Im, Jm, Km+g), blkm%vol(Im, Jm, Km+g))
          Mg(g)   = blkm % M(Im, Jm, Km+g)
          dlg(g)  = blkm % dl(Im, Jm, Km+g)
          volg(g) = blkm % vol(Im, Jm, Km+g)
        end select

    enddo ! End loop over ghost cell           

  end subroutine BC_Connect_Metrics

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function Is_Wall(type) result (ans)
      implicit none
      integer, intent(in) :: type
      logical :: ans

      select case (type)
        case (103, 301, 302)  ! wall BCs
          ans = .true.
        case default
          ans = .false.
      end select
    end function Is_Wall

end module ARES_Lib_Metrics