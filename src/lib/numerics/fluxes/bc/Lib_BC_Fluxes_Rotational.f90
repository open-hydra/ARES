module ARES_Lib_BC_Fluxes_Rotational
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public

contains

  subroutine BC_Rotational_Symmetry_Eul ( Im, Jm, Km, Fm, Blk )

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area
    real(R8) :: Prim(nprim), Flux(nprim)


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    normal = Blk % dir(dir) % f(Face_i,Face_j,Face_k) % n
    area = Blk % dir(dir) % f(Face_i,Face_j,Face_k) % a
    Prim = Blk % P(:,Im,Jm,Km)
    Flux(nv) = Prim(np) * Area * Normal(2)
    ! Residual update
    Blk % r(nv,Im,Jm,Km) = Blk % r(nv,Im,Jm,Km) + modFm2 * Flux(nv)

  end subroutine BC_Rotational_Symmetry_Eul


  subroutine BC_Rotational_Periodic_Eul ( Im, Jm, Km, Fm, Blk )
    use ARES_Lib_Reconstruction, only : Reconstruction
    use ARES_Lib_Metrics, only : delthe
    use ARES_Mod_Riemann
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: Blk
    ! Local
    integer :: Dir, Face_i, Face_j, Face_k, FluxSign
    integer :: Fs, Ks, K1, K2, K3, K4, modFm, Modfm1, ModFm2, ModFm3
    real(R8) :: Normal(3), Area, Prim(nprim), rho, Face_V, Face_W, Fmass, Flux(nprim)
    real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4
    real(R8), dimension(nprim) :: Prim_L, Prim_R
    real(R8) :: v_rot, w_rot, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
    real(R8) :: rho_L, rho_R, F_r, F_u, F_v, F_w, F_e
    real(R8) :: su, sel_L, sel_R


    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! boundary cell primitive variables
    Prim = Blk % P(:,Im,Jm,Km)
    rho = ph2vars( Prim(np), Prim(nh), rho_tab ) 

    if (ndir /= 3) then
      ! NB: for face 5 and 6: nx = 0 ; ny = +- sin(delthe/2) ; nz = cos(delthe/2).
      ! Get velocity at the interface by rotating it of an angle = delthe/2 in the y-z plane since the K=1 cell 
      ! spans from -delthe/2 to delthe/2. The component of (uu, vv, ww) on the normal to the face is w because of the rotation.
      ! For reference, the formula are: uf = u ; vf =  v*nz + w*ny ; wf = -v*ny + w*nz
      Face_V =   Prim(nv) * Normal(3) + Prim(nw) * Normal(2)
      Face_W = - Prim(nv) * Normal(2) + Prim(nw) * Normal(3)

      ! Except from f(np+2) and f(np+3), fluxes from face 5 and 6 are equal and opposite. No need to compute them
      ! (nor the transport of scalar quantities such as roi and mit). Only terms needed: pf*ny (axisymmetry) and fmass*w*ny 
      ! (centrifugal force) in Y direction, -fmass*v*ny (?) in Z direction.
      Fmass  = rho * Area * Prim(nw)     ! written like this because (uu*nx+vv*ny+ww*nz) = w
      Flux = 0d0
      Flux(nv) = Face_V * Fmass + Prim(np) * Area * Normal(2)
      Flux(nw) = Face_W * Fmass

      ! Residual update
      FluxSign = 1 - 2 * mod(Fm,2)
      Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + FluxSign * Flux
    else
      ! Same logic as standard connection
      if (Fm == 5) then
        Fs = 6
        Ks = blk % dim(3)
      elseif (Fm == 6) then
        Fs = 5
        Ks = 1
      else
        error stop ( ' Rotational periodic boundary condition is only for face 5-6' )
      end if

      K1 = modFm*( Ks + guide(Fs,3) ) + modFm1*( Km + guide(Fm,3) )
      K2 = modFm1*Km + modFm*Ks
      K3 = modFm1*Ks + modFm*Km
      K4 = modFm1*(Ks + guide(Fs,3)) + modFm*(Km + guide(Fm,3))

      ! Primitive variables for the stencil
      Prim1 = blk % P(:,Im,Jm,K1)
      Prim2 = blk % P(:,Im,Jm,K2)
      Prim3 = blk % P(:,Im,Jm,K3)
      Prim4 = blk % P(:,Im,Jm,K4)

      ! Rotation
      if (Fm == 5) then
        w_rot =  Prim1(nw) * cos(delthe) - Prim1(nv) * sin(delthe) 
        v_rot =  Prim1(nw) * sin(delthe) + Prim1(nv) * cos(delthe)
        Prim1(nv) = v_rot
        Prim1(nw) = w_rot
        w_rot =  Prim2(nw) * cos(delthe) - Prim2(nv) * sin(delthe) 
        v_rot =  Prim2(nw) * sin(delthe) + Prim2(nv) * cos(delthe)
        Prim2(nv) = v_rot
        Prim2(nw) = w_rot
      elseif (Fm == 6) then
        w_rot =  Prim3(nw) * cos(delthe) + Prim3(nv) * sin(delthe) 
        v_rot = -Prim3(nw) * sin(delthe) + Prim3(nv) * cos(delthe)
        Prim3(nv) = v_rot
        Prim3(nw) = w_rot
        w_rot =  Prim4(nw) * cos(delthe) + Prim4(nv) * sin(delthe) 
        v_rot = -Prim4(nw) * sin(delthe) + Prim4(nv) * cos(delthe)
        Prim4(nv) = v_rot
        Prim4(nw) = w_rot
      end if  

      Dl1 = blk % Dl(Im,Jm,K1) % c(3) * 0.5d0
      Dl2 = blk % Dl(Im,Jm,K2) % c(3) * 0.5d0
      Dl3 = blk % Dl(Im,Jm,K3) % c(3) * 0.5d0
      Dl4 = blk % Dl(Im,Jm,K4) % c(3) * 0.5d0

      Dist23 = Dl2 + Dl3
      Dist12 = Dl1 + Dl2
      Dist34 = Dl3 + Dl4
      
      ! Piecewise-linear MUSCL reconstruction of the L/R interface states (as in the standard connection BC)
      call Reconstruction ( Prim1, Prim2, Prim3, Prim4, Dist12, Dist23, Dist34, Dl2, Dl3, Prim_L, Prim_R )

      ! Auxiliary variables for the L/R states
      rho_L = ph2vars( Prim_L(np), Prim_L(nh), rho_tab )
      rho_R = ph2vars( Prim_R(np), Prim_R(nh), rho_tab )

      call Riemann ( Prim_L(np), Prim_L(nu), Prim_L(nv), Prim_L(nw), Prim_L(nh), &
                     Prim_R(np), Prim_R(nu), Prim_R(nv), Prim_R(nw), Prim_R(nh), &
                     normal(1),normal(2),normal(3), F_r, F_u, F_v, F_w, F_E )

      su = sign ( 0.5d0, F_r )
      Sel_L = 0.5d0 + su
      Sel_R = su - 0.5d0

      ! Fluxes
      Fmass = F_r * area
      Flux(np) = Fmass*( Sel_L - Sel_R )
      Flux(nu) = F_u * Area
      Flux(nv) = F_v * Area
      Flux(nw) = F_w * Area
      Flux(nh) = F_E * Area

      if (model==2) then
        Flux(nt:nprim) = Fmass*(sel_L*Prim_L(nt:nprim)/rho_L - sel_R*Prim_R(nt:nprim)/rho_R)
      end if
    
      ! Residual update
      if (Fm==6) then
        blk % r(:,Im,Jm,Km) = blk % r(:,Im,Jm,Km) + Flux
        ! Flux is in the face-6 frame (θ=delthe). Rotate momentum back to face-5
        ! frame (θ=0) before applying to Ks, so that (v,w) point in the right
        ! Cartesian directions for a cell sitting at θ=0.
        Flux(nv) = ( F_v * cos(delthe) + F_w * sin(delthe) ) * Area !!!!!!!!!!!!!!!!!!!!
        Flux(nw) = (-F_v * sin(delthe) + F_w * cos(delthe) ) * Area !!!!!!!!!!!!!!!!!!!!
        blk % r(:,Im,Jm,Ks) = blk % r(:,Im,Jm,Ks) - Flux
      end if
    end if

  end subroutine BC_Rotational_Periodic_Eul

end module ARES_Lib_BC_Fluxes_Rotational