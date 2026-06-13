module ARES_Lib_BC_Fluxes_Connection
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public

contains

  subroutine BC_Connection_Eul ( Im, Jm, Km, Fm, Blk )
    use ARES_Lib_Reconstruction, only : reconstruction
    use ARES_Mod_Riemann
    use FLINT_Lib_Thermodynamic
    
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: Blk
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dirm, Face_i, Face_j, Face_k
    integer :: I1, J1, K1, I2, J2, K2
    integer :: I3, J3, K3, I4, J4, K4
    real(R8), dimension(nprim) :: Prim1, Prim2, Prim3, Prim4
    real(R8), dimension(nprim) :: Prim_L, Prim_R, Flux
    real(R8) :: Normal(3), Area, Dl1, Dl2, Dl3, Dl4, Dist12, Dist23, Dist34
    real(R8) :: rho_L, rho_R, F_r, F_u, F_v, F_w, F_e, Fmass
    real(R8) :: su, sel_L, sel_R


    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, Dirm, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Stencil construction: 4 points for 2nd order reconstruction from both sides at boundary face
    I1 = modfm*( Im - guide(Fm,1) * gc) + modfm1*( Im + guide(Fm,1) )
    J1 = modFm*( Jm - guide(Fm,2) * gc) + modFm1*( Jm + guide(Fm,2) )
    K1 = modFm*( Km - guide(Fm,3) * gc) + modFm1*( Km + guide(Fm,3) )

    I2 = modFm1*Im + modFm*( Im - guide(Fm,1))
    J2 = modFm1*Jm + modFm*( Jm - guide(Fm,2))
    K2 = modFm1*Km + modFm*( Km - guide(Fm,3))

    I3 = modFm1*( Im - guide(Fm,1)) + modFm*Im
    J3 = modFm1*( Jm - guide(Fm,2)) + modFm*Jm
    K3 = modFm1*( Km - guide(Fm,3)) + modFm*Km

    I4 = modFm1*( Im - guide(Fm,1) * gc) + modFm*(Im + guide(Fm,1))
    J4 = modFm1*( Jm - guide(Fm,2) * gc) + modFm*(Jm + guide(Fm,2))
    K4 = modFm1*( Km - guide(Fm,3) * gc) + modFm*(Km + guide(Fm,3))

    ! Primitive variables for the stencil
    Prim1 = blk % P(:,I1,J1,K1)
    Prim2 = blk % P(:,I2,J2,K2)
    Prim3 = blk % P(:,I3,J3,K3)
    Prim4 = blk % P(:,I4,J4,K4)
    
    Normal = blk % dir(Dirm) % f(Face_i,Face_j,Face_k) % n
    Area = blk % dir(Dirm) % f(Face_i,Face_j,Face_k) % a
    Dl1 = blk % Dl(I1,J1,K1) % c(Dirm) * 0.5d0
    Dl2 = blk % Dl(I2,J2,K2) % c(Dirm) * 0.5d0
    Dl3 = blk % Dl(I3,J3,K3) % c(Dirm) * 0.5d0
    Dl4 = blk % Dl(I4,J4,K4) % c(Dirm) * 0.5d0

    Dist23 = Dl2 + Dl3
    Dist12 = Dl1 + Dl2
    Dist34 = Dl3 + Dl4

    call Reconstruction ( Prim1, Prim2, Prim3, Prim4, Dist12, Dist23, Dist34, Dl2, Dl3, Prim_L, Prim_R )

    ! Auxiliary variables for state (L)
    rho_L = ph2vars( Prim_L(np), Prim_L(nh), rho_tab ) 

    ! Auxiliary variables for state (R)
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

    if (model==2) &
      Flux(nt:nprim) = Fmass*(sel_L*Prim_L(nt:nprim)/rho_L - sel_R*Prim_R(nt:nprim)/rho_R)
    
    ! Residual update
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + Flux * modfm2

  end subroutine BC_Connection_Eul


  subroutine BC_Connection_Visc ( Im, Jm, Km, Fm, Blk, Mg, Pg, Sc, Sct, Prt, Prt_corr )
    use ARES_Lib_RANS
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Diffusive

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_tensor_3D_type), intent(in) :: Mg
    real(kind=8), intent(in)              :: Pg(nprim,6)
    real(R8), intent(in) :: Sc, Sct, Prt
    logical, intent(in)  :: Prt_corr
    type(ARES_block_type), intent(inout)  :: Blk
    ! Local
    integer :: dir, modfm2, Face_i, Face_j, Face_k
    integer :: Ig, Jg, Kg, i
    real(R8) :: Normal(3), Area, Waldis, M(3,3), Rough
    real(R8), dimension(nprim) :: Prim_loc, Prim_ghost, Visc_loc, Visc_ghost
    real(R8), dimension(nprim) :: Visc_ip, Visc_im, Visc_jp, Visc_jm, Visc_kp, Visc_km
    real(R8), dimension(nprim) :: Visc_ghost3, Visc_ghost4, Visc_ghost5, Visc_ghost6
    real(R8) :: Prim(nprim), Gradient_loc(nprim,3), Gradient_ghost(nprim,3), Flux_loc(nprim), Flux_ghost(nprim)


    ! Boundary face index
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    modfm2 = 1 - 2 * mod (Fm,2)

    ! Metric stuff
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)
    Normal = blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    area = blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    M = 0.5d0 * ( Blk % M(Im,Jm,Km) % c + Blk % M(Ig,Jg,Kg) % c )
    Waldis = 0.5d0 * ( Blk % yn(Im,Jm,Km) + Blk % yn(Ig,Jg,Kg) )
    Rough = 0.5d0 * ( Blk%ks(Im,Jm,Km) + Blk%ks(Ig,Jg,Kg) )

    ! Primitive/auxiliary variables and residual in the 2 connected cells
    Prim_loc = Blk % P(:,Im,Jm,Km)
    Prim_ghost = Pg (:,1)

    ! stencil building
    call Visc_Variables ( Prim_loc, Visc_loc )
    call Visc_Variables ( Blk % P (:,Im+1,Jm,Km), Visc_ip )
    call Visc_Variables ( Blk % P (:,Im-1,Jm,Km), Visc_im )
    call Visc_Variables ( Blk % P (:,Im,Jm+1,Km), Visc_jp )
    call Visc_Variables ( Blk % P (:,Im,Jm-1,Km), Visc_jm )
    call Visc_Variables ( Blk % P (:,Im,Jm,Km+1), Visc_kp )
    call Visc_Variables ( Blk % P (:,Im,Jm,Km-1), Visc_km )
    call Visc_Variables ( Prim_ghost, Visc_ghost )
    call Visc_Variables ( Pg(:,3), Visc_ghost3 )
    call Visc_Variables ( Pg(:,4), Visc_ghost4 )
    call Visc_Variables ( Pg(:,5), Visc_ghost5 )
    call Visc_Variables ( Pg(:,6), Visc_ghost6 )

    ! Gradient computations
    select case ( Fm )
      case ( 1 : 2 )
        Gradient_loc (:,1) = ( Visc_ghost - Visc_loc ) * modfm2
        Gradient_loc (:,2) = ( Visc_jp - Visc_jm ) * 0.5d0
        Gradient_loc (:,3) = ( Visc_kp - Visc_km ) * 0.5d0
        Gradient_ghost (:,1) = ( Visc_loc - Visc_ghost ) * (-modfm2)
        Gradient_ghost (:,2) = ( Visc_ghost4 - Visc_ghost3 ) * 0.5d0
        Gradient_ghost (:,3) = ( Visc_ghost6 - Visc_ghost5 ) * 0.5d0
      case ( 3 : 4 )
        Gradient_loc (:,1) = ( Visc_ip - Visc_im ) * 0.5d0
        Gradient_loc (:,2) = ( Visc_ghost - Visc_loc ) * modfm2
        Gradient_loc (:,3) = ( Visc_kp - Visc_km ) * 0.5d0
        Gradient_ghost (:,1) = ( Visc_ghost4 - Visc_ghost3 ) * 0.5d0
        Gradient_ghost (:,2) = ( Visc_loc - Visc_ghost ) * (-modfm2)
        Gradient_ghost (:,3) = ( Visc_ghost6 - Visc_ghost5 ) * 0.5d0
      case(5:6)
        Gradient_loc (:,1) = ( Visc_ip - Visc_im ) * 0.5d0
        Gradient_loc (:,2) = ( Visc_jp - Visc_jm ) * 0.5d0
        Gradient_loc (:,3) = ( Visc_ghost - Visc_loc ) * modfm2
        Gradient_ghost (:,1) = ( Visc_ghost4 - Visc_ghost3 ) * 0.5d0
        Gradient_ghost (:,2) = ( Visc_ghost6 - Visc_ghost5 ) * 0.5d0
        Gradient_ghost (:,3) = ( Visc_loc - Visc_ghost ) * (-modfm2)
    end select

    Gradient_loc   = matmul ( Gradient_loc, M )
    Gradient_ghost = matmul ( Gradient_ghost, M )

    Prim = ( Prim_loc + Prim_ghost ) * 0.5d0
    call Compute_Diffusive_Flux ( Prim, Gradient_loc, Area, Normal, Waldis, Rough, Flux_loc, Sc, Sct, Prt, Prt_corr )
    call Compute_Diffusive_Flux ( Prim, Gradient_ghost, Area, Normal, Waldis, Rough, Flux_ghost, Sc, Sct, Prt, Prt_corr )

    ! Residual update
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) - 0.5d0 * modfm2 * (Flux_loc + Flux_ghost)
    
    contains

      subroutine Visc_Variables ( Prim, Visc )

        implicit none
        real(R8), intent(in) :: Prim (nprim)
        real(R8), intent(out) :: Visc (nprim)
        ! Local
        real(R8) :: rho, T

        T = ph2vars( Prim(np), Prim(nh), T_tab ) 
        Visc(1:np) = 0.0_R8
        Visc(nu:nw) = Prim(nu:nw)
        Visc(nh) = T
        if (nprim>nh) then
          rho = ph2vars( Prim(np), Prim(nh), rho_tab ) 
          Visc(nt:nprim) = Prim(nt:nprim) / rho
        endif

      end subroutine Visc_Variables
        
  end subroutine BC_Connection_Visc

end module ARES_Lib_BC_Fluxes_Connection