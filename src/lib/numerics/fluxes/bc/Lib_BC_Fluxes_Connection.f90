module ARES_Lib_BC_Fluxes_Connection
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  public

contains

  !> Inviscid flux of a face between the cell (Im,Jm,Km) and its ghost across
  !> face Fm, from the four-cell stencil across the face; returns the share of
  !> the residual of the cell, dR (the caller adds it).
  subroutine BC_Connection_Eul ( Im, Jm, Km, Fm, Blk, dR )
    use ARES_Lib_Reconstruction, only : reconstruction
    use ARES_Mod_Riemann
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(in) :: Blk
    real(R8), intent(out) :: dR(nprim)
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

    ! Share of the residual of the cell
    dR = Flux * modfm2

  end subroutine BC_Connection_Eul


  !> Viscous flux of a face between the cell (Im,Jm,Km) and its ghost across
  !> face Fm; returns the share of the residual of the cell, dR (the caller adds it).
  subroutine BC_Connection_Visc ( Im, Jm, Km, Fm, Blk, Mg, Pg, Sc, Sct, Prt, Prt_corr, dR )
    use ARES_Lib_RANS
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_Diffusive

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_tensor_3D_type), intent(in) :: Mg
    real(kind=8), intent(in)              :: Pg(nprim,6)
    real(R8), intent(in) :: Sc, Sct, Prt
    logical, intent(in)  :: Prt_corr
    type(ARES_block_type), intent(in)     :: Blk
    real(R8), intent(out)                 :: dR(nprim)
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
    ! Metric tensor at the connection face: the mean of the two inverse metrics
    ! (cell edge vectors), inverted, so that the one-index difference across the
    ! face is divided by the centre-to-centre distance. Averaging M itself divides
    ! it by the mean of 1/h instead: with a large size jump across the connection
    ! (86:1 gives 22x) the normal gradient is far too steep and the explicit
    ! diffusion number on the big cell exceeds 1, although its own dt is fine.
    M = Face_Metric ( Blk % M(Im,Jm,Km) % c, Blk % M(Ig,Jg,Kg) % c )
    Waldis = 0.5d0 * ( Blk % yn(Im,Jm,Km) + Blk % yn(Ig,Jg,Kg) )
    Rough = 0.5d0 * ( Blk%ks(Im,Jm,Km) + Blk%ks(Ig,Jg,Kg) )

    ! Primitive/auxiliary variables and residual in the 2 connected cells
    Prim_loc = Blk % P(:,Im,Jm,Km)
    Prim_ghost = Pg (:,1)

    ! stencil building: the viscous variables of the cells of this block (Viscous_Variables), and those of the cells of
    ! the other block (Pg) computed here
    Visc_loc = Blk % Visc(:,Im,Jm,Km)
    Visc_ip = Blk % Visc(:,Im+1,Jm,Km)
    Visc_im = Blk % Visc(:,Im-1,Jm,Km)
    Visc_jp = Blk % Visc(:,Im,Jm+1,Km)
    Visc_jm = Blk % Visc(:,Im,Jm-1,Km)
    Visc_kp = Blk % Visc(:,Im,Jm,Km+1)
    Visc_km = Blk % Visc(:,Im,Jm,Km-1)
    call Visc_State ( Prim_ghost, Visc_ghost )
    call Visc_State ( Pg(:,3), Visc_ghost3 )
    call Visc_State ( Pg(:,4), Visc_ghost4 )
    call Visc_State ( Pg(:,5), Visc_ghost5 )
    call Visc_State ( Pg(:,6), Visc_ghost6 )

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

    ! Share of the residual of the cell
    dR = - 0.5d0 * modfm2 * (Flux_loc + Flux_ghost)

  end subroutine BC_Connection_Visc


  !> Viscous flux of a face between blocks (101, 201) with diffusive-metric = stencil: one face gradient, as on an
  !> interior face (Diffusive_Flux), with the cells in the order of growing index of this block (the cell and its ghost)
  !> and the face metric of the stencil (ARES_Lib_Stencil_Metrics), whose ghost cells are those of Pg(:,3:6). Where the
  !> index directions of the two blocks agree across the face and their nodes have the same bits, the operations are those
  !> of the same face in the unsplit mesh. Returns the share of the residual of the cell, dR (the caller adds it).
  subroutine BC_Connection_Visc_Stencil ( Im, Jm, Km, Fm, Blk, Pg, Sc, Sct, Prt, Prt_corr, dR )
    use ARES_Lib_Diffusive, only: Diffusive_Flux, Visc_State
    implicit none
    integer, intent(in)  :: Im, Jm, Km, Fm
    real(R8), intent(in) :: Sc, Sct, Prt
    logical, intent(in)  :: Prt_corr
    type(ARES_block_type), intent(in), target :: Blk
    real(R8), intent(in)              :: Pg(nprim,6)
    real(R8), intent(out)             :: dR(nprim)
    ! Local
    integer :: dir, Face_i, Face_j, Face_k, Ig, Jg, Kg, a, b, c
    integer :: e(3,3)
    real(R8), dimension(nprim) :: Rl, Rr
    real(R8) :: Vg(nprim,6)                      ! viscous variables of the cells of Pg (the other block)
    integer :: g
    real(R8), pointer, contiguous :: Mg(:,:,:)   ! metrics per group with a wall ghost cell in the stencil, else absent

    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)
    select case ( dir )
      case (1) ; a = 1 ; b = 2 ; c = 3
      case (2) ; a = 2 ; b = 1 ; c = 3
      case default ; a = 3 ; b = 1 ; c = 2
    end select
    e = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    Rl = 0d0
    Rr = 0d0

    ! viscous variables: of the cells of this block, those of Viscous_Variables; of the cells of the other block, here
    do g = 1, 6
      if ( g /= 2 ) call Visc_State ( Pg(:,g), Vg(:,g) )
    enddo
    nullify ( Mg )
    if ( Blk % dir(dir) % ig(Face_i,Face_j,Face_k) > 0 ) Mg => Blk % dir(dir) % Mg(:,:,:,Blk % dir(dir) % ig(Face_i,Face_j,Face_k))

    associate ( f => Blk % dir(dir) % f(Face_i,Face_j,Face_k), Mf => Blk % dir(dir) % Mf(Face_i,Face_j,Face_k) % c )
    if ( mod(Fm,2) == 0 ) then
      ! the cell is on the low-index side of the face: cell 1 = cell, cell 2 = ghost
      call Diffusive_Flux ( f % N, f % A, Blk % yn(Im,Jm,Km), Blk % yn(Ig,Jg,Kg),                     &
                            Blk % ks(Im,Jm,Km), Blk % ks(Ig,Jg,Kg),                                 &
                            Blk % P(:,Im,Jm,Km), Pg(:,1),                                           &
                            Blk % Visc(:,Im,Jm,Km), Vg(:,1),                                        &
                            Blk % Visc(:,Im-e(1,b),Jm-e(2,b),Km-e(3,b)),                            &
                            Blk % Visc(:,Im+e(1,b),Jm+e(2,b),Km+e(3,b)),                            &
                            Vg(:,3), Vg(:,4),                                                       &
                            Blk % Visc(:,Im-e(1,c),Jm-e(2,c),Km-e(3,c)),                            &
                            Blk % Visc(:,Im+e(1,c),Jm+e(2,c),Km+e(3,c)),                            &
                            Vg(:,5), Vg(:,6),                                                       &
                            Mf, Rl, Rr, a, b, c, Sc, Sct, Prt, Prt_corr, Mg = Mg )
      dR = Rl
    else
      ! cell 1 = ghost, cell 2 = cell
      call Diffusive_Flux ( f % N, f % A, Blk % yn(Ig,Jg,Kg), Blk % yn(Im,Jm,Km),                     &
                            Blk % ks(Ig,Jg,Kg), Blk % ks(Im,Jm,Km),                                 &
                            Pg(:,1), Blk % P(:,Im,Jm,Km),                                           &
                            Vg(:,1), Blk % Visc(:,Im,Jm,Km),                                        &
                            Vg(:,3), Vg(:,4),                                                       &
                            Blk % Visc(:,Im-e(1,b),Jm-e(2,b),Km-e(3,b)),                            &
                            Blk % Visc(:,Im+e(1,b),Jm+e(2,b),Km+e(3,b)),                            &
                            Vg(:,5), Vg(:,6),                                                       &
                            Blk % Visc(:,Im-e(1,c),Jm-e(2,c),Km-e(3,c)),                            &
                            Blk % Visc(:,Im+e(1,c),Jm+e(2,c),Km+e(3,c)),                            &
                            Mf, Rl, Rr, a, b, c, Sc, Sct, Prt, Prt_corr, Mg = Mg )
      dR = Rr
    end if
    end associate

  end subroutine BC_Connection_Visc_Stencil

end module ARES_Lib_BC_Fluxes_Connection