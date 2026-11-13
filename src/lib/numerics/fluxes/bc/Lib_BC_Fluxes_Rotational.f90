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
    use ARES_Lib_Metrics, only : cos_delthe, sin_delthe
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
        w_rot =  Prim1(nw) * cos_delthe - Prim1(nv) * sin_delthe 
        v_rot =  Prim1(nw) * sin_delthe + Prim1(nv) * cos_delthe
        Prim1(nv) = v_rot
        Prim1(nw) = w_rot
        w_rot =  Prim2(nw) * cos_delthe - Prim2(nv) * sin_delthe 
        v_rot =  Prim2(nw) * sin_delthe + Prim2(nv) * cos_delthe
        Prim2(nv) = v_rot
        Prim2(nw) = w_rot
      elseif (Fm == 6) then
        w_rot =  Prim3(nw) * cos_delthe + Prim3(nv) * sin_delthe 
        v_rot = -Prim3(nw) * sin_delthe + Prim3(nv) * cos_delthe
        Prim3(nv) = v_rot
        Prim3(nw) = w_rot
        w_rot =  Prim4(nw) * cos_delthe + Prim4(nv) * sin_delthe 
        v_rot = -Prim4(nw) * sin_delthe + Prim4(nv) * cos_delthe
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
        Flux(nv) = ( F_v * cos_delthe + F_w * sin_delthe ) * Area !!!!!!!!!!!!!!!!!!!!
        Flux(nw) = (-F_v * sin_delthe + F_w * cos_delthe ) * Area !!!!!!!!!!!!!!!!!!!!
        blk % r(:,Im,Jm,Ks) = blk % r(:,Im,Jm,Ks) - Flux
      end if
    end if

  end subroutine BC_Rotational_Periodic_Eul


  !> Viscous flux on a face (5 or 6) of type 200: the 2-D axisymmetric wedge (one cell in k) and the 3-D sector.
  !>
  !> The k ghost of the cell is the cell rotated by delthe about the x axis
  !> (Ghost_Axysimmetry), so the face is the face between the cell and that
  !> image, and its flux is computed as on an interior face (Diffusive_Flux):
  !> the k difference across the face, the i and j differences averaged over
  !> the row of the cell and its rotated image. The image row is rotated here
  !> with the same formula, because the k ghosts of the i and j ghost cells (the
  !> corners) are not filled. Face metric: with diffusive-metric = stencil the
  !> metric of the same differences of the positions, the image row the row of
  !> the source turned (Rotational_Face_Metric, ARES_Lib_Stencil_Metrics, stored on
  !> the face of the block, one per group of quantities with a wall ghost cell in
  !> the stencil), exact for the linear fields that the periodicity allows; with
  !> mean, Mf, the inverse of the mean of the two inverse metrics
  !> of the cell and its image, as on the connections, computed once at setup
  !> (Rotational_Periodic_Face_Metric, stored in the record: bc%rM).
  !>
  !> On the wedge only the momentum is updated. The flow there is axisymmetric by
  !> construction, so, as for the inviscid flux above, the mass, energy and
  !> scalar fluxes of faces 5 and 6 are equal and opposite and cancel; the
  !> momentum fluxes do not. Without swirl their sum is the normal viscous
  !> stress of the faces, tau_thetatheta = 2 mu u_r/r + lambda div(u), acting
  !> along -r: the -tau_thetatheta/r of the axisymmetric equations, next to the
  !> +p/r of the inviscid flux. With swirl the faces also carry tau_rtheta and
  !> tau_xtheta.
  !>
  !> On a 3-D sector (ndir = 3) the image of the cell next to face 6 is the cell
  !> next to face 5 (the source, Ks) turned forward by the sector angle, and the
  !> face is computed once, by the record of face 6, as the inviscid flux: the
  !> image row is the source row turned, the image metric the source metric
  !> turned (its rows are vectors), wall distance and roughness are the mean of
  !> the cell and the source, as on an interior face; the cell gets -F, the
  !> source the same flux with its momentum turned back. A translated sector
  !> (delthe = 0) is the same without the turn. Not handled: a planar 2-D mesh
  !> with type 200 on k (delthe = 0), where the image is the cell itself and the
  !> fluxes of faces 5 and 6 are the same and cancel.
  subroutine BC_Rotational_Periodic_Visc ( Im, Jm, Km, Fm, Blk, Mf, Sc, Sct, Prt, Prt_corr )
    use ARES_Lib_Metrics, only : delthe, cos_delthe, sin_delthe
    use ARES_Lib_Diffusive, only : Compute_Diffusive_Flux, Diffusive_Flux
    use ARES_Config_Types_m, only : obj_space_scheme
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout), target :: Blk
    real(R8), intent(in) :: Mf(3,3)      ! face metric of mean (Rotational_Periodic_Face_Metric)
    real(R8), intent(in) :: Sc, Sct, Prt
    logical, intent(in)  :: Prt_corr
    ! Local
    integer :: Dir, Face_i, Face_j, Face_k, Kg, Ks, modfm2, v
    real(R8) :: Normal(3), Area, Prim(nprim), Flux(nprim), Gradient(nprim,3), g1, g2, g3
    real(R8), dimension(nprim) :: Rot_im, Rot_ip, Rot_jm, Rot_jp, Rl, Rr
    real(R8), pointer, contiguous :: Mg(:,:,:)   ! stencil: metrics per group with a wall ghost cell in the stencil, else absent

    ! No wedge on a 2-D mesh: delthe exactly zero, written as two ordered comparisons (same result, also for -0 and NaN)
    if (ndir < 3 .and. (delthe >= 0d0 .and. delthe <= 0d0)) return
    ! 3-D sector: the record of face 6 computes the face for both cells
    if (ndir == 3 .and. Fm /= 6) return

    call Face_Index ( Fm, Dir, Im, Jm, Km, Face_i, Face_j, Face_k )
    modfm2 = 1 - 2 * mod(Fm,2)     ! +1 on face 6 (image at k+1), -1 on face 5 (image at k-1)
    Kg = Km - guide(Fm,3)          ! k of the image of the cell (its ghost)
    Ks = Km                        ! k of the source of the image: the cell itself on a wedge,
    if (ndir == 3) Ks = Blk % dim(3) - Km + 1   ! the cell next to the opposite k face on a sector

    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area   = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! The viscous variables of the cells (Viscous_Variables); those of the image row, turned
    Rot_im = Rotated ( Blk % Visc(:,Im-1,Jm,Ks) )
    Rot_ip = Rotated ( Blk % Visc(:,Im+1,Jm,Ks) )
    Rot_jm = Rotated ( Blk % Visc(:,Im,Jm-1,Ks) )
    Rot_jp = Rotated ( Blk % Visc(:,Im,Jm+1,Ks) )

    if ( obj_space_scheme % stencil_metric ) then
      ! diffusive-metric = stencil: one face gradient as on an interior face, with the cells in the order of growing k (on
      ! face 6 the cell, then its image; on face 5 the image, then the cell; the rows along the face in either order, as
      ! their differences are averaged), the metric of the face stored at setup and, with a wall ghost cell in the stencil,
      ! the metrics per group. The image has the wall distance and roughness of its source, as below.
      nullify ( Mg )
      if ( Blk % dir(Dir) % ig(Face_i,Face_j,Face_k) > 0 ) &
        Mg => Blk % dir(Dir) % Mg(:,:,:,Blk % dir(Dir) % ig(Face_i,Face_j,Face_k))
      Rl = 0d0
      Rr = 0d0
      if ( Fm == 6 ) then
        call Diffusive_Flux ( Normal, Area, Blk % yn(Im,Jm,Km), Blk % yn(Im,Jm,Ks),                          &
                              Blk % ks(Im,Jm,Km), Blk % ks(Im,Jm,Ks),                                      &
                              Blk % P(:,Im,Jm,Km), Blk % P(:,Im,Jm,Kg),                                    &
                              Blk % Visc(:,Im,Jm,Km), Blk % Visc(:,Im,Jm,Kg),                              &
                              Blk % Visc(:,Im-1,Jm,Km), Blk % Visc(:,Im+1,Jm,Km), Rot_im, Rot_ip,          &
                              Blk % Visc(:,Im,Jm-1,Km), Blk % Visc(:,Im,Jm+1,Km), Rot_jm, Rot_jp,          &
                              Blk % dir(Dir) % Mf(Face_i,Face_j,Face_k) % c, Rl, Rr, 3, 1, 2,              &
                              Sc, Sct, Prt, Prt_corr, Mg = Mg )
        Flux = -Rl                     ! Res1 (the cell) takes -flux
      else
        call Diffusive_Flux ( Normal, Area, Blk % yn(Im,Jm,Ks), Blk % yn(Im,Jm,Km),                          &
                              Blk % ks(Im,Jm,Ks), Blk % ks(Im,Jm,Km),                                      &
                              Blk % P(:,Im,Jm,Kg), Blk % P(:,Im,Jm,Km),                                    &
                              Blk % Visc(:,Im,Jm,Kg), Blk % Visc(:,Im,Jm,Km),                              &
                              Rot_im, Rot_ip, Blk % Visc(:,Im-1,Jm,Km), Blk % Visc(:,Im+1,Jm,Km),          &
                              Rot_jm, Rot_jp, Blk % Visc(:,Im,Jm-1,Km), Blk % Visc(:,Im,Jm+1,Km),          &
                              Blk % dir(Dir) % Mf(Face_i,Face_j,Face_k) % c, Rl, Rr, 3, 1, 2,              &
                              Sc, Sct, Prt, Prt_corr, Mg = Mg )
        Flux = Rr                      ! Res2 (the cell) takes +flux
      end if
    else
      ! Gradients in computational space: k across the face, i and j along it
      Gradient(:,3) = ( Blk % Visc(:,Im,Jm,Kg) - Blk % Visc(:,Im,Jm,Km) ) * modfm2
      Gradient(:,1) = ( Blk % Visc(:,Im+1,Jm,Km) - Blk % Visc(:,Im-1,Jm,Km) + Rot_ip - Rot_im ) * 0.25d0
      Gradient(:,2) = ( Blk % Visc(:,Im,Jm+1,Km) - Blk % Visc(:,Im,Jm-1,Km) + Rot_jp - Rot_jm ) * 0.25d0

      ! ... and in physical space
      do v = 1, nprim
        g1 = Gradient(v,1) ; g2 = Gradient(v,2) ; g3 = Gradient(v,3)
        Gradient(v,1) = g1*Mf(1,1) + g2*Mf(2,1) + g3*Mf(3,1)
        Gradient(v,2) = g1*Mf(1,2) + g2*Mf(2,2) + g3*Mf(3,2)
        Gradient(v,3) = g1*Mf(1,3) + g2*Mf(2,3) + g3*Mf(3,3)
      end do

      ! The image has the wall distance and roughness of its source: the face takes their mean with the cell, as an
      ! interior face (on a wedge the source is the cell)
      Prim = 0.5d0 * ( Blk % P(:,Im,Jm,Km) + Blk % P(:,Im,Jm,Kg) )
      call Compute_Diffusive_Flux ( Prim, Gradient, Area, Normal, &
                                    0.5d0 * ( Blk % yn(Im,Jm,Km) + Blk % yn(Im,Jm,Ks) ), &
                                    0.5d0 * ( Blk % ks(Im,Jm,Km) + Blk % ks(Im,Jm,Ks) ), &
                                    Flux, Sc, Sct, Prt, Prt_corr )
    end if

    if (ndir == 3) then
      ! Residual update of the two cells of the face: the cell next to face 6, and the source next to face 5 with the
      ! momentum turned back, as the inviscid flux
      Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) - Flux
      g2 = Flux(nv) ; g3 = Flux(nw)
      Flux(nv) =  g2 * cos_delthe + g3 * sin_delthe
      Flux(nw) = -g2 * sin_delthe + g3 * cos_delthe
      Blk % r(:,Im,Jm,Ks) = Blk % r(:,Im,Jm,Ks) + Flux
    else
      ! Residual update: momentum only (see above)
      Blk % r(nu:nw,Im,Jm,Km) = Blk % r(nu:nw,Im,Jm,Km) - modfm2 * Flux(nu:nw)
    end if

    contains

      !> The same values in the image row: v and w rotated as in Ghost_Axysimmetry
      function Rotated ( Visc ) result ( Rot )
        implicit none
        real(R8), intent(in) :: Visc(nprim)
        real(R8) :: Rot(nprim)

        Rot = Visc
        Rot(nv) = Visc(nv) * cos_delthe + guide(Fm,3) * Visc(nw) * sin_delthe
        Rot(nw) = Visc(nv) * sin_delthe * (-guide(Fm,3)) + Visc(nw) * cos_delthe

      end function Rotated

  end subroutine BC_Rotational_Periodic_Visc


  !> Face metric of the viscous flux of a face of type 200 (BC_Rotational_Periodic_Visc) with diffusive-metric =
  !> mean (with stencil the face takes Rotational_Face_Metric of ARES_Lib_Stencil_Metrics, which takes
  !> from this one the columns of the index directions with one cell): the inverse of the mean of
  !> the two inverse metrics (Face_Metric) of the cell and of its image. It depends on the metrics only, so
  !> Setup_Metrics computes it once, for the records of type 200 of the blocks of the rank, and stores it in the record
  !> (bc%rM). On a wedge the metric of the image is that of the ghost cell, the exact rotation of the cell
  !> (BC_Extrapolate_Metrics); on a 3-D sector it is the metric of the source (the cell next to the opposite k face)
  !> turned. Zero where the flux is not computed: a planar 2-D mesh, and face 5 of a sector (the record of face 6
  !> computes the face).
  function Rotational_Periodic_Face_Metric ( Im, Jm, Km, Fm, Blk ) result ( M )
    use ARES_Lib_Metrics, only : delthe, cos_delthe, sin_delthe
    use ARES_Lib_Diffusive, only : Face_Metric
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(in) :: Blk
    real(R8) :: M(3,3)
    ! Local
    integer :: Kg, Ks, d
    real(R8) :: Mg(3,3), g2, g3

    M = 0d0
    ! The faces whose flux BC_Rotational_Periodic_Visc computes (the same two tests)
    if (ndir < 3 .and. (delthe >= 0d0 .and. delthe <= 0d0)) return
    if (ndir == 3 .and. Fm /= 6) return

    Kg = Km - guide(Fm,3)          ! k of the image of the cell (its ghost)
    Ks = Km                        ! k of the source of the image: the cell itself on a wedge,
    if (ndir == 3) Ks = Blk % dim(3) - Km + 1   ! the cell next to the opposite k face on a sector

    if (ndir == 3) then
      ! Metric of the image, the source turned: x' = R x gives grad xi' = R grad xi, each row turned as a vector
      Mg = Blk % M(Im,Jm,Ks) % c
      do d = 1, 3
        g2 = Mg(d,2) ; g3 = Mg(d,3)
        Mg(d,2) = g2 * cos_delthe + guide(Fm,3) * g3 * sin_delthe
        Mg(d,3) = g2 * sin_delthe * (-guide(Fm,3)) + g3 * cos_delthe
      end do
    else
      Mg = Blk % M(Im,Jm,Kg) % c
    end if
    M = Face_Metric ( Blk % M(Im,Jm,Km) % c, Mg )

  end function Rotational_Periodic_Face_Metric

end module ARES_Lib_BC_Fluxes_Rotational