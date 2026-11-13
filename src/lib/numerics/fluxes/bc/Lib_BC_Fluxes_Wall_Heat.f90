module ARES_Lib_BC_Fluxes_Wall_Heat
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_Heat ( Im, Jm, Km, Fm, Blk, Heat_Flux, Ovar, rec )
    use ARES_Lib_Fluid
    use ARES_Lib_RANS
    use FLINT_Lib_Thermodynamic
    use ARES_Config_Types_m, only: obj_rans
    use ARES_Lib_Prt_Correction, only: delta_Prt

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: Heat_Flux
    type(ARES_block_type), intent(inout), target :: Blk
    real(R8), optional, dimension(8), intent(inout) :: Ovar
    ! the record of the face: its wall stencil (diffusive-metric = stencil)
    type(ARES_bc_type), optional, intent(in), target :: rec
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k
    real(R8) :: Normal(3), Area, Dist, Prim(nprim), rho, temp, Gradient(nprim,3)
    real(R8), pointer, contiguous :: M(:,:), Mc(:,:)   ! metrics of the face and of the cell, read in place
    logical  :: wst
    real(R8) :: mil, kl, Stress(3), Flux(nprim), Prim_Wall(nprim), rho_wall, Twall, dl, entalpy, pressure
    real(R8) :: hs, mit, cp, kt ! roughness variables
    real(R8) :: dPrt ! Turbulent Prandtl correction for rough wall

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = 1d-20 ! since the flux is computed at the wall
    Mc => Blk % M (Im,Jm,Km) % c
    M => Mc
    ! diffusive-metric = stencil: velocity and turbulence take the metric of the wall face (Wall_Face_Stencils in
    ! ARES_Lib_Stencil_Metrics), with no difference of their wall values along the wall: the velocity of a no-slip wall
    ! is uniform along it, while the turbulence values that vary along it (omega of the k-omega models, from the wall
    ! distance of the cell; nu~ on a rough wall) leave out their derivative along the wall (docs/theory/numerics.md); the
    ! wall temperature and the wall values of the turbulence variables below are estimated with the metric of the cell
    wst = .false.
    if ( present(rec) ) wst = rec % wstencil
    if ( wst ) M => rec % wM

    ! Roughness
    if ( obj_rans%rough ) then 
      hs = Blk%ks(Im,Jm,Km)
    else  
      hs = 0d0
    end if

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    pressure = Prim(np)
    rho  = ph2vars( pressure, Prim(nh), rho_tab ) 
    temp = ph2vars( pressure, Prim(nh), T_tab   ) 

    ! Approximation: Twall==Tcell. Therefore, mil(wall) = mil(cell) and mit = 0 -> Strictly valid for null HeatFlux
    mil = ph2vars( pressure, Prim(nh), mi_tab2D )
    kl  = ph2vars( pressure, Prim(nh), k_tab2D  )
    cp  = ph2vars( pressure, Prim(nh), hT_tab   )
    
    ! Initialization
    Gradient = 0d0 ! concentration gradient is zero for one specie fluid
    Prim_Wall = Prim        ! -------------------------------- ! 

    ! Velocity gradient
    Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3

    ! RANS stuff
    if (model==2) then 
      dl = Blk%yn(Im,Jm,Km)
      call RANS_Set_Wall_Values( mil_wall=mil, rans_wall=Prim_Wall(nt:nprim), rans_cell=Prim(nt:nprim), &
                                 metric=dot_product( Mc(Dir,:), normal), dist=dl, ks=hs)
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho ) * modfm3
      call Eddy_Viscosity ( mut=mit, rans_variables=Prim_Wall(nt:nprim), mul=mil, rho=rho, vel_gradient=Gradient(nu:nw,:), &
                            walldist=0d0, ks=hs)   ! wall face: y = 0 (dl stays the cell distance for Set_Wall_Values)
    else
      mit = 0d0
    end if

    ! Temperature gradient
    if ( obj_rans%Prt_correction ) then
      ! Evaluated at the wall face (y=0). The rough-wall value mitilde_w = vk*u_tau*rho*0.03*hs
      ! yields the same u_tau as the cell value taken at y=dl, but leaves the near-wall damping
      ! exp(-y/hs) equal to 1, as required for a flux computed on the wall face itself.
      dPrt = delta_Prt ( mil, cp, kl, Prim_Wall(nt), hs, 0d0 )
    else
      dPrt = 0.00d0
    end if
    kt = mit*cp/(obj_rans%Prt+dPrt)
    ! Heat_Flux is positive into the fluid on every face: the wall is hotter than the cell when it heats the fluid
    Gradient(nh,Dir) = -Heat_Flux / ( 2d0*(kl+kt)*dot_product( Mc(Dir,:), normal ) )

    ! Wall variables
    Twall = temp - Gradient(nh,Dir)
    entalpy   = pT2h( pressure, Twall ) 
    rho_wall  = ph2vars( pressure, entalpy, rho_tab )  
    Prim_Wall(nh) = entalpy ! -------------------------------- ! 

    Gradient = matmul ( Gradient, M )

    ! Stress vector (turbulent variables at the wall face: Rij = 0 and k = 0)
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, mit, Prim_Wall(nt:) )

    ! Fluxes
    Flux(np) = 0d0 
    Flux(nu:nw) = Stress * Area
    Flux(nh) = Area * Heat_Flux * modfm2   ! k grad T . n: Heat_Flux on faces 2, 4, 6, -Heat_Flux on 1, 3, 5

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dist )
      ! The ghost cells of the wall, turbulence variables included, are those of Fill_Ghost_Cell (Ghost_Wall), from the
      ! same wall values: one set of ghost cells for the fluxes and the source terms
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    ! qw given as k grad T . n along the face normal of the solver, as the other walls (Heat_Flux is positive into the
    ! gas): the wall output writes Heat_Flux itself on every face
    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), into=-modfm2, &
                                                    qw=modfm2 * Heat_Flux,                    &
                                                    y=blk%dl(im,jm,km)%c(dir)*0.5d0,          &
                                                    Tw=Twall, rhow=rho_wall, mu=mil, exit_array=Ovar)

  end subroutine BC_Wall_Heat

end module ARES_Lib_BC_Fluxes_Wall_Heat