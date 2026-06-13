module ARES_Lib_BC_Fluxes_Wall_Temperature
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm, Compute_Wall_Properties

  implicit none
  public

contains

  subroutine BC_Wall_Temperature ( Im, Jm, Km, Fm, Blk, T_wall, Ovar )
    use ARES_Lib_Fluid
    use ARES_Lib_RANS
    use FLINT_Lib_Thermodynamic
    use ARES_Config_Types_m, only: obj_rans
    use ARES_Lib_Prt_Correction, only: delta_Prt

    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    real(R8), intent(in) :: T_wall
    type(ARES_block_type), intent(inout) :: Blk
    real(R8), optional, dimension(8), intent(out) :: Ovar
    ! Local
    integer :: modfm, modfm1, modfm2, modfm3, Dir, Face_i, Face_j, Face_k, ig, jg, kg
    real(R8) :: Normal(3), Area, Dist, M(3,3), Prim(nprim), rho, T, rho_wall, dl
    real(R8) :: Prim_wall(nprim), mil, kl, Gradient(nprim,3), Stress(3), Flux(nprim), pressure, entalpy
    real(R8) :: hs, mit, cp, kt ! roughness variables
    real(R8) :: dPrt ! Turbulent Prandtl correction for rough wall

    call Compute_Modfm ( fm, modfm, modfm1, modfm2, modfm3 )
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a
    Dist = 1d-20 ! since the flux is computed at the wall
    M = Blk % M (Im,Jm,Km) % c

    ! Roughness
    if ( obj_rans%rough ) then 
      hs = Blk%ks(Im,Jm,Km)
    else  
      hs = 0d0
    end if

    ! boundary cell variables
    Prim = Blk % P(:,Im,Jm,Km)
    rho = ph2vars( Prim(np), Prim(nh), rho_tab )  ! cell density
    T   = ph2vars( Prim(np), Prim(nh), T_tab   ) ! cell Temperature

    ! Initialization
    Gradient = 0d0

    ! Boundary layer: p(wall) = p(cell). Also, non catalytic wall, so R(wall) = R(cell) and ci(wall) = ci(cell), but not roi
    pressure  = Prim(np) ! wall pressure = cell pressure
    entalpy   = pT2h(pressure, T_wall ) ! Wall entalpy
    rho_wall  = ph2vars(pressure, entalpy, rho_tab )   ! Wall density

    Prim_Wall = Prim        ! -------------------------------- ! 
    Prim_Wall(nh) = entalpy ! -------------------------------- !     

    ! Laminar viscosity/conductivity at interface computed at Twall
    mil = ph2vars( pressure, entalpy, mi_tab2D )
    kl  = ph2vars( pressure, entalpy, k_tab2D  )

    ! Difference in Across-face direction: symmetry; tangential direction: null
    Gradient(nu:nw,Dir) = Prim(nu:nw) * modfm3                  ! Velocity
    Gradient(nh,Dir) = ( T - T_wall ) * modfm3                  ! Temperature

    if (model==2)  then
      dl = blk % yn(im,jm,km)
      call RANS_Set_Wall_Values( mil_wall=mil, rans_wall=Prim_Wall(nt:nprim), rans_cell=Prim(nt:nprim), &
                                 metric=dot_product( M(Dir,:), normal ), dist=dl, ks=hs )
      Gradient(nt:nprim,Dir) = ( Prim(nt:nprim)/rho - Prim_Wall(nt:nprim)/rho_wall ) * modfm3 
      call Eddy_Viscosity ( mut=mit, rans_variables=Prim_Wall(nt:nprim), &
                            mul=mil, rho=rho_wall, vel_gradient=Gradient(nu:nw,:), &
                            walldist=dl, ks=hs)
    else 
      mit = 0d0
    end if

    ! Transformation
    Gradient = matmul ( Gradient, M )

    ! Stress vector
    Stress = Stress_Vector ( Gradient(nu:nw,:), Normal, mil, mit, Prim(nt:) )
    cp = ph2vars(pressure, entalpy, hT_tab ) 
    if ( obj_rans%Prt_correction ) then
      dPrt = delta_Prt ( mil, cp, kl, Prim(nt), hs, dl )
    else
      dPrt = 0.00d0
    end if
    kt = mit*cp/(obj_rans%Prt+dPrt)

    ! Fluxes
    Flux(1:np) = 0d0 
    Flux(nu:nw) = Stress * Area
    Flux(nh) = Area * (kl+kt) * dot_product ( Gradient(nh,:), Normal )

    if (model==2) then
      call RANS_Diffusive_Flux ( flux=Flux(nt:nprim), &
                                 rans_variables=Prim_Wall(nt:nprim), &
                                 vel_gradient=Gradient(nu:nw,:), &
                                 rans_gradient=Gradient(nt:nprim,:), &
                                 mul=mil, rho=rho_wall, &
                                 area=area, normal=normal, dist=dist )
      ! Ghost-cell extrapolation of RANS variables
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)
      call RANS_Extrapolate_Wall ( Prim(nt:nprim), Prim_Wall(nt:nprim), &
                                   rho, rho_wall, Blk % P(nt:nprim,Ig,Jg,Kg) )
    endif

    ! Residual update
    Blk % r (:,Im,Jm,Km) = Blk % r (:,Im,Jm,Km) - modfm2 * Flux

    if (present(Ovar)) call Compute_Wall_Properties(stress=Stress, pw=prim(np), Tw=T_wall,          &
                                                    rhow=rho_wall, mu=mil,                          &
                                                    y=blk%dl(im,jm,km)%c(dir)*0.5d0,                &
                                                    qw=(kl+kt) * dot_product ( Gradient(nh,:), Normal ), &
                                                    exit_array=Ovar)

  end subroutine BC_Wall_Temperature
  
end module ARES_Lib_BC_Fluxes_Wall_Temperature