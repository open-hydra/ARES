module ARES_Lib_Newstate
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Newstate_Preconditioned, Newstate_Primitive

contains


  subroutine Newstate_Primitive ( domain, irk )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_time_scheme, obj_irs
    use ARES_Global_m
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_RK
    use ARES_Lib_RANS
    use ARES_Lib_IRS
    use ARES_Mod_MPI, only: is_local_block

    implicit none
    type(ARES_domain_type), intent(inout) :: domain 
    integer, intent(in)                   :: irk
    ! Local
    integer :: i, j, k, s, b
    integer :: n_rk
    logical :: irs_enabled
    real(R8) :: irs_beta

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle

      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        ! Primitive variables residuals computation from conservation form residuals. 
        ! ResCons => ResPrim => ResPrim * dt/V. If no IRS, primitive variables
        ! are also updated. 
        call compute ( domain % blk(b) % P(:,i,j,k), domain % blk(b) % PO(:,i,j,k),    &
                       domain % blk(b) % r(:,i,j,k), domain % blk(b) % dtlocal(i,j,k), &
                       domain % blk(b) % vol(i,j,k), irk, irs_enabled )

        if ( .not.irs_enabled ) call check_and_fix_state(domain%blk(b)%P(:,i,j,k), b, i, j, k)

      enddo; enddo; enddo

    enddo

    if (irs_enabled) then

      call Residual_Smoothing ( domain, irs_beta )
      
      do b = 1, domain % nb ! Loop over blocks
        if (.not. is_local_block(b)) cycle

        !$omp do collapse(3)
        do k = 1, domain % blk(b) % dim(3)
        do j = 1, domain % blk(b) % dim(2)
        do i = 1, domain % blk(b) % dim(1)  

          !% Apply irk-RK stage to get newstate in primitive variables form.
          domain % blk(b) % P(:,i,j,k) = RK_Stage ( irk, n_rk, &
                                                        domain % blk(b) % P (:,i,j,k), &
                                                        domain % blk(b) % PO(:,i,j,k), &
                                                        domain % blk(b) % r (:,i,j,k)  )
          
          call check_and_fix_state(domain%blk(b)%P(:,i,j,k), b, i, j, k)

        enddo; enddo; enddo

      enddo ! end of loop over blocks

    endif

    contains

      subroutine compute ( prim, primold, residual, dt, volume, irk, irs_enabled )

        implicit none
        integer, intent(in) :: irk
        logical, intent(in) :: irs_enabled
        real(kind=8), intent(in) :: primold(nprim), dt, volume
        real(kind=8), intent(inout) :: prim(nprim), residual(nprim)
        ! Local
        real(kind=8) :: Tdiff, inv_volume, inv_rho, rho, velocity(3)
        real(kind=8) :: drdp, p, h, drh, Denom, Res(nprim)

        rho = ph2vars( prim(np), prim(nh), rho_tab ) 
        inv_volume  = 1d0 / volume
        inv_rho     = 1d0 / rho
        velocity = prim (nu:nw)

        !% Primitive residuals computation.
        p = prim(np)
        h = prim(nh)
        drdp = ph2vars( p, h, rp_tab ) ! Non preconditioned drho/dp
        drh  = ph2vars( p, h, rh_tab ) !  drho/dh
        Denom = 1.0_R8 / ( rho*drdp + drh )

        Res(np)  = Denom * ( ( (h - 0.5_R8*Norm2(velocity)**2)*drh + rho )*residual(np) + sum(velocity*drh*residual(nu:nw)) - drh*residual(nh) )
        Res(nu)   = -velocity(1)*inv_rho*residual(np) + inv_rho*residual(nu)
        Res(nv)   = -velocity(2)*inv_rho*residual(np) + inv_rho*residual(nv)
        Res(nw)   = -velocity(3)*inv_rho*residual(np) + inv_rho*residual(nw)
        Res(nh)   = Denom * ( ( (-h + 0.5_R8*Norm2(velocity)**2)*drdp + 1.0_R8 ) *residual(np) - sum(velocity*drdp*residual(nu:nw)) + drdp*residual(nh) )

        ! Turbulence (to do: include prod*vol into residual in spalart ...done)
        if (nprim > nh) Res (nt:nprim) = residual(nt:nprim)

        !% Time residual computation 
        residual = - Res * dt * inv_volume
        if ( .not.irs_enabled ) prim = RK_Stage ( irk, n_rk, prim, primold, residual )

      end subroutine compute 

  end subroutine Newstate_Primitive

  subroutine Newstate_Preconditioned ( domain, irk )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_time_scheme, obj_irs
    use ARES_Global_m
    use FLINT_Lib_Thermodynamic
    use ARES_Lib_RK
    use ARES_Lib_RANS
    use ARES_Lib_IRS
    use ARES_Lib_Preconditioning, only : Comp_Ur
    use ARES_Mod_MPI, only: is_local_block

    implicit none
    type(ARES_domain_type), intent(inout) :: domain 
    integer, intent(in)                   :: irk
    ! Local
    integer :: i, j, k, s, b
    integer :: n_rk
    logical :: irs_enabled
    real(R8) :: irs_beta

    n_rk        = obj_time_scheme%n_RK
    irs_enabled = obj_irs%enabled
    irs_beta    = obj_irs%beta

    do b = 1, domain % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle

      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        call compute ( domain % blk(b) % P(:,i,j,k), domain % blk(b) % PO(:,i,j,k),    &
                       domain % blk(b) % r(:,i,j,k), domain % blk(b) % dtlocal(i,j,k), &
                       domain % blk(b) % vol(i,j,k), irk, irs_enabled )

        if ( .not.irs_enabled ) call check_and_fix_state(domain%blk(b)%P(:,i,j,k), b, i, j, k)

      enddo; enddo; enddo

    enddo

    if (irs_enabled) then

      call Residual_Smoothing ( domain, irs_beta )
      
      do b = 1, domain % nb ! Loop over blocks
        if (.not. is_local_block(b)) cycle

        !$omp do collapse(3)
        do k = 1, domain % blk(b) % dim(3)
        do j = 1, domain % blk(b) % dim(2)
        do i = 1, domain % blk(b) % dim(1)  

          !% Apply irk-RK stage to get newstate in primitive variables form.
          domain % blk(b) % P(:,i,j,k) = RK_Stage ( irk, n_rk, &
                                                        domain % blk(b) % P (:,i,j,k), &
                                                        domain % blk(b) % PO(:,i,j,k), &
                                                        domain % blk(b) % r (:,i,j,k)  )
          
          call check_and_fix_state(domain%blk(b)%P(:,i,j,k), b, i, j, k)

        enddo; enddo; enddo

      enddo ! end of loop over blocks

    endif

    contains

      subroutine compute ( prim, primold, residual, dt, volume, irk, irs_enabled )

        implicit none
        integer, intent(in) :: irk
        logical, intent(in) :: irs_enabled
        real(kind=8), intent(in) :: primold(nprim), dt, volume
        real(kind=8), intent(inout) :: prim(nprim), residual(nprim)
        ! Local
        real(kind=8) :: Tdiff, inv_volume, inv_rho, rho, velocity(3)
        real(kind=8) :: p, h, drh, Denom, Res(nprim), cp, rho_T, sound, Ur, Theta

        p = prim(np)
        h = prim(nh)
        rho   = ph2vars( p, h, rho_tab    ) ! density
        sound = ph2vars( p, h, sound_tab  ) ! speed of sound
        drh   = ph2vars( p, h, rh_tab     ) ! drho/dh
        rho_T = ph2vars( p, h, dT_tab     ) ! drho/dT  
        cp    = ph2vars( p, h, hT_tab     ) ! dh/dT
        
        inv_volume  = 1d0 / volume
        inv_rho     = 1d0 / rho
        velocity = prim (nu:nw)
        
        Ur = Comp_Ur( Norm2(prim(nu:nw)), sound ) ! reference velocity
        Theta = 1.0_R8 / Ur**2 - rho_T / ( rho * cp ) ! Preconditioned drho/dp
        Denom = 1.0_R8 / ( rho*Theta + drh )
        
        Res(np)  = Denom * ( ( (h - 0.5_R8*Norm2(velocity)**2)*drh + rho )*residual(np) + sum(velocity*drh*residual(nu:nw)) - drh*residual(nh) )
        Res(nu)   = -velocity(1)*inv_rho*residual(np) + inv_rho*residual(nu)
        Res(nv)   = -velocity(2)*inv_rho*residual(np) + inv_rho*residual(nv)
        Res(nw)   = -velocity(3)*inv_rho*residual(np) + inv_rho*residual(nw)
        Res(nh)   = Denom * ( ( (-h + 0.5_R8*Norm2(velocity)**2)*Theta + 1.0_R8 ) *residual(np) - sum(velocity*Theta*residual(nu:nw)) + Theta*residual(nh) )

        ! Turbulence (to do: include prod*vol into residual in spalart ...done)
        if (nprim > nh) Res (nt:nprim) = residual(nt:nprim)

        !% Time residual computation 
        residual = - Res * dt * inv_volume
        if ( .not.irs_enabled ) prim = RK_Stage ( irk, n_rk, prim, primold, residual )

      end subroutine compute 
  end subroutine Newstate_Preconditioned

  subroutine check_and_fix_state ( prim, b, i, j, k )
    use ARES_Global_m
    use ARES_Lib_RANS
    use ARES_Mod_MPI, only: mpi_abort_all

    implicit none
    real(kind=8), dimension(:), intent(inout) :: prim
    integer, intent(in) :: b, i, j, k
    ! local


    ! Check integration
    if (isnan(product(prim(:))) .or. prim(np)<0d0 ) then
      write(*,*) "Integration failed in"
      write(*,*) b, i, j, k
      write(*,*) prim(1:nprim)
      call mpi_abort_all("NaN or p<0 detected")
    endif

    ! Check turbulence variables
    if (model==2) call RANS_Enforce_Realizability ( prim(nt:) )

  end subroutine check_and_fix_state

end module ARES_Lib_Newstate
