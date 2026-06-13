module ARES_Mod_Explicit
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Explicit_Step

contains

  subroutine Explicit_Step ( domain, External_Function )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m
    use ARES_Global_m
    use ARES_Mod_dt,         only: Set_Global_dt, Compute_dt
    use ARES_Lib_Ghost,      only: Fill_Ghost_Cell, Ghost_Wall_Extrapolation
    use ARES_Mod_Fluxes,     only: Fluxes
    use ARES_Mod_BC_Fluxes,  only: BC_Fluxes
    use ARES_Lib_RANS,       only: RANS_Source_Terms
    use ARES_Mod_Newstate,   only: RK_Newstate
    use ARES_Mod_Diagnostic, only: Compute_Residual
    !use ARES_Lib_RotatingFrame, only: RotatingFrame_Source_Terms
    use ARES_Mod_MPI, only: is_local_block, mpi_reduce_sum_r8_array, &
                           mpi_is_root, mpi_bcast_logical, mpi_bcast_integer
    implicit none
    type(ARES_domain_type), intent(inout) :: domain(obj_multigrid%MGL)
    external :: External_Function
    ! Local
    logical  :: endsim, iosim, endmg
    integer  :: i_rk, i_strang, b, level
    real(R8) :: average(nres)

    level = obj_multigrid%MG_level
    obj_multigrid%change_MG = .false.
    domain(level) % iter = domain(level) % iter + 1
    obj_sim_param%iter_from_call = obj_sim_param%iter_from_call + 1
    obj_sim_param%iter_general   = obj_sim_param%iter_general + 1

    if (obj_sim_param%HYDRA_time_accurate) then
      call Set_Global_dt ( domain(level) )
      domain(level) % time = domain(level) % time + domain(level) % dtglobal
    else
      domain(level) % dtglobal = 1d5
      call Compute_dt ( domain(level), obj_time_scheme%cfl, obj_time_scheme%vnn, obj_time_scheme%rampa_cfl_iter )  ! Compute local and minimum time step
      if ( obj_time_scheme%time_accurate ) then
        call Set_Global_dt ( domain(level) )  ! Time-accurate: apply global minimum time step
        domain(level) % time = domain(level) % time + domain(level) % dtglobal 
      endif
    endif

    !$omp parallel

    call Copy_State ( domain(level) )

    rk: do i_rk = 1, obj_time_scheme%n_RK

      call Fill_Ghost_Cell ( domain(level) )               ! Fill ghost cells
      call Fluxes ( domain(level) )                        ! Convective and Diffusive fluxes
      call BC_Fluxes ( domain(level) )                     ! Boundary fluxes

      call External_Function ( domain(level) )             ! External function (e.g. source terms)

      if (model==2) then
        call Ghost_Wall_Extrapolation ( domain(level) )    ! Fill ghost cells for RHS of turbulence eqns
        call RANS_Source_Terms ( domain(level) )           ! RHS of turbulence equations
      endif

      call RK_Newstate ( domain(level), i_rk )             ! State update

    enddo rk

    !$omp end parallel

    ! Save residuals
    if ( (obj_sim_param%iter_from_call == 1) .or. (mod (domain(level) % iter, obj_io%res_diter) == 0) ) then
      obj_sim_param%residuotot = 0.d0
      do b = 1, domain(level) % nb
        if (.not. is_local_block(b)) cycle
        call Compute_Residual ( new=domain(level)%blk(b)%P, &
                                old=domain(level)%blk(b)%PO, &
                                dt=domain(level)%blk(b)%dtlocal, &
                                n=domain(level)%blk(b)%dim, &
                                average=average, &
                                total=obj_sim_param%residuotot )
      enddo
      call mpi_reduce_sum_r8_array(obj_sim_param%residuotot, nres)
      if (mpi_is_root) obj_sim_param%residuotot = sqrt ( obj_sim_param%residuotot ) ! L2 norm time derivative
    endif

    ! Determine simulation control flags on root, then broadcast to all ranks
    if (mpi_is_root) then
      endsim = .false.
      endmg  = .false.
      if (level == 1) then
        endsim = ( obj_sim_param%iter_from_call >= domain(1) % itermax ) &
            .or. ( obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold ) &
            .or. ( domain(1) % time >= obj_sim_param%time_threshold )
      else
        endmg = ( obj_sim_param%iter_from_call >= domain(level) % itermax ) &
           .or. ( obj_sim_param%residuotot(1) <= obj_sim_param%res_threshold ) &
           .or. ( domain(level) % time >= obj_sim_param%time_threshold )
      endif

      iosim  = ( mod (domain(level) % iter, obj_io%sol_diter) == 0) &
          .or. ( domain(level) % time >= obj_sim_param%time_from_call + obj_io%sol_dtime )

      if ( endsim ) then
        obj_sim_param%TODO = 3
      elseif ( iosim .or. endmg ) then
        obj_sim_param%TODO = 2
        if (endmg .and. obj_multigrid%MGL > 1) then
          obj_multigrid%change_MG = .true.
        endif
      else
        obj_sim_param%TODO = 1
      endif
    end if

    ! Broadcast simulation control from root to all ranks
    call mpi_bcast_integer(obj_sim_param%TODO)
    call mpi_bcast_logical(obj_multigrid%change_MG)
      
  end subroutine Explicit_Step


  subroutine Copy_State ( domain )
    use ARES_Advanced_Types_m
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, j, k, b

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        
        domain % blk(b) % PO (:,i,j,k) = domain % blk(b) % P (:,i,j,k)
      
      enddo; enddo; enddo
      !$omp end do
    enddo

  end subroutine Copy_State

end module ARES_Mod_Explicit