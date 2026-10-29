module ARES_Mod_BC_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  
  implicit none
  private
  public :: BC_Fluxes

contains

  subroutine BC_Fluxes ( domain )
    !---------------------------------------------------------------------------
    ! Thin wrapper: delegates to BC_Fluxes_core with explicit-shape arrays.
    ! This avoids Intel ifx 2023 -check bounds descriptor corruption when
    ! accessing allocatable components of a derived type simultaneously.
    !---------------------------------------------------------------------------
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(inout), target :: domain

    call BC_Fluxes_core ( domain%bc, domain%blk, domain%local_bf_first, domain%local_bf_count, &
                          domain%local_bf_idx, domain%nb, size(domain%bc), size(domain%blk), &
                          size(domain%local_bf_idx) )

  end subroutine BC_Fluxes


  subroutine BC_Fluxes_core ( bc, blk, bf_first, bf_count, bf_idx, nb, nbc, nblk, nidx )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_rans, obj_sim_param, obj_space_scheme, obj_riemann
    use ARES_Global_m, only: model, np, nrans
    use ARES_Mod_MPI, only: is_local_block
    use ARES_Lib_BC_Fluxes_Connection
    use ARES_Lib_BC_Fluxes_Rotational
    use ARES_Lib_BC_Fluxes_Symmetry
    use ARES_Lib_BC_Fluxes_Inflow
    use ARES_Lib_BC_Fluxes_Outflow
    use ARES_Lib_BC_Fluxes_Extrapolation
    use ARES_Lib_BC_Fluxes_Wall_Heat
    use ARES_Lib_BC_Fluxes_Wall_Temperature
    implicit none
    integer,               intent(in)    :: nb, nbc, nblk, nidx
    type(ARES_bc_type),    intent(inout) :: bc(nbc)
    type(ARES_block_type), intent(inout) :: blk(nblk)
    integer,               intent(in)    :: bf_first(nb,6), bf_count(nb,6), bf_idx(nidx)
    ! Local
    integer  :: f, lower, upper, i, ii, b
    integer  :: error
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Fs
    real(R8) :: T0, g, pstat
    real(R8) :: r_fc(3)
    real(R8) :: Sc, Sct, Prt
    logical  :: Prt_corr

    Sc  = obj_sim_param%Sc
    Sct = obj_rans%Sct
    Prt = obj_rans%Prt
    Prt_corr = obj_rans%Prt_correction

    ! BC fluxes are computed block by block and, within a block, face by face, over the
    ! records of the blocks this rank owns only (lists built once by build_local_bc_index):
    ! no rank scans the records of the whole domain. The records of one face of one block
    ! write distinct target cells (the periodic BC also writes the partner cell on the
    ! opposite face, one per record), so one worksharing loop over them is race-free; its
    ! implicit barrier keeps the face contributions of a corner cell summed in the order 1..6.
    blocks: do b = 1, nb
      if (.not. is_local_block(b)) cycle blocks
      faces: do f = 1, 6

        lower = bf_first(b,f)
        upper = lower + bf_count(b,f) - 1

        !$omp do schedule ( dynamic ) private(ii, i, Bm, Im, Jm, Km, Fm, error, Bs, Fs, T0, g, pstat, r_fc)
        do ii = lower, upper
          i  = bf_idx(ii)
          Bm = bc(i) % b
          Im = bc(i) % i
          Jm = bc(i) % j
          Km = bc(i) % k
          Fm = bc(i) % f
          select case ( bc(i) % type )

            case (101,102,201) ! connection & chimera (101=block connect, 102=chimera)
              call BC_Connection_Eul ( Im, Jm, Km, Fm, blk(Bm) )
              if (model>0) &
                call BC_Connection_Visc ( Im, Jm, Km, Fm, blk(Bm), bc(i) % Mg(1), bc(i) % Pg, &
                                          Sc, Sct, Prt, Prt_corr )
            case (103,104) ! multi-Solver coupling (104 = with chimera)
              call BC_Symmetry_Eul ( Im, Jm, Km, Fm, blk(Bm) )
              blk(Bm) % R(:,Im,Jm,Km) = blk(Bm) % R(:,Im,Jm,Km) + bc(i) % ext_flux

            case (200) ! periodic
              call BC_Rotational_Periodic_Eul ( Im, Jm, Km, Fm, blk(Bm) )

            case (300) ! symmetry
              call BC_Symmetry_Eul ( Im, Jm, Km, Fm, blk(Bm) )
              if (model>0) &
                call BC_Symmetry_Visc ( Im, Jm, Km, Fm, blk(Bm) )

            case (301) ! wall: prescribed heat flux
              call BC_Symmetry_Eul ( Im, Jm, Km, Fm, blk(Bm) )
              if (model>0) then
                call BC_Wall_Heat ( Im, Jm, Km, Fm, blk(Bm), bc(i) % qw )
              endif

            case (302) ! wall: prescribed temperature
              call BC_Symmetry_Eul ( Im, Jm, Km, Fm, blk(Bm) )
              if (model>0) then
                call BC_Wall_Temperature ( Im, Jm, Km, Fm, blk(Bm), bc(i) % Tw )
              endif

            case (400) ! extrapolation
              call BC_Extrapolation ( Im, Jm, Km, Fm, blk(Bm) )

            case (404) ! inlet: static temperature T + prescribed mass flux
              call BC_Inlet_MassFlux_T ( Bm, Im, Jm, Km, Fm, blk(Bm), &
                                          bc(i) % T0, bc(i) % mdot, &
                                          bc(i) % rel_fac, &
                                          bc(i) % alpha, bc(i) % beta, &
                                          bc(i) % ci(1:nrans), &
                                          error )
              if (error == 1) call BC_Symmetry_Eul ( Im, Jm, Km, Fm, blk(Bm) )

            case (405) ! inlet: supersonic — Mach + static temperature and static pressure
              call BC_Inlet_Supersonic_Static ( Bm, Im, Jm, Km, Fm, blk(Bm), &
                                                bc(i) % mach, &
                                                bc(i) % T0, bc(i) % pamb, &
                                                bc(i) % rel_fac, &
                                                bc(i) % alpha, bc(i) % beta, &
                                                bc(i) % ci(1:nrans), &
                                                bc(i) % mdot, error )
              if (error == 1) call BC_Symmetry_Eul ( Im, Jm, Km, Fm, blk(Bm) )

            case (406) ! outlet
              call BC_Outflow ( Bm, Im, Jm, Km, Fm, blk(Bm),bc(i) % pAmb, bc(i) % rel_fac, &
                                       bc(i) % mdot, error )
              if (error == 1) then
                call BC_Symmetry_Eul ( Im, Jm, Km, Fm, blk(Bm) )
              endif

          end select

        enddo

      enddo faces
    enddo blocks

  end subroutine BC_Fluxes_core

end module ARES_Mod_BC_Fluxes