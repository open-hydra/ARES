!> The manifold (BC 501): the faces of a group joined by a lumped plenum. The gas that leaves the domain through the
!> outlet faces of the group (role 1) comes back through its inlet faces (role 2). Once per stage, before the boundary
!> fluxes, Manifold_Gather sums over the faces of each group on every rank and sets the state of the group; BC_Fluxes
!> applies it with BC_Inlet_MassFlux_H0 on the inlet faces and the outlet 406 on the outlet faces. The roles are those
!> of the records: no sign of the flow changes them. The groups are built at set-up, per grid level, from the records
!> every rank has read (Manifold_Setup); the faces of a group may lie on any ranks and be cut by the decomposition.
module ARES_Lib_BC_Fluxes_Manifold
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index, Compute_Modfm

  implicit none
  private
  public :: Manifold_Setup, Manifold_Gather, MANIFOLD_OUTLET, MANIFOLD_INLET

  integer, parameter :: MANIFOLD_OUTLET = 1   !< role of a face where the gas leaves the domain into the manifold
  integer, parameter :: MANIFOLD_INLET  = 2   !< role of a face where it comes back in

contains

  !─────────────────────────────────────────────────────────────────────────────
  !> The groups of the manifold records of a domain (one grid level), from the records every rank has read: one group
  !> per id, in the order the ids first appear; its parameters from its first record, which every other record of the
  !> group must repeat (ATLAS writes them group-wide); at least one record of each role in the whole domain. Each record
  !> gets the index of its group (man_idx). The run stops on every rank on a bad file.
  subroutine Manifold_Setup ( domain, level )
    use ARES_Mod_MPI, only: mpi_stop_all, mpi_is_root
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    integer,                intent(in)    :: level
    ! Local
    integer :: i, m, n501, nman
    integer, allocatable :: gids(:)
    character(len=400) :: msg

    domain % nman = 0
    if (allocated(domain % man)) deallocate(domain % man)
    n501 = 0
    do i = 1, domain % nbound
      if (domain % bc(i) % type == 501) n501 = n501 + 1
    end do
    allocate ( domain % man(0) )
    if (n501 == 0) return

    ! The distinct ids, in the order of the file (the values of each record were checked when read: Need_Manifold)
    allocate ( gids(n501) )
    nman = 0
    do i = 1, domain % nbound
      if (domain % bc(i) % type /= 501) cycle
      if ( .not. any ( gids(1:nman) == domain % bc(i) % man_gid ) ) then
        nman = nman + 1
        gids(nman) = domain % bc(i) % man_gid
      end if
    end do

    deallocate ( domain % man )
    allocate ( domain % man(nman) )
    domain % nman = nman
    do m = 1, nman
      domain % man(m) % gid = gids(m)
      allocate ( domain % man(m) % ci_in(max(nrans,1)) )
      domain % man(m) % ci_in = 0d0
    end do

    ! The records of each group: the parameters of the first, repeated by the others; the roles counted
    do i = 1, domain % nbound
      if (domain % bc(i) % type /= 501) cycle
      do m = 1, nman
        if (gids(m) == domain % bc(i) % man_gid) exit
      end do
      domain % bc(i) % man_idx = m
      associate ( g => domain % man(m), r => domain % bc(i) )
        if (g % n_out + g % n_in == 0) then
          g % K = r % man_K ; g % rf = r % rel_fac ; g % mratio = r % man_mratio ; g % p0ratio = r % man_p0ratio
          g % dh0 = r % man_dh0
          g % first = i
        else
          call Same ( 'K',          r % man_K,       g % K )
          call Same ( 'rel_fac',    r % rel_fac,     g % rf )
          call Same ( 'mdot_ratio', r % man_mratio,  g % mratio )
          call Same ( 'p0_ratio',   r % man_p0ratio, g % p0ratio )
          call Same ( 'dh0',        r % man_dh0,     g % dh0 )
        end if
        if (r % man_role == MANIFOLD_OUTLET) then
          g % n_out = g % n_out + 1
        else
          g % n_in = g % n_in + 1
        end if
      end associate
    end do

    do m = 1, nman
      if (domain % man(m) % n_out == 0 .or. domain % man(m) % n_in == 0) then
        write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] manifold group ', gids(m), ' of grid level ', level, ' has ', &
          domain % man(m) % n_out, ' outlet and ', domain % man(m) % n_in, ' inlet records: a manifold joins at least '// &
          'one face where the gas leaves the domain (role 1) and one where it comes back (role 2)'
        call mpi_stop_all ( trim(msg) )
      end if
    end do

    if (mpi_is_root) then
      write(*,'(A,I0,A)') ' Manifold groups of grid level ', level, ':'
      do m = 1, nman
        associate ( g => domain % man(m) )
          write(*,'(A,I0,A,I0,A,I0,A,ES10.3,A,F6.3,A,ES10.3,A,ES10.3,A,ES12.5,A)') '   group ', g % gid, ': ', g % n_out, &
            ' outlet and ', g % n_in, ' inlet records; K ', g % K, ', rel_fac ', g % rf, ', mdot-ratio ', g % mratio, &
            ', p0-ratio ', g % p0ratio, ', dh0 ', g % dh0, ' J/kg'
        end associate
      end do
    end if

  contains

    subroutine Same ( what, v, v1 )
      character(len=*), intent(in) :: what
      real(R8),         intent(in) :: v, v1
      if ( abs(v - v1) <= 1d-12 * max(1d0, abs(v1)) ) return
      write(msg,'(A,I0,A,I0,A,I0,3A,ES15.8,A,I0,A,ES15.8)') '[ERROR] manifold (501) record ', i, ' of grid level ', &
        level, ' (group ', domain % bc(i) % man_gid, '): its ', trim(what), ' is ', v, ' while record ', &
        domain % man(m) % first, ' of the group has ', v1
      call mpi_stop_all ( trim(msg)//'; the parameters of a manifold are those of its section, the same on every record' )
    end subroutine Same

  end subroutine Manifold_Setup


  !─────────────────────────────────────────────────────────────────────────────
  !> The state the faces of each manifold group take at this stage, from the cells of the group's faces on every rank.
  !> Called by every thread of the parallel region at the top of BC_Fluxes: one thread sums over the records of this rank,
  !> reduces the sums over the ranks (MPI from the primary thread, as the connection exchange) and sets the targets; the
  !> others wait at the barrier.
  !>
  !> Outlet faces, the cells whose gas leaves: the mass flow sum g dA (g = rho u.n), and with it the total enthalpy
  !> h0 = h + v^2/2, the dynamic pressure qd = rho v^2/2 and the specific turbulence variables, each as sum g phi dA; every
  !> cell gives its area. Inlet faces, every cell: the area, the static pressure sum p dA, and the dynamic pressure
  !> sum g qd dA with its own g = rho |u.n|. Then, for the group, the mass-flow-weighted means: the inlet faces take the
  !> mass flux g_in = mdot-ratio mdot / A_in, the total enthalpy h0_in = h0_out - dh0 and the mean turbulence; the outlet
  !> faces take the static pressure p_out = p0_out - q_out of the total pressure p0_out = p0_in / p0-ratio + K q_out, with
  !> p0_in = mean p of the inlet cells + their mass-weighted qd: the total pressure is the same on the two sides of a
  !> loss-free plenum (K = 0, p0-ratio = 1), and mass and total enthalpy cross it exactly (mdot-ratio = 1, dh0 = 0).
  subroutine Manifold_Gather ( domain )
    use ARES_Mod_MPI,        only: mpi_size_, mpi_allreduce_sum_r8_array, mpi_is_root
    use ARES_Config_Types_m, only: obj_io
    use FLINT_Lib_Thermodynamic
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer  :: nsum, m, ii, i, Bm, Im, Jm, Km, Fm, dir, Face_i, Face_j, Face_k, modfm, modfm1, modfm2, modfm3, k
    integer  :: iAout, imdot, igh0, igq, iAin, ipA, igin, igqin, iturb
    real(R8), allocatable :: sums(:,:), flat(:)
    real(R8) :: A, n(3), Prim(nprim), rho, Un, g, q, qd
    real(R8) :: A_in, q_out, q_in, p_in

    if (domain % nman == 0) return

    !$omp master
    ! The slots of the sums of a group
    iAout = 1           ! area of the outlet faces
    imdot = 2           ! mass flow leaving through the outlet cells, sum g dA
    igh0  = 3           ! sum g h0 dA over them
    igq   = 4           ! sum g qd dA over them, qd = rho v^2/2
    iAin  = 5           ! area of the inlet faces
    ipA   = 6           ! sum p dA over the inlet cells
    igin  = 7           ! sum g dA over the inlet cells, g = rho |u.n|
    igqin = 8           ! sum g qd dA over the inlet cells
    iturb = 8           ! then the nrans turbulence variables, sum g phi dA over the leaving outlet cells (phi specific)
    nsum  = iturb + nrans
    allocate ( sums(nsum, domain % nman) )
    sums = 0d0

    do ii = 1, domain % n_local_bc
      i = domain % local_bc_idx(ii)
      if (domain % bc(i) % type /= 501) cycle
      m  = domain % bc(i) % man_idx
      Bm = domain % bc(i) % b ; Im = domain % bc(i) % i ; Jm = domain % bc(i) % j ; Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      call Compute_Modfm ( Fm, modfm, modfm1, modfm2, modfm3 )
      call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )
      associate ( blk => domain % blk(Bm) )
        A    = blk % dir(dir) % f(Face_i,Face_j,Face_k) % a
        n    = blk % dir(dir) % f(Face_i,Face_j,Face_k) % n
        Prim = blk % P(:,Im,Jm,Km)
      end associate
      rho = ph2vars ( Prim(np), Prim(nh), rho_tab )
      Un  = real(modfm2, R8) * dot_product ( Prim(nu:nw), n )   ! normal velocity out of the domain
      q   = 0.5d0 * sum ( Prim(nu:nw)**2 )                      ! specific kinetic energy
      qd  = rho * q                                             ! dynamic pressure
      select case ( domain % bc(i) % man_role )
      case ( MANIFOLD_OUTLET )
        sums(iAout,m) = sums(iAout,m) + A
        if (Un > 0d0) then
          g = rho * Un
          sums(imdot,m) = sums(imdot,m) + g * A
          sums(igh0,m)  = sums(igh0,m)  + g * ( Prim(nh) + q ) * A
          sums(igq,m)   = sums(igq,m)   + g * qd * A
          do k = 1, nrans
            sums(iturb+k,m) = sums(iturb+k,m) + g * Prim(nt+k-1) / rho * A
          end do
        end if
      case ( MANIFOLD_INLET )
        g = rho * abs(Un)
        sums(iAin,m)  = sums(iAin,m)  + A
        sums(ipA,m)   = sums(ipA,m)   + Prim(np) * A
        sums(igin,m)  = sums(igin,m)  + g * A
        sums(igqin,m) = sums(igqin,m) + g * qd * A
      end select
    end do

    if (mpi_size_ > 1) then
      flat = reshape ( sums, [nsum * domain % nman] )
      call mpi_allreduce_sum_r8_array ( flat, nsum * domain % nman )
      sums = reshape ( flat, [nsum, domain % nman] )
    end if

    ! The state of each group
    do m = 1, domain % nman
      associate ( gr => domain % man(m) )
        A_in = sums(iAin,m)
        p_in = 0d0
        if (A_in > 0d0) p_in = sums(ipA,m) / A_in
        gr % mdot_out = sums(imdot,m)
        gr % active   = gr % mdot_out > 0d0 .and. A_in > 0d0
        if (gr % active) then
          ! The mass-flow-weighted means of the gas that leaves
          gr % h0_out = sums(igh0,m) / gr % mdot_out
          q_out       = sums(igq,m)  / gr % mdot_out
          ! The inlet faces: the mass flux, the total enthalpy carried over (less dh0) and the turbulence of the gas
          ! that leaves
          gr % g_in   = gr % mratio * gr % mdot_out / A_in
          gr % h0_in  = gr % h0_out - gr % dh0
          do k = 1, nrans
            gr % ci_in(k) = sums(iturb+k,m) / gr % mdot_out
          end do
          ! The outlet faces: the static pressure of the total pressure that the inlet side gives the group, through the
          ! loss and the ratio (the dynamic pressures of the two sides, mass-weighted)
          q_in = 0d0
          if (sums(igin,m) > 0d0) q_in = sums(igqin,m) / sums(igin,m)
          gr % p0_in  = p_in + q_in
          gr % p0_out = gr % p0_in / gr % p0ratio + gr % K * q_out
          gr % p_out  = gr % p0_out - q_out
          gr % mdot_in = gr % g_in * A_in
        else
          ! Nothing leaves yet: the inlet faces are slip walls, the outlet faces take the pressure of the inlet side
          gr % g_in  = 0d0
          gr % p_out = p_in
          gr % mdot_in = 0d0 ; gr % h0_out = 0d0 ; gr % h0_in = 0d0 ; gr % p0_out = p_in ; gr % p0_in = p_in
        end if
      end associate
    end do

    ! One line per group at the shell interval, once per iteration (the first stage that sees it)
    if ( mpi_is_root .and. obj_io % shell_diter > 0 ) then
      if ( mod ( domain % iter, obj_io % shell_diter ) == 0 ) then
        do m = 1, domain % nman
          associate ( gr => domain % man(m) )
            if ( gr % last_line == domain % iter ) cycle
            gr % last_line = domain % iter
            write(*,'(A,I0,A,I0,A,I0,A,2ES17.9,A,2ES17.9,A,2ES17.9,A,ES11.3,A)') ' manifold ', gr % gid, ', grid level ', &
              domain % mg_level, ', iteration ', domain % iter, ': mdot out/in', gr % mdot_out, gr % mdot_in, &
              ' kg/s, h0 out/in', gr % h0_out, gr % h0_in, ' J/kg, p0 out/in', gr % p0_out, gr % p0_in, ' Pa, taken out', &
              gr % mdot_out * gr % dh0, ' W'
          end associate
        end do
      end if
    end if

    deallocate ( sums )
    !$omp end master
    !$omp barrier

  end subroutine Manifold_Gather

end module ARES_Lib_BC_Fluxes_Manifold
