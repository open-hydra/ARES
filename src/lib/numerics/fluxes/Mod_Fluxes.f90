module ARES_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Fluxes, Setup_Quasi2D_W, Zero_W_Residual, Zero_V_Residual

contains

    subroutine Fluxes ( domain )
      use ARES_Advanced_Types_m
      use ARES_Mod_MPI, only: is_local_block
      use ARES_Config_Types_m, only: obj_riemann, obj_rans, obj_sim_param

      implicit none
      type(ARES_domain_type), intent(inout) :: domain
      ! Local
      integer :: b
      real(R8) :: Sc, Sct, Prt
      logical :: Prt_corr

      Sc  = obj_sim_param%Sc
      Sct = obj_rans%Sct
      Prt = obj_rans%Prt
      Prt_corr = obj_rans%Prt_correction

      do b = 1, domain % nb ! Loop over blocks
        if (.not. is_local_block(b)) cycle

        call Fluxes_blk ( domain % blk(b) % P,   &
                          domain % blk(b) % r,   &
                          domain % blk(b) % dir, &
                          domain % blk(b) % dl,  &
                          domain % blk(b) % yn,  &
                          domain % blk(b) % ks,  &
                          domain % blk(b) % m,   &
                          domain % blk(b) % dim, &
                          Sc, Sct, Prt, Prt_corr )
      end do

    end subroutine Fluxes


    !> Planar quasi-2D case (one cell in k, ndir = 2, no wedge angle) whose k faces are all symmetry planes
    !> (300): the two mirror planes one cell apart give w the mode dw/dt = -2 a w / dz (a the speed of
    !> sound), far outside the stability region of the explicit schemes at the time step of the in-plane
    !> cells, so round-off in w grows until the run fails. As in planar codes, which do not solve the third
    !> direction, the z velocity is then set to zero and the z momentum is not integrated: w = 0 at every
    !> step, also when the initial field has w /= 0 (a field restarted from a run hit by the mode, or a 3-D
    !> field), which a frozen rho w would carry as an error on p and u. Null (0) or periodic (200) k faces,
    !> and axisymmetric wedges (delthe /= 0, where the cell-centre w of a cell off the wedge mid-plane
    !> carries part of the radial velocity), are left as they are. The k faces must also be planes normal to z
    !> (the x and y components of the unit normal of every k face at most tol_plane): the planar flow has
    !> w = 0. A layer whose thickness varies in the plane (k faces not normal to z) gives the gas a z velocity
    !> next to an inclined face, and freezing it would change the flow (the mass flow, the pressure): w is then
    !> integrated, with a warning, and its mode across the layer can limit the time step.
    !> A 1-D case (ndir = 1, one cell in j and in k) has the same mirror planes in both directions:
    !> symmetry on all the k faces gives w the mode above, and symmetry on all the j faces gives v the
    !> mode dv/dt = -2 a v / dy. The velocity across the planes, w and v, is then set to zero and not
    !> integrated either (freeze_w, freeze_v), with the same condition on the faces: the j faces planes normal
    !> to y for v, the k faces planes normal to z for w. A channel whose section varies along x, one cell
    !> across with the faces along its walls, keeps the velocity across, which is integrated, with a warning.
    !> Called after the initial field is in place and before the ghost cells and the coarse multigrid levels
    !> are filled from it.
    subroutine Setup_Quasi2D_W ( domain )
      use ARES_Advanced_Types_m
      use ARES_Global_m,    only: ndir, freeze_w, freeze_v, nv, nw
      use ARES_Lib_Metrics, only: delthe
      use ARES_Mod_MPI,     only: mpi_is_root, is_local_block, mpi_allreduce_max_r8
      implicit none
      type(ARES_domain_type), intent(inout) :: domain
      ! Local
      real(R8), parameter :: tol_plane = 1.0e-10_R8   ! Largest off-axis component of the unit normal of a face of a plane
      integer  :: i, j, k, b, nj, nk, nsymj, nsym, n(3)
      real(R8) :: offj, offk, offj_all, offk_all

      nj = 0
      nk = 0
      nsymj = 0
      nsym = 0
      do i = 1, domain % nbound
        select case ( domain % bc(i) % f )
        case (3, 4)
          nj = nj + 1
          if ( domain % bc(i) % type == 300 ) nsymj = nsymj + 1
        case (5, 6)
          nk = nk + 1
          if ( domain % bc(i) % type == 300 ) nsym = nsym + 1
        end select
      enddo
      ! No wedge: delthe exactly zero, written as two ordered comparisons (same result, also for -0 and NaN)
      freeze_w = ( ( ndir == 2 .or. ndir == 1 ) .and. delthe >= 0d0 .and. delthe <= 0d0 .and. nk > 0 .and. nsym == nk )
      freeze_v = ( ndir == 1 .and. nj > 0 .and. nsymj == nj )

      ! The symmetry faces must be planes normal to their axis: the j faces to y (freeze_v), the k faces to z (freeze_w);
      ! the largest values over the blocks of every rank, so that every rank takes the same decision
      if ( freeze_v .or. freeze_w ) then
        offj = 0.0_R8
        offk = 0.0_R8
        do b = 1, domain % nb
          if (.not. is_local_block(b)) cycle
          n = domain % blk(b) % dim
          do k = 1, n(3); do j = 0, n(2); do i = 1, n(1)
            offj = max ( offj, abs(domain % blk(b) % dir(2) % f(i,j,k) % N(1)), abs(domain % blk(b) % dir(2) % f(i,j,k) % N(3)) )
          enddo; enddo; enddo
          do k = 0, n(3); do j = 1, n(2); do i = 1, n(1)
            offk = max ( offk, abs(domain % blk(b) % dir(3) % f(i,j,k) % N(1)), abs(domain % blk(b) % dir(3) % f(i,j,k) % N(2)) )
          enddo; enddo; enddo
        enddo
        call mpi_allreduce_max_r8 ( offj, offj_all )
        call mpi_allreduce_max_r8 ( offk, offk_all )
        if ( freeze_v .and. .not. offj_all <= tol_plane ) then
          if ( mpi_is_root ) write(*,'(A,ES9.2,A)') &
            ' [WARNING] 1-D case with symmetry (300) j faces that are not planes normal to y (unit normal off y by up to', &
            offj_all, '): the y velocity is integrated, and its mode across the cells can limit the time step.'
          freeze_v = .false.
        endif
        if ( freeze_w .and. .not. offk_all <= tol_plane ) then
          if ( mpi_is_root .and. ndir == 2 ) write(*,'(A,ES9.2,A)') &
            ' [WARNING] Quasi-2D case with symmetry (300) k faces that are not planes normal to z (unit normal off z by up to', &
            offk_all, '): the z velocity is integrated, and its mode across the layer can limit the time step.'
          if ( mpi_is_root .and. ndir == 1 ) write(*,'(A,ES9.2,A)') &
            ' [WARNING] 1-D case with symmetry (300) k faces that are not planes normal to z (unit normal off z by up to', &
            offk_all, '): the z velocity is integrated, and its mode across the cells can limit the time step.'
          freeze_w = .false.
        endif
      endif

      if ( .not. ( freeze_w .or. freeze_v ) ) return

      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle
        if ( freeze_w ) domain % blk(b) % P(nw,:,:,:) = 0.0_R8
        if ( freeze_v ) domain % blk(b) % P(nv,:,:,:) = 0.0_R8
      enddo

      if ( mpi_is_root ) then
        if ( freeze_w .and. ndir == 2 ) write(*,'(A)') &
          ' [WARNING] Quasi-2D case with symmetry (300) k faces: the z velocity is set to zero and not integrated.'
        if ( freeze_w .and. ndir == 1 ) write(*,'(A)') &
          ' [WARNING] 1-D case with symmetry (300) k faces: the z velocity is set to zero and not integrated.'
        if ( freeze_v ) write(*,'(A)') &
          ' [WARNING] 1-D case with symmetry (300) j faces: the y velocity is set to zero and not integrated.'
      endif

    end subroutine Setup_Quasi2D_W


    !> With freeze_w, zero the z-momentum residual before the stage update (see Setup_Quasi2D_W).
    subroutine Zero_W_Residual ( domain )
      use ARES_Advanced_Types_m
      use ARES_Global_m, only: nw
      use ARES_Mod_MPI,  only: is_local_block
      implicit none
      type(ARES_domain_type), intent(inout) :: domain
      ! Local
      integer :: b, j, k, n(3)

      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle

        n = domain % blk(b) % dim

        !$omp do collapse (2)
        do k = 1, n(3)
        do j = 1, n(2)
          domain % blk(b) % R(nw,1:n(1),j,k) = 0.0_R8
        enddo; enddo

      enddo

    end subroutine Zero_W_Residual


    !> With freeze_v, zero the y-momentum residual before the stage update (see Setup_Quasi2D_W).
    subroutine Zero_V_Residual ( domain )
      use ARES_Advanced_Types_m
      use ARES_Global_m, only: nv
      use ARES_Mod_MPI,  only: is_local_block
      implicit none
      type(ARES_domain_type), intent(inout) :: domain
      ! Local
      integer :: b, j, k, n(3)

      do b = 1, domain % nb
        if (.not. is_local_block(b)) cycle

        n = domain % blk(b) % dim

        !$omp do collapse (2)
        do k = 1, n(3)
        do j = 1, n(2)
          domain % blk(b) % R(nv,1:n(1),j,k) = 0.0_R8
        enddo; enddo

      enddo

    end subroutine Zero_V_Residual


    subroutine Fluxes_blk ( Prim, Res, Dir, dl, yn, ks, M, n, Sc, Sct, Prt, Prt_corr)
      use ARES_Global_m
      use ARES_Lib_Convective
      use ARES_Lib_Diffusive
      use ARES_Base_Types_m

      implicit none
      integer, dimension(3), intent(in) :: n
      real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
      real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(out) :: Res
      type(d_metrics_type), dimension(3), intent(in) :: Dir
      type(ARES_vector_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: dl
      real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: yn
      real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: ks
      type(ARES_tensor_3D_type), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: M
      real(R8), intent(in) :: Sc, Sct, Prt
      logical, intent(in)  :: Prt_corr
      ! Local
      integer :: i, j, k


      ! Reset residuals to zero
      !$omp do
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1)
        Res(:,i,j,k) = 0d0
      end do ; enddo ; enddo

      !$omp do collapse (2)
      do k = 1, n(3)
      do j = 1, n(2)
      do i = 1, n(1) - 1
      
        call Convective_Flux ( dl(i-1:i+2,j,k) % c(1), &
                               Dir(1) % f(i,j,k) % n,  &
                               Dir(1) % f(i,j,k) % a,  &
                               Prim(:,i-1:i+2,j,k),    &
                               Res (:,i:i+1,j,k)       )

        if (model>0)  then
          call Diffusive_Flux ( Dir(1) % f(i,j,k) % n,  &
                                Dir(1) % f(i,j,k) % a,  &
                                yn(i  ,j,k),            &
                                yn(i+1,j,k),            &
                                ks(i  ,j,k),            &
                                ks(i+1,j,k),            &
                                Prim(:,i  ,j,k),        &
                                Prim(:,i+1,j,k),        &
                                Prim(:,i  ,j-1,k),      &
                                Prim(:,i  ,j+1,k),      &
                                Prim(:,i+1,j-1,k),      &
                                Prim(:,i+1,j+1,k),      &
                                Prim(:,i  ,j,k-1),      &
                                Prim(:,i  ,j,k+1),      &
                                Prim(:,i+1,j,k-1),      &
                                Prim(:,i+1,j,k+1),      &
                                M(i  ,j,k) % c,         &
                                M(i+1,j,k) % c,         &
                                Res (:,i  ,j,k),        &
                                Res (:,i+1,j,k),        &
                                1, 2, 3, Sc, Sct, Prt, Prt_corr )
        end if
                        

      enddo; enddo; enddo

      !$omp do collapse (2)
      do k = 1, n(3)
      do i = 1, n(1)
      do j = 1, n(2) - 1
      
        call Convective_Flux ( dl(i,j-1:j+2,k) % c(2), &
                               Dir(2) % f(i,j,k) % n,  &
                               Dir(2) % f(i,j,k) % a,  &
                               Prim(:,i,j-1:j+2,k),    &
                               Res (:,i,j:j+1,k)       )

        if (model>0) then
          call Diffusive_Flux ( Dir(2) % f(i,j,k) % n,  &
                                Dir(2) % f(i,j,k) % a,  &
                                yn(i,j  ,k),            &
                                yn(i,j+1,k),            &
                                ks(i,j  ,k),            &
                                ks(i,j+1,k),            &
                                Prim(:,i,j  ,k),        &
                                Prim(:,i,j+1,k),        &
                                Prim(:,i-1,j  ,k),      &
                                Prim(:,i+1,j  ,k),      &
                                Prim(:,i-1,j+1,k),      &
                                Prim(:,i+1,j+1,k),      &
                                Prim(:,i,j  ,k-1),      &
                                Prim(:,i,j  ,k+1),      &
                                Prim(:,i,j+1,k-1),      &
                                Prim(:,i,j+1,k+1),      &
                                M(i,j  ,k) % c,         &
                                M(i,j+1,k) % c,         &
                                Res (:,i,j  ,k),        &
                                Res (:,i,j+1,k),        &
                                2, 1, 3, Sc, Sct, Prt, Prt_corr )
        end if

      enddo; enddo; enddo

      !$omp do collapse (2)
      do j = 1, n(2)
      do i = 1, n(1)
      do k = 1, n(3) - 1
      
        call Convective_Flux ( dl(i,j,k-1:k+2) % c(3), &
                               Dir(3) % f(i,j,k) % n,  &
                               Dir(3) % f(i,j,k) % a,  &
                               Prim(:,i,j,k-1:k+2),    &
                               Res (:,i,j,k:k+1)       )

        if (model>0) then
          call Diffusive_Flux ( Dir(3) % f(i,j,k) % n,  &
                                Dir(3) % f(i,j,k) % a,  &
                                yn(i,j,k  ),            &
                                yn(i,j,k+1),            &
                                ks(i,j,k  ),            &
                                ks(i,j,k+1),            &
                                Prim(:,i,j,k  ),        &
                                Prim(:,i,j,k+1),        &
                                Prim(:,i-1,j,k  ),      &
                                Prim(:,i+1,j,k  ),      &
                                Prim(:,i-1,j,k+1),      &
                                Prim(:,i+1,j,k+1),      &
                                Prim(:,i,j-1,k  ),      &
                                Prim(:,i,j+1,k  ),      &
                                Prim(:,i,j-1,k+1),      &
                                Prim(:,i,j+1,k+1),      &
                                M(i,j,k  ) % c,         &
                                M(i,j,k+1) % c,         &
                                Res (:,i,j,k  ),        &
                                Res (:,i,j,k+1),        &
                                3, 1, 2, Sc, Sct, Prt, Prt_corr )
        end if

      enddo; enddo; enddo

    end subroutine Fluxes_blk

end module ARES_Mod_Fluxes