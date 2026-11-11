module ARES_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

contains

  subroutine Setup_Metrics ( domain )
    use ARES_Base_Types_m
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Lib_Metrics
    use ARES_Lib_BC_Fluxes_Rotational, only: Rotational_Periodic_Face_Metric
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k, g
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(ARES_vector_3D_type) :: N1, N2, N3, N4, N5, N6, N7, N8

    call Check_Mesh_Type ( domain )

    ! A cell folded at one of its corners, or turned against its block, stops the run here
    call Check_Folded_Cells ( domain )

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      ! Compute metric tensor and cell dimension across i,j,k. => block % M, & block % dl
      !$omp do collapse(3) private(i, j, k, N1, N2, N3, N4, N5, N6, N7, N8)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        N1 % c = domain % blk(b) % node(i-1,j-1,k-1) % c
        N2 % c = domain % blk(b) % node(i-1,j-1,k  ) % c
        N3 % c = domain % blk(b) % node(i-1,j  ,k-1) % c
        N4 % c = domain % blk(b) % node(i-1,j  ,k  ) % c
        N5 % c = domain % blk(b) % node(i  ,j-1,k-1) % c
        N6 % c = domain % blk(b) % node(i  ,j-1,k  ) % c
        N7 % c = domain % blk(b) % node(i  ,j  ,k-1) % c
        N8 % c = domain % blk(b) % node(i  ,j  ,k  ) % c
        call Compute_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, domain % blk(b) % M(i,j,k), domain % blk(b) % dl(i,j,k), domain % blk(b) % vol(i,j,k), &
                                     loc = [b, i, j, k, 0] )
      enddo; enddo; enddo
      !$omp end parallel
      ! A degenerate or inverted cell (Compute_Metric_Tensor wrote the first one): no further block here; every rank
      ! stops together below
      if ( bad_metric ) exit

      ! Compute Normal & Area
      call Compute_Norm_Area ( domain % blk(b) )
    enddo
    call Stop_On_Bad_Metric ( )

    ! A viscous run on the 2-D axisymmetric wedge or on a 3-D sector stops when the axis is exactly at r = 0
    ! (Check_Wedge_Axis)
    if ( model > 0 ) call Check_Wedge_Axis ( domain )

    if (model>1) then

      ! Compute wall minimum distance => block%yn
      call Compute_Yn ( domain % blk(:), domain % bc, domain % nb, domain % nbound )

      !$omp parallel
      ! Set wall-distance in connected cells
      !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks)
      do i = 1, domain % nbound
        select case ( domain % bc(i) % type )
          case(101,201) ! block connection
            Bm = domain % bc(i) % b
            Im = domain % bc(i) % i
            Jm = domain % bc(i) % j
            Km = domain % bc(i) % k
            Fm = domain % bc(i) % f
            Bs = domain % bc(i) % bs
            Is = domain % bc(i) % is
            Js = domain % bc(i) % js
            Ks = domain % bc(i) % ks
            call Yn_Connection ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                 Is, Js, Ks, domain % blk(Bs) )
          case(102) ! chimera: no donor copy, extend the receiver
            Bm = domain % bc(i) % b
            if (.not. is_local_block(Bm)) cycle
            Im = domain % bc(i) % i
            Jm = domain % bc(i) % j
            Km = domain % bc(i) % k
            Fm = domain % bc(i) % f
            call Yn_Chimera ( Im, Jm, Km, Fm, domain % blk(Bm) )
        end select
      enddo
      !$omp end parallel
    endif

    !$omp parallel private(g)
    ! Create the nodes for gc layers of ghost cell (of the symmetry faces the first layer only: below)
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, domain % nbound
      Bm = domain % bc(i) % b
      ! Ghost metrics land in blk(Bm) and bc(i), read only by the rank owning Bm, which alone has the metrics of the
      ! cells of Bm that they start from
      if (.not. is_local_block(Bm)) cycle
      Im = domain % bc(i) % i
      Jm = domain % bc(i) % j
      Km = domain % bc(i) % k
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type )
        case(101,201) ! block connection
          Bs = domain % bc(i) % bs
          Is = domain % bc(i) % is
          Js = domain % bc(i) % js
          Ks = domain % bc(i) % ks
          Fs = domain % bc(i) % fs
          d11s = domain % bc(i) % d11
          d12s = domain % bc(i) % d12
          d21s = domain % bc(i) % d21
          d22s = domain % bc(i) % d22
          call BC_Connect_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                    Is, Js, Ks, Fs, domain % blk(Bs), d11s, d12s, d21s, d22s, &
                                    domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, &
                                    periodic = (domain % bc(i) % type == 201) )
        case(300) ! symmetry
          call BC_Symmetry_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, 1 )
        case default  ! Rientra qui anche chimera
          call BC_Extrapolate_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, Bm, &
                                        chimera = ( domain % bc(i) % type == 102 ) )
      end select
    enddo
    !$omp end do
    ! The further layers of the symmetry faces copy the cells g - 1 behind the face, which in a direction with fewer
    ! cells than g are ghost cells of the opposite face: one layer at a time, after every face has written the previous
    ! one (the barrier at the end of each loop), so that no face reads a ghost cell that another one is writing
    do g = 2, gc
      !$omp do schedule (dynamic) private(i, Bm)
      do i = 1, domain % nbound
        if ( domain % bc(i) % type /= 300 ) cycle
        Bm = domain % bc(i) % b
        if (.not. is_local_block(Bm)) cycle
        call BC_Symmetry_Metrics ( domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, &
                                   domain % bc(i) % f, domain % blk(Bm), &
                                   domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg, g )
      enddo
      !$omp end do
    enddo
    !$omp end parallel
    ! A degenerate or inverted ghost cell (Compute_Metric_Tensor wrote the first one): every rank stops together
    call Stop_On_Bad_Metric ( )

    ! Type 200: the face metric of the viscous flux depends on the metrics of the cell and of its image (above) only:
    ! once here, into the record, for the faces of the blocks of this rank (Rotational_Periodic_Face_Metric)
    if ( model > 0 ) then
      do i = 1, domain % nbound
        if ( domain % bc(i) % type /= 200 ) cycle
        Bm = domain % bc(i) % b
        if (.not. is_local_block(Bm)) cycle
        domain % bc(i) % rM = Rotational_Periodic_Face_Metric ( domain % bc(i) % i, domain % bc(i) % j, &
                                domain % bc(i) % k, domain % bc(i) % f, domain % blk(Bm) )
      enddo
    end if

  end subroutine Setup_Metrics


  !> A viscous run on the 2-D axisymmetric wedge (type 200 on the k faces, one cell in k), or on a 3-D sector closed by
  !> type 200 on its k faces, stops at setup when a boundary face of a block so closed has zero area, that is when the
  !> axis is exactly at r = 0, or when its nodes of the two k planes coincide (axis nodes at a small radius that are not
  !> turned with the other nodes of their k plane). Behind a face of zero area the ghost cell is a copy of the cell
  !> (Ghost_Symmetry with a zero normal). The faces of type 200 (BC_Rotational_Periodic_Visc) take the derivative in r
  !> across the axis with the ghost cell one cell away: they see half of d u_r/dr, and the hoop stress of the first row of
  !> cells is wrong by 1/6 of it in a Homann flow (on the wedge, 32 and 64 cells, and on a sector of four cells of 1
  !> degree), whatever the cell size. Put the axis nodes at a small radius, 1e-6 of the radial size of the first cell,
  !> turned with the other nodes of their k plane: there the face on the axis has an area, its ghost cell is the mirror
  !> image of the cell, and these derivatives are right.
  subroutine Check_Wedge_Axis ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: ndir
    use ARES_Lib_Metrics, only: delthe
    use ARES_Lib_BC_Fluxes, only: Face_Index
    use ARES_Mod_MPI, only: is_local_block, mpi_abort_all
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    ! Local
    integer :: r, b, dir, fi, fj, fk
    logical, allocatable :: closed(:)

    ! The blocks closed by type 200: none when delthe is exactly zero (two ordered comparisons, as in
    ! BC_Rotational_Periodic_Visc); in 2-D every block is then a wedge, in 3-D delthe is not zero only with sectors turned
    ! about the x axis (Check_Mesh_Type), the blocks with a type-200 record on faces 5 and 6, as there
    if ( delthe >= 0d0 .and. delthe <= 0d0 ) return
    allocate ( closed(domain % nb) )
    closed = ( ndir < 3 )
    if ( ndir == 3 ) then
      do r = 1, domain % nbound
        if ( domain % bc(r) % type == 200 .and. domain % bc(r) % f >= 5 ) closed(domain % bc(r) % b) = .true.
      enddo
    end if
    do r = 1, domain % nbound
      b = domain % bc(r) % b
      if ( .not. is_local_block(b) .or. .not. closed(b) ) cycle
      call Face_Index ( domain % bc(r) % f, dir, domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k, fi, fj, fk )
      if ( domain % blk(b) % dir(dir) % f(fi,fj,fk) % a > 0d0 ) cycle
      if ( ndir == 3 ) then
        write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] block ', b, ', cell', domain % bc(r) % i, domain % bc(r) % j, &
                                              domain % bc(r) % k, ', face ', domain % bc(r) % f, &
                                              ': a face of zero area on the 3-D sector (type 200 on the k faces):'
        write(*,'(A)') ' [ERROR]   the axis nodes are exactly at r = 0, or those of the two k planes coincide (not turned'
        write(*,'(A)') ' [ERROR]   with the other nodes of their k plane). The ghost cells behind the faces of zero area copy'
        write(*,'(A)') ' [ERROR]   the cells, and the hoop stress of the first row of cells would be wrong (by 1/6 of it in a'
        write(*,'(A)') ' [ERROR]   Homann flow). Put the axis nodes at a small radius, 1e-6 of the radial size of the first'
        write(*,'(A)') ' [ERROR]   cell, turned with the other nodes of their k plane.'
        call mpi_abort_all ( 'axis of the 3-D sector at r = 0 (see above)' )
      else
        write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] block ', b, ', cell', domain % bc(r) % i, domain % bc(r) % j, &
                                              domain % bc(r) % k, ', face ', domain % bc(r) % f, &
                                              ': a face of zero area on the 2-D axisymmetric wedge (type 200):'
        write(*,'(A)') ' [ERROR]   the axis nodes are exactly at r = 0, or those of the two k planes coincide (not turned'
        write(*,'(A)') ' [ERROR]   with the other nodes of their k plane). The ghost cells behind the faces of zero area copy'
        write(*,'(A)') ' [ERROR]   the cells, and the hoop stress of the first row of cells would be wrong (by 1/6 of it in a'
        write(*,'(A)') ' [ERROR]   Homann flow). Put the axis nodes at a small radius, 1e-6 of the radial size of the first'
        write(*,'(A)') ' [ERROR]   cell, turned with the other nodes of their k plane.'
        call mpi_abort_all ( 'axis of the 2-D axisymmetric wedge at r = 0 (see above)' )
      end if
    enddo

  end subroutine Check_Wedge_Axis


  !> A cell is the image of the unit cube by the trilinear map of its eight nodes, and a finite volume needs that map
  !> one-to-one: its Jacobian keeps the sign of the block in the whole cell. At a corner the Jacobian is
  !> J = e_i . (e_j x e_k), with e_i, e_j and e_k the three edges from that corner along i, j and k. The run stops at
  !> setup when, at a corner of a cell of a local block, s J <= 1e-12 |e_i| |e_j| |e_k| (Hadamard's bound, the measure of
  !> the determinant check of Compute_Metric_Tensor), with s the sign of the determinant of the mean edges of the first
  !> cell of the block, the cell from which Compute_Norm_Area takes the orientation of the normals; a degenerate first
  !> cell, which gives the block no orientation, stops the run there, as that determinant check would. A corner with an
  !> edge of zero length (nodes that coincide on an axis) is not checked. On a 2-D mesh, plane or a wedge whose nodes are
  !> off the axis, the corners decide for the whole cell; in 3-D they are a necessary condition (docs/theory/numerics.md).
  !> The volume of a cell, the integral of that Jacobian, has the sign of the block when the cell is not turned over, and
  !> Compute_Metric_Tensor stops a cell whose volume has not the sign of the determinant at its centre. Setup_Metrics calls
  !> this on every multigrid level, whose nodes are every other node of the level above: a stop there names the level.
  subroutine Check_Folded_Cells ( domain )
    use ARES_Advanced_Types_m
    use ARES_Mod_MPI, only: is_local_block, mpi_abort_all
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    ! Local
    integer :: b, i, j, k, ca, cb, cc, n(3)
    real(R8) :: s, det, hm, ei(3), ej(3), ek(3), h, Jv, Am(3,3)
    character(len=32) :: lev

    lev = ':'
    if ( domain % mg_level > 1 ) write(lev,'(A,I0,A)') ' of multigrid level ', domain % mg_level, ':'

    do b = 1, domain % nb
      if ( .not. is_local_block(b) ) cycle
      n = domain % blk(b) % dim
      associate ( node => domain % blk(b) % node )
        ! Orientation of the block: the determinant of the mean edges of its first cell, as Compute_Metric_Tensor; a
        ! degenerate first cell gives none, and stops the run here as it would there
        Am(1,:) = node(1,0,0)%c - node(0,0,0)%c + node(1,0,1)%c - node(0,0,1)%c &
                + node(1,1,0)%c - node(0,1,0)%c + node(1,1,1)%c - node(0,1,1)%c
        Am(2,:) = node(0,1,0)%c - node(0,0,0)%c + node(0,1,1)%c - node(0,0,1)%c &
                + node(1,1,0)%c - node(1,0,0)%c + node(1,1,1)%c - node(1,0,1)%c
        Am(3,:) = node(0,0,1)%c - node(0,0,0)%c + node(0,1,1)%c - node(0,1,0)%c &
                + node(1,0,1)%c - node(1,0,0)%c + node(1,1,1)%c - node(1,1,0)%c
        Am = 0.25d0 * Am
        det = Triple ( Am(1,:), Am(2,:), Am(3,:) )
        hm  = norm2 ( Am(1,:) ) * norm2 ( Am(2,:) ) * norm2 ( Am(3,:) )
        if ( .not. ( abs(det) > 1d-12 * hm ) ) then
          write(*,'(A,I0,A,ES12.4,A,ES12.4)') ' [ERROR] Mesh'//trim(lev)//' block ', b, &
            ', cell 1 1 1: degenerate, so the block has no orientation; det = ', det, ', |A1| |A2| |A3| = ', hm
          call mpi_abort_all ( 'degenerate cell in the mesh (see above)' )
        end if
        s = sign ( 1d0, det )
        do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)
          do cc = 0, 1; do cb = 0, 1; do ca = 0, 1
            ei = node(i,        j-1+cb, k-1+cc)%c - node(i-1,      j-1+cb, k-1+cc)%c
            ej = node(i-1+ca,   j,      k-1+cc)%c - node(i-1+ca,   j-1,    k-1+cc)%c
            ek = node(i-1+ca,   j-1+cb, k     )%c - node(i-1+ca,   j-1+cb, k-1   )%c
            h = norm2 ( ei ) * norm2 ( ej ) * norm2 ( ek )
            if ( h == 0d0 ) cycle   ! an edge of zero length: nodes that coincide on an axis
            Jv = Triple ( ei, ej, ek )
            if ( s * Jv > 1d-12 * h ) cycle
            write(*,'(A,I0,A,3(1X,I0),A,3(1X,I0),A)') ' [ERROR] Mesh'//trim(lev)//' block ', b, ', cell', i, j, k, &
              ', corner node', i-1+ca, j-1+cb, k-1+cc, &
              ': the cell is folded or turned against its block'
            write(*,'(A,ES12.4,A,ES12.4,A,I0)') ' [ERROR] e_i . (e_j x e_k) = ', Jv, ', |e_i| |e_j| |e_k| = ', h, &
              ', sign of the block (from its first cell): ', nint(s)
            write(*,'(A)') ' [ERROR] At every corner of a cell the edges along i, j and k must keep the orientation of the block.'
            write(*,'(A,3(1X,ES14.6))') ' [ERROR] corner node:', node(i-1+ca, j-1+cb, k-1+cc)%c
            call mpi_abort_all ( 'folded cell in the mesh (see above)' )
          enddo; enddo; enddo
        enddo; enddo; enddo
      end associate
    enddo

  contains

    pure function Triple ( u, v, w ) result ( t )
      real(R8), intent(in) :: u(3), v(3), w(3)
      real(R8) :: t
      t = u(1) * ( v(2)*w(3) - v(3)*w(2) ) + u(2) * ( v(3)*w(1) - v(1)*w(3) ) + u(3) * ( v(1)*w(2) - v(2)*w(1) )
    end function Triple

  end subroutine Check_Folded_Cells

end module ARES_Mod_Metrics
