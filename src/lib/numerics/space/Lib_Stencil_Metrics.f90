!> Face metrics of the viscous fluxes with diffusive-metric = stencil.
!>
!> The face gradient takes, for each quantity, one difference across the face (the two cells) and one per tangential
!> direction (the mean of the central differences of the two cells), D phi = (D_a, D_b, D_c), and maps them to Cartesian
!> components with a face metric: grad phi = M_f^T D phi. Here M_f comes from the same differences of the cell positions:
!> J^T = [D_a x, D_b x, D_c x] (columns), M_f = (J^T)^-1. For a linear field phi = phi0 + g . x every difference is
!> D phi = J g, so the face gradient is g on any mesh: skewed, stretched, with jumps of size, folded at the faces between
!> blocks. "In order to obtain the identity tensor numerically ... the two difference approximations must be the same"
!> (M. Vinokur, An analysis of finite-difference and finite-volume formulations of conservation laws, NASA CR-177416,
!> 1986, p. 26; J. Comput. Phys. 81, 1989). The inverse-mean metric of the two cells, [1/2 (A_L + A_R)]^-T, has
!> the right row across the face but not along it: its tangential rows are the mean edges of the two cells, not the
!> differences of the centres that D_b and D_c take, so on skewed, non-smooth meshes the face gradient of a linear field
!> is wrong.
!>
!> Rows kept from the inverse-mean metric, so that its bits do not change where nothing has to: the row across the face
!> (1/2 (a_L + a_R) is x_R - x_L with the centre as the mean of the eight nodes), and the rows of the index directions with
!> one cell (a 1-D case, the k direction of a 2-D mesh or of a wedge). On a 1-D case and a plane 2-D mesh the field does
!> not vary along them, D phi = 0. On a wedge the k neighbours are the cell turned by -delthe and +delthe, so the velocity
!> differs across them, and that row, the mean edge of the cells, is the difference of the turned positions divided by
!> cos^2(delthe/2): the error of order delthe^2 of the azimuthal derivatives (docs/theory/numerics.md, Axisymmetric wedge).
!>
!> Positions: the centre of a cell is the mean of its eight nodes; the first layer of ghost cells takes the position at
!> which its value is right for a linear field: the centre of the cell across a face between blocks (plus the translation of
!> a periodic pair), behind a wall (103, 104, 301, 302) 2 x_f - x1 (x_f the centre of the face) for the quantities with a
!> wall value and x1 + s (x1 - x2) for the extrapolated ones, behind a 400 face x1 + s (x1 - x2) (Ghost_Wall; s = l0/l1,
!> x1 with one cell across the face), 3 x1 - 3 x2 + x3 behind the faces whose ghost value is 3 u1 - 3 u2 + u3
!> (Ghost_Extrapolate; 2 x1 - x2, or x1, in a block with two cells, or one, across the face), the cell itself behind the
!> zero-gradient faces (0, 404-406), its mirror image across a symmetry face (300; the cell itself on a face of zero area,
!> an axis, where the ghost value is a copy of the cell), and on the k faces of a wedge or of a 3-D sector (200) the source
!> of the ghost cell, the cell next to the opposite k face of the block (on a wedge the cell itself), rotated by the angle
!> of the wedge or sector, or, on a translated sector, translated from its outer face to this face as across a periodic
!> pair. Chimera (102) faces take the extrapolated position too, an approximation: their values come from donor cells.
!> Behind a wall the two positions of a ghost cell give the faces with that ghost cell in their stencil one metric per group
!> of quantities (1 unused: a single species, 2 temperature, 3 velocity, 4 turbulence variables: dir(d)%Mg,
!> Diffusive_Flux). The centres of every block are computed on every rank, from the nodes every rank holds until the end of
!> the setup of the metrics, so no message is needed; the face metrics are stored for the blocks of the rank only
!> (blk%dir(d)%Mf, on the faces of the block, interior and boundary).
!>
!> Walls: the gradient at a wall face takes the difference between the cell and the wall value across the face and the
!> differences of the wall values along it, with the metric of the same differences of the positions (cell centre, face
!> centres): Wall_Face_Stencils. With mean the differences along the wall are zero and the metric is that of the cell,
!> wrong by (dphi_w/ds) tan(theta) on a wall whose value varies along it, on cells skewed by theta.
module ARES_Lib_Stencil_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Stencil_Metrics, Source_Metrics, METRIC_MEAN, METRIC_STENCIL

  integer, parameter :: METRIC_MEAN = 0, METRIC_STENCIL = 2

  !> Cell centres of one block with one layer of ghost cells: x(1:3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc); behind the walls
  !> and the 400 faces also a second position, xe (that of the extrapolated quantities), and the groups of quantities that
  !> take it, em (bit g-1 for group g: 1 species and scalars, unused with a single species, 2 temperature, 3 velocity,
  !> 4 turbulence variables)
  type :: centres_t
    real(R8), allocatable :: x(:,:,:,:), xe(:,:,:,:)
    integer, allocatable  :: em(:,:,:)
  end type centres_t

  integer, parameter :: NG = 4, ALL_EXTRAPOLATED = 15   ! groups of quantities of the viscous fluxes; em with every bit set

contains

  !> Centre of cell (i,j,k) of a block: the mean of its eight nodes.
  pure function Cell_Centre ( node, i, j, k ) result ( x )
    use ARES_Base_Types_m, only: ARES_vector_3D_type
    implicit none
    type(ARES_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer, intent(in) :: i, j, k
    real(R8) :: x(3)

    x = 0.125d0 * ( node(i-1,j-1,k-1) % c + node(i,j-1,k-1) % c + node(i-1,j,k-1) % c + node(i,j,k-1) % c &
                  + node(i-1,j-1,k  ) % c + node(i,j-1,k  ) % c + node(i-1,j,k  ) % c + node(i,j,k  ) % c )

  end function Cell_Centre


  !> Centre (mean of the four nodes) and unit normal (direction of growing index) of face Fm of cell (i,j,k). The normal
  !> of a face of zero area (its nodes on an axis) is zero, as in Compute_Norm_Area: Mirror then gives the point itself,
  !> as Ghost_Symmetry, which reflects the velocity with that normal, copies the cell there.
  pure subroutine Face_Centre_Normal ( node, i, j, k, Fm, xf, nf )
    use ARES_Base_Types_m, only: ARES_vector_3D_type
    implicit none
    type(ARES_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer, intent(in) :: i, j, k, Fm
    real(R8), intent(out) :: xf(3), nf(3)
    ! Local
    real(R8) :: p1(3), p2(3), p3(3), p4(3), d1(3), d2(3), v(3)
    integer :: a

    select case (Fm)
      case (1:2)
        a = i - mod(Fm,2)
        p1 = node(a,j-1,k-1) % c ; p2 = node(a,j,k-1) % c ; p3 = node(a,j,k) % c ; p4 = node(a,j-1,k) % c
      case (3:4)
        a = j - mod(Fm,2)
        p1 = node(i-1,a,k-1) % c ; p2 = node(i,a,k-1) % c ; p3 = node(i,a,k) % c ; p4 = node(i-1,a,k) % c
      case default
        a = k - mod(Fm,2)
        p1 = node(i-1,j-1,a) % c ; p2 = node(i,j-1,a) % c ; p3 = node(i,j,a) % c ; p4 = node(i-1,j,a) % c
    end select
    xf = 0.25d0 * ( p1 + p2 + p3 + p4 )
    d1 = p3 - p1
    d2 = p4 - p2
    v = [ d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1) ]
    if ( norm2(v) > 0d0 ) then
      nf = v / norm2(v)
    else
      nf = 0d0
    end if

  end subroutine Face_Centre_Normal


  !> Mirror image of x across the plane of a face (centre xf, unit normal nf).
  pure function Mirror ( x, xf, nf ) result ( y )
    implicit none
    real(R8), intent(in) :: x(3), xf(3), nf(3)
    real(R8) :: y(3)

    y = x - 2d0 * dot_product ( x - xf, nf ) * nf

  end function Mirror


  !> Positions of the cell centres (every block) and face metrics (blocks of this rank) for diffusive-metric = stencil.
  !> Called by Setup_Metrics for every multigrid level, after the ghost metrics.
  subroutine Setup_Stencil_Metrics ( domain )
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, model
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    type(centres_t), allocatable :: xc(:)
    integer :: b, i, j, k, n(3), d, r
    integer, allocatable :: nfill(:,:)

    allocate ( xc(domain % nb) )

    ! 1. Centres of the cells of every block
    do b = 1, domain % nb
      n = domain % blk(b) % dim
      allocate ( xc(b) % x(3, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc) )
      allocate ( xc(b) % xe, mold = xc(b) % x )
      allocate ( xc(b) % em(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc) )
      xc(b) % x = 0d0
      xc(b) % xe = 0d0
      xc(b) % em = 0
      do k = 1, n(3) ; do j = 1, n(2) ; do i = 1, n(1)
        xc(b) % x(:,i,j,k) = Cell_Centre ( domain % blk(b) % node, i, j, k )
      enddo ; enddo ; enddo
    enddo

    ! 2. Ghost cells of the faces that are not between blocks, every block
    do r = 1, domain % nbound
      select case ( domain % bc(r) % type )
        case (101, 201)
          cycle
        case default
          call Ghost_Position ( domain, xc, r )
      end select
    enddo

    ! 3. Ghost cells of the faces between blocks: the centre of the cell across (computed in 1, not a ghost)
    do r = 1, domain % nbound
      select case ( domain % bc(r) % type )
        case (101, 201)
          call Connection_Ghost_Position ( domain, xc, r )
      end select
    enddo

    ! 4. Face metrics of the blocks of this rank; with a ghost cell at two positions in the stencil (behind a wall), one
    !    per group of quantities
    allocate ( nfill(3, domain % nb) )
    nfill = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n = domain % blk(b) % dim
      if ( .not. allocated(domain % blk(b) % dir(1) % Mf) ) then
        allocate ( domain % blk(b) % dir(1) % Mf(0:n(1), n(2), n(3)) )
        allocate ( domain % blk(b) % dir(2) % Mf(n(1), 0:n(2), n(3)) )
        allocate ( domain % blk(b) % dir(3) % Mf(n(1), n(2), 0:n(3)) )
        allocate ( domain % blk(b) % dir(1) % ig(0:n(1), n(2), n(3)) )
        allocate ( domain % blk(b) % dir(2) % ig(n(1), 0:n(2), n(3)) )
        allocate ( domain % blk(b) % dir(3) % ig(n(1), n(2), 0:n(3)) )
      end if
      do d = 1, 3
        domain % blk(b) % dir(d) % Mf(:,:,:) = ARES_tensor_3D_type(0d0)
        domain % blk(b) % dir(d) % ig = 0
        if ( allocated(domain % blk(b) % dir(d) % Mg) ) deallocate ( domain % blk(b) % dir(d) % Mg )
        allocate ( domain % blk(b) % dir(d) % Mg(3, 3, NG, Group_Faces(domain, xc, b, d)) )
      enddo
      call Interior_Face_Metrics ( b, domain % blk(b), xc(b), nfill(:,b) )
    enddo
    do r = 1, domain % nbound
      b = domain % bc(r) % b
      if (.not. is_local_block(b)) cycle
      select case ( domain % bc(r) % type )
        case (101, 201)
          call Connection_Face_Metric ( domain, xc, r, nfill(:,b) )
        case (200)
          call Rotational_Face_Metric ( domain, xc, r, nfill(:,b) )
        case (300)
          call Symmetry_Face_Metric ( domain % blk(b), xc(b), domain % bc(r) )
      end select
    enddo

    ! 5. Wall faces of the blocks of this rank
    call Wall_Face_Stencils ( domain, xc )

    ! 6. Cell metrics of the turbulence source terms (RANS runs)
    if ( model > 1 ) call Cell_Source_Metrics ( domain, xc )

    deallocate ( xc )

  end subroutine Setup_Stencil_Metrics


  !> Position of the first ghost cell of boundary record r (not a face between blocks).
  subroutine Ghost_Position ( domain, xc, r )
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: ndir
    use ARES_Parameters_m, only: guide
    use ARES_Lib_Metrics, only: delthe
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    type(centres_t), intent(inout) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, Im, Jm, Km, Fm, Ig, Jg, Kg, d, nd, Ks
    integer :: c1(3), c2(3), c3(3)
    real(R8) :: x1(3), x2(3), x3(3), xf(3), nf(3), y(3), ang

    b  = domain % bc(r) % b
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    Ig = Im - guide(Fm,1) ; Jg = Jm - guide(Fm,2) ; Kg = Km - guide(Fm,3)
    d  = (Fm + 1) / 2
    nd = domain % blk(b) % dim(d)
    c1 = [Im, Jm, Km]
    c2 = c1 + guide(Fm,:)
    c3 = c2 + guide(Fm,:)
    x1 = xc(b) % x(:,c1(1),c1(2),c1(3))

    select case ( domain % bc(r) % type )

      case (0, 404:406)                         ! zero gradient: the ghost value is the cell value
        xc(b) % x(:,Ig,Jg,Kg) = x1

      case (300)                                ! symmetry: mirror image across the face
        call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Km, Fm, xf, nf )
        xc(b) % x(:,Ig,Jg,Kg) = Mirror ( x1, xf, nf )

      case (200)                                ! k faces of a wedge or of a 3-D sector, turned or translated
        if ( Fm >= 5 .and. ( ndir == 3 .or. .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) ) then
          ! The source of the ghost cell (Ghost_Axysimmetry), the cell next to the opposite k face of the block (on a wedge,
          ! one cell in k, the cell itself), moved by the periodicity
          Ks = domain % blk(b) % dim(3) - Km + 1
          y = xc(b) % x(:,Im,Jm,Ks)
          if ( delthe >= 0d0 .and. delthe <= 0d0 ) then
            ! Translated sector (no angle: Check_Mesh_Type found the last k node plane of the block to be the first one
            ! translated, node by node): the cell of the other end plus the translation from its outer face (11 - Fm) to
            ! this face, the same operations as for a periodic pair (Connection_Ghost_Position)
            call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Km, Fm, xf, nf )
            call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Ks, 11 - Fm, x2, nf )
            xc(b) % x(:,Ig,Jg,Kg) = y + ( xf - x2 )
          else
            ang = delthe
            if ( Fm == 5 ) ang = -delthe
            xc(b) % x(:,Ig,Jg,Kg) = [ y(1), y(2)*cos(ang) - y(3)*sin(ang), y(2)*sin(ang) + y(3)*cos(ang) ]
          end if
        else
          xc(b) % x(:,Ig,Jg,Kg) = x1             ! a plane 2-D mesh, where these faces carry no viscous flux
        end if

      case (103, 104, 301, 302, 400)            ! walls and 400 (Ghost_Wall): x1 + s (x1 - x2) for the extrapolated
        y = x1                                  ! quantities (x1 with one cell), 2 x_f - x1 for those with a wall value
        if ( nd >= 2 ) y = x1 + domain % bc(r) % sx * ( x1 - xc(b) % x(:,c2(1),c2(2),c2(3)) )
        xc(b) % xe(:,Ig,Jg,Kg) = y
        xc(b) % em(Ig,Jg,Kg) = Extrapolated_Groups ( domain % bc(r) % type )
        if ( xc(b) % em(Ig,Jg,Kg) == ALL_EXTRAPOLATED ) then
          xc(b) % x(:,Ig,Jg,Kg) = y
        else
          call Face_Centre_Normal ( domain % blk(b) % node, Im, Jm, Km, Fm, xf, nf )
          xc(b) % x(:,Ig,Jg,Kg) = 2d0 * xf - x1
        end if

      case default                              ! 3 u1 - 3 u2 + u3 (Ghost_Extrapolate), chimera
        if ( nd >= 3 ) then
          x2 = xc(b) % x(:,c2(1),c2(2),c2(3))
          x3 = xc(b) % x(:,c3(1),c3(2),c3(3))
          xc(b) % x(:,Ig,Jg,Kg) = 3d0 * x1 - 3d0 * x2 + x3
        elseif ( nd == 2 ) then
          x2 = xc(b) % x(:,c2(1),c2(2),c2(3))
          xc(b) % x(:,Ig,Jg,Kg) = 2d0 * x1 - x2
        else
          xc(b) % x(:,Ig,Jg,Kg) = x1
        end if
    end select

  end subroutine Ghost_Position


  !> Position of the first ghost cell of a face between blocks (101) or of a periodic pair (201): the centre of the source
  !> cell, plus the translation from the source face to this face for 201.
  subroutine Connection_Ghost_Position ( domain, xc, r )
    use ARES_Advanced_Types_m
    use ARES_Parameters_m, only: guide
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    type(centres_t), intent(inout) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, bs, Ig, Jg, Kg

    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    Ig = domain % bc(r) % i - guide(domain % bc(r) % f,1)
    Jg = domain % bc(r) % j - guide(domain % bc(r) % f,2)
    Kg = domain % bc(r) % k - guide(domain % bc(r) % f,3)
    xc(b) % x(:,Ig,Jg,Kg) = xc(bs) % x(:,domain % bc(r) % is, domain % bc(r) % js, domain % bc(r) % ks) &
                          + Periodic_Shift ( domain, r )

  end subroutine Connection_Ghost_Position


  !> Groups of quantities extrapolated behind a face of type t (Ghost_Wall), as the bits of em (1 species and scalars, unused
  !> with a single species, 2 temperature, 3 velocity, 4 turbulence variables); the others take the wall value.
  pure integer function Extrapolated_Groups ( t )
    implicit none
    integer, intent(in) :: t

    select case (t)
      case (103, 104, 301)
        Extrapolated_Groups = 3                  ! species, temperature (103, 104: the walls coupled to another solver)
      case (302)
        Extrapolated_Groups = 1                  ! species
      case default
        Extrapolated_Groups = ALL_EXTRAPOLATED   ! 400
    end select

  end function Extrapolated_Groups


  !> Position of the value of group g in cell (i,j,k) of a block: the centre, behind a wall the position of the group.
  pure function Pos ( xc, i, j, k, g ) result ( y )
    implicit none
    type(centres_t), intent(in) :: xc
    integer, intent(in) :: i, j, k, g
    real(R8) :: y(3)

    if ( btest(xc % em(i,j,k), g-1) ) then
      y = xc % xe(:,i,j,k)
    else
      y = xc % x(:,i,j,k)
    end if

  end function Pos


  !> A ghost cell whose groups of quantities sit at two positions: behind a wall (a wall ghost cell below).
  pure logical function Wall_Ghost ( xc, i, j, k )
    implicit none
    type(centres_t), intent(in) :: xc
    integer, intent(in) :: i, j, k

    Wall_Ghost = xc % em(i,j,k) /= 0 .and. xc % em(i,j,k) /= ALL_EXTRAPOLATED

  end function Wall_Ghost


  !> Interior face of direction d between cell c and the next one along d: a wall ghost cell among the tangential
  !> neighbours that its face gradient takes.
  pure logical function Interior_Mixed ( xc, n, d, c ) result ( mixed )
    implicit none
    type(centres_t), intent(in) :: xc
    integer, intent(in) :: n(3), d, c(3)
    ! Local
    integer :: t, s, e(3), c2(3), a(3)

    mixed = .false.
    c2 = c
    c2(d) = c(d) + 1
    do t = 1, 3
      if ( t == d .or. n(t) == 1 ) cycle
      e = 0
      e(t) = 1
      do s = -1, 1, 2
        a = c + s*e
        if ( Wall_Ghost(xc, a(1), a(2), a(3)) ) mixed = .true.
        a = c2 + s*e
        if ( Wall_Ghost(xc, a(1), a(2), a(3)) ) mixed = .true.
      enddo
    enddo

  end function Interior_Mixed


  !> Face of record r between blocks: a wall ghost cell among the tangential neighbours its face gradient takes, in this
  !> block or among the cells of the source block of Pg(:,3:6).
  logical function Connection_Mixed ( domain, xc, r ) result ( mixed )
    use ARES_Advanced_Types_m
    use ARES_Lib_BC_Fluxes, only: Face_Index
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, bs, c(3), dir, fi, fj, fk, q, t, tb(2), sc(3,4), n(3), e(3)

    mixed = .false.
    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    c  = [domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k]
    n  = domain % blk(b) % dim
    call Face_Index ( domain % bc(r) % f, dir, c(1), c(2), c(3), fi, fj, fk )
    select case (dir)
      case (1) ; tb = [2, 3]
      case (2) ; tb = [1, 3]
      case default ; tb = [1, 2]
    end select
    call Source_Neighbours ( domain % bc(r), sc )
    do q = 1, 2
      t = tb(q)
      if ( n(t) == 1 ) cycle
      e = 0
      e(t) = 1
      if ( Wall_Ghost(xc(b), c(1)+e(1), c(2)+e(2), c(3)+e(3)) ) mixed = .true.
      if ( Wall_Ghost(xc(b), c(1)-e(1), c(2)-e(2), c(3)-e(3)) ) mixed = .true.
      if ( Wall_Ghost(xc(bs), sc(1,2*q-1), sc(2,2*q-1), sc(3,2*q-1)) ) mixed = .true.
      if ( Wall_Ghost(xc(bs), sc(1,2*q), sc(2,2*q), sc(3,2*q)) ) mixed = .true.
    enddo

  end function Connection_Mixed


  !> Record r is a face of type 200 whose viscous flux BC_Rotational_Periodic_Visc computes: faces 5 and 6 of a wedge (a
  !> 2-D mesh turned about the x axis), face 6 of a 3-D sector, turned or translated (it computes the face for the two
  !> cells). The same two tests as there and as Rotational_Periodic_Face_Metric.
  logical function Rotational_Flux_Face ( domain, r ) result ( flux )
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: ndir
    use ARES_Lib_Metrics, only: delthe
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    integer, intent(in) :: r

    flux = .false.
    if ( domain % bc(r) % type /= 200 .or. domain % bc(r) % f < 5 ) return
    if ( ndir < 3 .and. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) return
    if ( ndir == 3 .and. domain % bc(r) % f /= 6 ) return
    flux = .true.

  end function Rotational_Flux_Face


  !> Face of type 200 of record r: a wall ghost cell among the neighbours along the face that its gradient takes, in the row
  !> of the cell or in the row of its source (the cell next to the opposite k face, the cell itself on a wedge), whose
  !> values, turned, are those of the row of the image (BC_Rotational_Periodic_Visc).
  logical function Rotational_Mixed ( domain, xc, r ) result ( mixed )
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: r
    ! Local
    integer :: b, Im, Jm, Km, Ks, t, s, n(3), e(3)

    mixed = .false.
    b  = domain % bc(r) % b
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    n  = domain % blk(b) % dim
    Ks = n(3) - Km + 1
    do t = 1, 2
      if ( n(t) == 1 ) cycle
      e = 0
      e(t) = 1
      do s = -1, 1, 2
        if ( Wall_Ghost(xc(b), Im+s*e(1), Jm+s*e(2), Km) ) mixed = .true.
        if ( Wall_Ghost(xc(b), Im+s*e(1), Jm+s*e(2), Ks) ) mixed = .true.
      enddo
    enddo

  end function Rotational_Mixed


  !> Number of the faces of direction d of block b (interior, between blocks, and of type 200) with a wall ghost cell in
  !> their stencil.
  integer function Group_Faces ( domain, xc, b, d ) result ( nf )
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: b, d
    ! Local
    integer :: i, j, k, r, f, ii, n(3), e(3)

    nf = 0
    n = domain % blk(b) % dim
    e = 0
    e(d) = 1
    do k = 1, n(3) - e(3) ; do j = 1, n(2) - e(2) ; do i = 1, n(1) - e(1)
      if ( Interior_Mixed(xc(b), n, d, [i,j,k]) ) nf = nf + 1
    enddo ; enddo ; enddo
    do f = 2*d - 1, 2*d                                                 ! the records of the two faces of direction d
      do ii = domain % local_bf_first(b,f), domain % local_bf_first(b,f) + domain % local_bf_count(b,f) - 1
        r = domain % local_bf_idx(ii)
        select case ( domain % bc(r) % type )
          case (101, 201)
            if ( Connection_Mixed(domain, xc, r) ) nf = nf + 1
          case (200)
            if ( .not. Rotational_Flux_Face(domain, r) ) cycle
            if ( Rotational_Mixed(domain, xc, r) ) nf = nf + 1
        end select
      enddo
    enddo

  end function Group_Faces


  !> Translation of a periodic pair (201): centre of this face minus centre of the source face; zero for 101.
  function Periodic_Shift ( domain, r ) result ( t )
    use ARES_Advanced_Types_m
    implicit none
    type(ARES_domain_type), intent(in) :: domain
    integer, intent(in) :: r
    real(R8) :: t(3)
    ! Local
    real(R8) :: xfm(3), xfs(3), nf(3)
    integer :: b, bs

    t = 0d0
    if ( domain % bc(r) % type /= 201 ) return
    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    call Face_Centre_Normal ( domain % blk(b) % node, domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k, &
                              domain % bc(r) % f, xfm, nf )
    call Face_Centre_Normal ( domain % blk(bs) % node, domain % bc(r) % is, domain % bc(r) % js, domain % bc(r) % ks, &
                              domain % bc(r) % fs, xfs, nf )
    t = xfm - xfs

  end function Periodic_Shift


  !> Transposed face Jacobian of the inverse-mean metric of two cells, (A_L^T + A_R^T)/2 with A^T = M^-1: its columns
  !> are the mean edges; the column across the face is x_R - x_L.
  function Mean_Edges ( M1, M2 ) result ( JT )
    use ARES_Lib_Diffusive, only: Inverse3
    implicit none
    real(R8), intent(in) :: M1(3,3), M2(3,3)
    real(R8) :: JT(3,3)

    JT = 0.5d0 * ( Inverse3 ( M1 ) + Inverse3 ( M2 ) )

  end function Mean_Edges


  !> Metric (J^T)^-1 of face f of cell (i,j,k) of block b, or of the cell itself when f = 0 (kind: which one). The run
  !> stops when the three differences of the stencil do not span the space, |det J^T| <= 1e-12 |D_a x| |D_b x| |D_c x|:
  !> a degenerate face or cell, where the metric is not defined. A declared limit of diffusive-metric = stencil, which
  !> mean does not have, although a viscous run with mean stops on a 3-D mesh at a boundary face of zero area
  !> (Check_Wedge_Axis). (A wall on a face of zero area, as on the axis of a 3-D mesh around it, stops
  !> before, in Wall_Face_Stencils.)
  function Face_Inverse ( JT, b, i, j, k, f, kind ) result ( Mf )
    use ARES_Lib_Diffusive, only: Inverse3
    use ARES_Mod_MPI, only: mpi_abort_all
    implicit none
    real(R8), intent(in) :: JT(3,3)
    integer, intent(in) :: b, i, j, k, f
    character(len=*), intent(in) :: kind
    real(R8) :: Mf(3,3)
    ! Local
    real(R8) :: det
    integer :: t

    det = JT(1,1) * ( JT(2,2)*JT(3,3) - JT(2,3)*JT(3,2) ) &
        - JT(1,2) * ( JT(2,1)*JT(3,3) - JT(2,3)*JT(3,1) ) &
        + JT(1,3) * ( JT(2,1)*JT(3,2) - JT(2,2)*JT(3,1) )
    if ( .not. ( abs(det) > 1d-12 * norm2(JT(:,1)) * norm2(JT(:,2)) * norm2(JT(:,3)) ) ) then
      if ( f > 0 ) then
        write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] diffusive-metric = stencil: block ', b, ', cell', i, j, k, &
                                              ', face ', f, ' ('//kind//')'
      else
        write(*,'(A,I0,A,3(1X,I0),A)') ' [ERROR] diffusive-metric = stencil: block ', b, ', cell', i, j, k, ' ('//kind//')'
      end if
      write(*,'(A)') ' [ERROR]   degenerate face or cell: the metric of the stencil is not defined there. The differences'
      write(*,'(A)') ' [ERROR]   of the cell centres it takes:'
      do t = 1, 3
        write(*,'(A,I0,A,3(1X,ES14.6))') ' [ERROR]   difference ', t, ':', JT(:,t)
      enddo
      write(*,'(A)') ' [ERROR] A declared limit of diffusive-metric = stencil; diffusive-metric = mean does not'
      write(*,'(A)') ' [ERROR] need this metric (but on a 3-D mesh it stops at a boundary face of zero area).'
      call mpi_abort_all ( 'metric of diffusive-metric = stencil not defined (see above)' )
    end if
    Mf = Inverse3 ( JT )

  end function Face_Inverse


  !> Face metrics of the interior faces of block b; with a wall ghost cell in the stencil also one per group of quantities
  !> (dir(d)%Mg, column dir(d)%ig, the next of nfill(d)).
  subroutine Interior_Face_Metrics ( b, blk, xc, nfill )
    use ARES_Advanced_Types_m, only: ARES_block_type
    implicit none
    integer, intent(in) :: b
    type(ARES_block_type), intent(inout) :: blk
    type(centres_t), intent(in) :: xc
    integer, intent(inout) :: nfill(3)
    ! Local
    integer :: i, j, k, n(3), d, t, i2, j2, k2, g
    integer :: e(3,3)
    real(R8) :: JT(3,3), JTg(3,3)

    n = blk % dim
    e = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    do d = 1, 3
      do k = 1, n(3) - e(3,d) ; do j = 1, n(2) - e(2,d) ; do i = 1, n(1) - e(1,d)
        i2 = i + e(1,d) ; j2 = j + e(2,d) ; k2 = k + e(3,d)
        JT = Mean_Edges ( blk % M(i,j,k) % c, blk % M(i2,j2,k2) % c )
        do t = 1, 3
          if ( t == d .or. n(t) == 1 ) cycle
          JT(:,t) = 0.25d0 * ( ( xc % x(:,i +e(1,t),j +e(2,t),k +e(3,t)) - xc % x(:,i -e(1,t),j -e(2,t),k -e(3,t)) ) &
                             + ( xc % x(:,i2+e(1,t),j2+e(2,t),k2+e(3,t)) - xc % x(:,i2-e(1,t),j2-e(2,t),k2-e(3,t)) ) )
        enddo
        blk % dir(d) % Mf(i,j,k) % c = Face_Inverse ( JT, b, i, j, k, 2*d, 'interior' )
        if ( .not. Interior_Mixed(xc, n, d, [i,j,k]) ) cycle
        nfill(d) = nfill(d) + 1
        blk % dir(d) % ig(i,j,k) = nfill(d)
        do g = 1, NG
          JTg = JT
          do t = 1, 3
            if ( t == d .or. n(t) == 1 ) cycle
            JTg(:,t) = 0.25d0 * ( ( Pos(xc, i +e(1,t), j +e(2,t), k +e(3,t), g) - Pos(xc, i -e(1,t), j -e(2,t), k -e(3,t), g) ) &
                                + ( Pos(xc, i2+e(1,t), j2+e(2,t), k2+e(3,t), g) - Pos(xc, i2-e(1,t), j2-e(2,t), k2-e(3,t), g) ) )
          enddo
          blk % dir(d) % Mg(:,:,g,nfill(d)) = Face_Inverse ( JTg, b, i, j, k, 2*d, 'interior' )
        enddo
      enddo ; enddo ; enddo
    enddo

  end subroutine Interior_Face_Metrics


  !> Face metric of the face of boundary record r between blocks, stored on the face of block b. The two cells in the
  !> order of growing index of block b (L, R) and their tangential neighbours: in block b, and for the ghost cell the
  !> cells of the source block that Fill_BC_Ghost_Connection takes for Pg(:,3:6).
  subroutine Connection_Face_Metric ( domain, xc, r, nfill )
    use ARES_Advanced_Types_m
    use ARES_Parameters_m, only: guide
    use ARES_Lib_Diffusive, only: Inverse3
    use ARES_Lib_BC_Fluxes, only: Face_Index
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: r
    integer, intent(inout) :: nfill(3)
    ! Local
    integer :: b, bs, Im, Jm, Km, Fm, Ig, Jg, Kg, dir, fi, fj, fk, t, tb(2), q, n(3), g
    integer :: e(3,3), sc(3,4)
    real(R8) :: JT(3,3), xg(3,4), shift(3), dloc(3), dgh(3), JTg(3,3)

    b  = domain % bc(r) % b
    bs = domain % bc(r) % bs
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    Ig = Im - guide(Fm,1) ; Jg = Jm - guide(Fm,2) ; Kg = Km - guide(Fm,3)
    n  = domain % blk(b) % dim
    e  = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    call Face_Index ( Fm, dir, Im, Jm, Km, fi, fj, fk )
    select case (dir)
      case (1) ; tb = [2, 3]
      case (2) ; tb = [1, 3]
      case default ; tb = [1, 2]
    end select

    ! Cells of the source block that hold Pg(:,3:6) (Fill_BC_Ghost_Connection), and their positions
    call Source_Neighbours ( domain % bc(r), sc )
    shift = Periodic_Shift ( domain, r )
    do q = 1, 4
      xg(:,q) = xc(bs) % x(:,sc(1,q),sc(2,q),sc(3,q)) + shift
    enddo

    ! Face Jacobian: the cells in the order L, R of block b
    if ( mod(Fm,2) == 0 ) then
      JT = Mean_Edges ( domain % blk(b) % M(Im,Jm,Km) % c, domain % blk(b) % M(Ig,Jg,Kg) % c )
    else
      JT = Mean_Edges ( domain % blk(b) % M(Ig,Jg,Kg) % c, domain % blk(b) % M(Im,Jm,Km) % c )
    end if
    do q = 1, 2
      t = tb(q)
      if ( n(t) == 1 ) cycle
      dloc = xc(b) % x(:,Im+e(1,t),Jm+e(2,t),Km+e(3,t)) - xc(b) % x(:,Im-e(1,t),Jm-e(2,t),Km-e(3,t))
      dgh  = xg(:,2*q) - xg(:,2*q-1)
      if ( mod(Fm,2) == 0 ) then
        JT(:,t) = 0.25d0 * ( dloc + dgh )
      else
        JT(:,t) = 0.25d0 * ( dgh + dloc )
      end if
    enddo
    domain % blk(b) % dir(dir) % Mf(fi,fj,fk) % c = Face_Inverse ( JT, b, Im, Jm, Km, Fm, 'between blocks' )

    ! A wall ghost cell in the stencil: one metric per group of quantities, with the positions of the group
    if ( .not. Connection_Mixed(domain, xc, r) ) return
    nfill(dir) = nfill(dir) + 1
    domain % blk(b) % dir(dir) % ig(fi,fj,fk) = nfill(dir)
    do g = 1, NG
      JTg = JT
      do q = 1, 2
        t = tb(q)
        if ( n(t) == 1 ) cycle
        dloc = Pos(xc(b), Im+e(1,t), Jm+e(2,t), Km+e(3,t), g) - Pos(xc(b), Im-e(1,t), Jm-e(2,t), Km-e(3,t), g)
        dgh  = Pos(xc(bs), sc(1,2*q), sc(2,2*q), sc(3,2*q), g) - Pos(xc(bs), sc(1,2*q-1), sc(2,2*q-1), sc(3,2*q-1), g)
        if ( mod(Fm,2) == 0 ) then
          JTg(:,t) = 0.25d0 * ( dloc + dgh )
        else
          JTg(:,t) = 0.25d0 * ( dgh + dloc )
        end if
      enddo
      domain % blk(b) % dir(dir) % Mg(:,:,g,nfill(dir)) = Face_Inverse ( JTg, b, Im, Jm, Km, Fm, 'between blocks' )
    enddo

  end subroutine Connection_Face_Metric


  !> Indices in the source block of the cells of Pg(:,3:6) of a connection record, as in Fill_BC_Ghost_Connection.
  subroutine Source_Neighbours ( rec, sc )
    use ARES_Advanced_Types_m, only: ARES_bc_type
    implicit none
    type(ARES_bc_type), intent(in) :: rec
    integer, intent(out) :: sc(3,4)
    ! Local
    integer :: Is, Js, Ks

    Is = rec % is ; Js = rec % js ; Ks = rec % ks
    select case ( rec % fs )
      case (1:2)
        sc(:,1) = [Is, Js - rec % d11, Ks - rec % d12]
        sc(:,2) = [Is, Js + rec % d11, Ks + rec % d12]
        sc(:,3) = [Is, Js - rec % d21, Ks - rec % d22]
        sc(:,4) = [Is, Js + rec % d21, Ks + rec % d22]
      case (3:4)
        sc(:,1) = [Is - rec % d11, Js, Ks - rec % d12]
        sc(:,2) = [Is + rec % d11, Js, Ks + rec % d12]
        sc(:,3) = [Is - rec % d21, Js, Ks - rec % d22]
        sc(:,4) = [Is + rec % d21, Js, Ks + rec % d22]
      case default
        sc(:,1) = [Is - rec % d11, Js - rec % d12, Ks]
        sc(:,2) = [Is + rec % d11, Js + rec % d12, Ks]
        sc(:,3) = [Is - rec % d21, Js - rec % d22, Ks]
        sc(:,4) = [Is + rec % d21, Js + rec % d22, Ks]
    end select

  end subroutine Source_Neighbours


  !> Face metric of a symmetry face (300) of block b: the cell and its mirror image across the face, D_a = (mirror - cell)
  !> along the direction of growing index, and the tangential differences of the cell and of the mirror images of its
  !> neighbours (BC_Symmetry_Visc). Faces across an index direction with one cell (the k faces of a 2-D mesh) keep the cell
  !> metric: nothing to store. A face of zero area (on an axis), whose viscous flux is zero, keeps the cell metric too: its
  !> mirror image is the cell itself, and the difference across it would be zero.
  subroutine Symmetry_Face_Metric ( blk, xc, rec )
    use ARES_Advanced_Types_m
    use ARES_Lib_Diffusive, only: Inverse3
    use ARES_Lib_BC_Fluxes, only: Face_Index
    implicit none
    type(ARES_block_type), intent(inout) :: blk
    type(centres_t), intent(in) :: xc
    type(ARES_bc_type), intent(in) :: rec
    ! Local
    integer :: Im, Jm, Km, Fm, dir, fi, fj, fk, t, n(3), modfm2
    integer :: e(3,3)
    real(R8) :: JT(3,3), xf(3), nf(3), x1(3), xp(3), xm(3)

    Im = rec % i ; Jm = rec % j ; Km = rec % k ; Fm = rec % f
    n = blk % dim
    e = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    call Face_Index ( Fm, dir, Im, Jm, Km, fi, fj, fk )
    if ( n(dir) == 1 ) return
    if ( .not. ( blk % dir(dir) % f(fi,fj,fk) % a > 0d0 ) ) then
      blk % dir(dir) % Mf(fi,fj,fk) % c = blk % M(Im,Jm,Km) % c
      return
    end if
    modfm2 = 1 - 2 * mod(Fm,2)
    call Face_Centre_Normal ( blk % node, Im, Jm, Km, Fm, xf, nf )
    x1 = xc % x(:,Im,Jm,Km)
    JT = Inverse3 ( blk % M(Im,Jm,Km) % c )
    JT(:,dir) = ( Mirror ( x1, xf, nf ) - x1 ) * modfm2
    do t = 1, 3
      if ( t == dir .or. n(t) == 1 ) cycle
      xp = Pos ( xc, Im+e(1,t), Jm+e(2,t), Km+e(3,t), 3 )        ! the velocity (BC_Symmetry_Visc takes no other gradient)
      xm = Pos ( xc, Im-e(1,t), Jm-e(2,t), Km-e(3,t), 3 )
      JT(:,t) = 0.25d0 * ( ( xp - xm ) + ( Mirror ( xp, xf, nf ) - Mirror ( xm, xf, nf ) ) )
    enddo
    blk % dir(dir) % Mf(fi,fj,fk) % c = Face_Inverse ( JT, rec % b, Im, Jm, Km, Fm, 'symmetry' )

  end subroutine Symmetry_Face_Metric


  !> Face metric of a face of type 200 of record r whose viscous flux BC_Rotational_Periodic_Visc computes (Rotational_Flux_Face:
  !> the k faces of a wedge, face 6 of a 3-D sector), stored on the face of the block as for the faces between blocks. The
  !> gradient takes the cell and its image across the face, the ghost cell, whose values are those of its source (the cell
  !> next to the opposite k face of the block, the cell itself on a wedge) turned by the angle of the periodicity (no turn on
  !> a translated sector), and along the face the central differences of the row of the cell and of the row of the image,
  !> which are those of the row of the source turned in the same way. The metric takes the same differences of the
  !> positions: across the face the image, at the position of Ghost_Position, minus the cell along growing k; along the face
  !> the differences of the row of the cell and those of the row of the source turned. So the face gradient is exact for the
  !> linear fields that the periodicity allows, u(R x) = R u(x) for a turn R about x, a velocity periodic along the period of a
  !> translated sector (docs/theory/numerics.md), on any mesh. An index direction with one cell keeps the column of the mean
  !> edges of the cell and of its image (Inverse3 of the inverse-mean metric of the record, bc%rM). With a wall ghost cell
  !> among the neighbours along the face, in the row of the cell or of the source, one metric per group of quantities, as at
  !> the interior faces.
  subroutine Rotational_Face_Metric ( domain, xc, r, nfill )
    use ARES_Advanced_Types_m
    use ARES_Parameters_m, only: guide
    use ARES_Lib_Metrics, only: delthe
    use ARES_Lib_Diffusive, only: Inverse3
    use ARES_Lib_BC_Fluxes, only: Face_Index
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    type(centres_t), intent(in) :: xc(:)
    integer, intent(in) :: r
    integer, intent(inout) :: nfill(3)
    ! Local
    integer :: b, Im, Jm, Km, Fm, Kg, Ks, dir, fi, fj, fk, t, g, n(3), modfm2
    integer :: e(3,3)
    real(R8) :: JT(3,3), JTg(3,3), Rt(3,3), ang

    if ( .not. Rotational_Flux_Face(domain, r) ) return
    b  = domain % bc(r) % b
    Im = domain % bc(r) % i ; Jm = domain % bc(r) % j ; Km = domain % bc(r) % k
    Fm = domain % bc(r) % f
    n  = domain % blk(b) % dim
    e  = reshape ( [1,0,0, 0,1,0, 0,0,1], [3,3] )
    Kg = Km - guide(Fm,3)
    Ks = n(3) - Km + 1
    modfm2 = 1 - 2 * mod(Fm,2)                   ! +1 on face 6 (the cell on the low side of the face), -1 on face 5
    call Face_Index ( Fm, dir, Im, Jm, Km, fi, fj, fk )
    ! The turn of the periodicity behind face Fm, as in Ghost_Position: +delthe behind face 6, -delthe behind face 5
    ang = delthe
    if ( Fm == 5 ) ang = -delthe
    Rt = reshape ( [1d0, 0d0, 0d0, 0d0, cos(ang), sin(ang), 0d0, -sin(ang), cos(ang)], [3,3] )

    JT = Inverse3 ( domain % bc(r) % rM )
    JT(:,3) = ( xc(b) % x(:,Im,Jm,Kg) - xc(b) % x(:,Im,Jm,Km) ) * modfm2
    do t = 1, 2
      if ( n(t) == 1 ) cycle
      JT(:,t) = 0.25d0 * ( ( xc(b) % x(:,Im+e(1,t),Jm+e(2,t),Km) - xc(b) % x(:,Im-e(1,t),Jm-e(2,t),Km) ) &
                         + matmul ( Rt, xc(b) % x(:,Im+e(1,t),Jm+e(2,t),Ks) - xc(b) % x(:,Im-e(1,t),Jm-e(2,t),Ks) ) )
    enddo
    domain % blk(b) % dir(dir) % Mf(fi,fj,fk) % c = Face_Inverse ( JT, b, Im, Jm, Km, Fm, 'type 200' )

    ! A wall ghost cell in the stencil: one metric per group of quantities, with the positions of the group
    if ( .not. Rotational_Mixed(domain, xc, r) ) return
    nfill(dir) = nfill(dir) + 1
    domain % blk(b) % dir(dir) % ig(fi,fj,fk) = nfill(dir)
    do g = 1, NG
      JTg = JT
      do t = 1, 2
        if ( n(t) == 1 ) cycle
        JTg(:,t) = 0.25d0 * ( ( Pos(xc(b), Im+e(1,t), Jm+e(2,t), Km, g) - Pos(xc(b), Im-e(1,t), Jm-e(2,t), Km, g) ) &
                            + matmul ( Rt, Pos(xc(b), Im+e(1,t), Jm+e(2,t), Ks, g) - Pos(xc(b), Im-e(1,t), Jm-e(2,t), Ks, g) ) )
      enddo
      domain % blk(b) % dir(dir) % Mg(:,:,g,nfill(dir)) = Face_Inverse ( JTg, b, Im, Jm, Km, Fm, 'type 200' )
    enddo

  end subroutine Rotational_Face_Metric


  !> Wall faces with diffusive-metric = stencil (types 301 and 302; not the coupled walls 103 and 104, whose wall flux
  !> comes from the other solver, so that no wall routine takes their gradient): the face metric of the wall gradient and
  !> the tangential differences of the imposed wall temperature, stored in the record (bc%wstencil, wM, wdT) and used by
  !> BC_Wall_Temperature and BC_Wall_Heat. The difference across the face is the one of the wall routines,
  !> modfm3 (phi_c - phi_w), with the position difference modfm3 (x_c - x_f); along the wall the differences of the wall
  !> values at the centres (mean of the four nodes) of the neighbouring faces of the same block face and of the same type,
  !> central, or one-sided where the wall ends; along an index direction with one cell, or with no such neighbour, the
  !> edge of the cell and no difference of the field. Exact for a linear field whose wall values are its values at the
  !> face centres: T_w of the 302 records, and the velocity, zero on the wall.
  subroutine Wall_Face_Stencils ( domain, xc )
    use ARES_Advanced_Types_m
    use ARES_Lib_Diffusive, only: Inverse3
    use ARES_Mod_MPI, only: is_local_block, mpi_abort_all
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    type(centres_t), intent(in) :: xc(:)
    ! Local
    integer, allocatable :: map(:,:)
    integer :: b, f, ii, r, rn(2), q, s, dir, tb(2), n(3), c(3), cn(3), typ, nside
    real(R8) :: JT(3,3), xf(3), nf(3), xfn(3,2), Ac(3,3), dT

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n = domain % blk(b) % dim
      do f = 1, 6
        if ( domain % local_bf_count(b,f) == 0 ) cycle
        dir = (f + 1) / 2
        select case (dir)
          case (1) ; tb = [2, 3]
          case (2) ; tb = [1, 3]
          case default ; tb = [1, 2]
        end select
        ! Records of the walls of this block face, by their two tangential cell indices
        allocate ( map(n(tb(1)), n(tb(2))) )
        map = 0
        do ii = domain % local_bf_first(b,f), domain % local_bf_first(b,f) + domain % local_bf_count(b,f) - 1
          r = domain % local_bf_idx(ii)
          if ( .not. Is_Wall(domain % bc(r) % type) ) cycle
          c = [domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k]
          map(c(tb(1)), c(tb(2))) = r
        enddo
        do ii = domain % local_bf_first(b,f), domain % local_bf_first(b,f) + domain % local_bf_count(b,f) - 1
          r = domain % local_bf_idx(ii)
          typ = domain % bc(r) % type
          if ( .not. Is_Wall(typ) ) cycle
          c = [domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k]
          call Face_Centre_Normal ( domain % blk(b) % node, c(1), c(2), c(3), f, xf, nf )
          if ( .not. ( norm2(nf) > 0d0 ) ) then
            ! a wall on a face of zero area (an axis): along the axis the centres of its neighbours coincide
            write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] diffusive-metric = stencil: block ', b, ', cell', c, ', face ', &
                                                  f, ' (wall): a wall on a face of zero area (an axis).'
            write(*,'(A)') ' [ERROR]   Its gradient is not defined there: use the symmetry condition (300) on the axis.'
            call mpi_abort_all ( 'wall on a face of zero area (see above)' )
          end if
          Ac = Inverse3 ( domain % blk(b) % M(c(1),c(2),c(3)) % c )        ! columns: the edges of the cell
          JT = 0d0
          JT(:,dir) = Modfm3_Of(f) * ( xc(b) % x(:,c(1),c(2),c(3)) - xf )
          domain % bc(r) % wdT = 0d0
          do q = 1, 2
            ! neighbouring wall faces along tb(q), same block face and same type: rn(1) at +1, rn(2) at -1
            rn = 0
            do s = 1, 2
              cn = c
              cn(tb(q)) = c(tb(q)) + 3 - 2*s
              if ( cn(tb(q)) < 1 .or. cn(tb(q)) > n(tb(q)) ) cycle
              rn(s) = map(cn(tb(1)), cn(tb(2)))
              if ( rn(s) /= 0 ) then
                if ( domain % bc(rn(s)) % type /= typ ) rn(s) = 0
              end if
              if ( rn(s) /= 0 ) then
                call Face_Centre_Normal ( domain % blk(b) % node, cn(1), cn(2), cn(3), f, xfn(:,s), nf )
              end if
            enddo
            nside = count ( rn /= 0 )
            if ( n(tb(q)) == 1 .or. nside == 0 ) then
              JT(:,tb(q)) = Ac(:,tb(q))                                      ! the edge, no difference
              cycle
            end if
            dT = 0d0
            if ( nside == 2 ) then
              JT(:,tb(q)) = 0.5d0 * ( xfn(:,1) - xfn(:,2) )
              if ( typ == 302 ) dT = 0.5d0 * ( domain % bc(rn(1)) % Tw - domain % bc(rn(2)) % Tw )
            elseif ( rn(1) /= 0 ) then
              JT(:,tb(q)) = xfn(:,1) - xf
              if ( typ == 302 ) dT = domain % bc(rn(1)) % Tw - domain % bc(r) % Tw
            else
              JT(:,tb(q)) = xf - xfn(:,2)
              if ( typ == 302 ) dT = domain % bc(r) % Tw - domain % bc(rn(2)) % Tw
            end if
            domain % bc(r) % wdT(q) = dT
          enddo
          domain % bc(r) % wM = Face_Inverse ( JT, b, c(1), c(2), c(3), f, 'wall' )
          domain % bc(r) % wstencil = .true.
        enddo
        deallocate ( map )
      enddo
    enddo

  contains

    logical function Is_Wall ( t )
      integer, intent(in) :: t
      Is_Wall = ( t == 301 .or. t == 302 )
    end function Is_Wall

    !> modfm3 of Compute_Modfm: +2 on the faces 1, 3, 5 (the wall at the low index of the cell), -2 on 2, 4, 6.
    integer function Modfm3_Of ( face )
      integer, intent(in) :: face
      Modfm3_Of = 2 * ( 2 * mod(face,2) - 1 )
    end function Modfm3_Of

  end subroutine Wall_Face_Stencils


  !> Cell metrics of the source terms of the turbulence models and of the eddy viscosity of the output, with
  !> diffusive-metric = stencil: M = (J_c^T)^-1, J_c^T(:,d) = (x_{c+e_d} - x_{c-e_d})/2, the same central differences of
  !> the positions that those terms take of the values (along an index direction with one cell, the edge of the cell and
  !> no difference). Exact for a linear field whose ghost cells sit where their value is right: the ghost cells of the
  !> fluxes, at the positions of Ghost_Position, behind a wall that of the group of the quantity (Ghost_Wall: 2 x_f - x_c
  !> for the quantities with a wall value, x_c + s (x_c - x_2) for the extrapolated ones). One metric per group: the
  !> velocity (blk%Ms); the density, whose ghost value follows the temperature behind a wall (blk%Msc: exact where the
  !> pressure does not vary across the wall cell, as in a boundary layer); the turbulence variables (blk%Mst, only on the
  !> blocks where their position is not that of the velocity, which no boundary type of ARES gives: allocated when one
  !> does). The velocity
  !> gradients of SA-RC (Compute_Velocity_Gradient): blk%Msa, behind the faces between blocks (101, 201) and of a sector
  !> turned about x (200) the position of the cell across, whose gradient the ghost cell holds (Connection_Ghost_Position,
  !> Ghost_Position), behind the other faces the cell itself, whose gradient is copied there.
  subroutine Cell_Source_Metrics ( domain, xc )
    use ARES_Advanced_Types_m
    use ARES_Global_m, only: gc, ndir
    use ARES_Parameters_m, only: guide
    use ARES_Lib_Metrics, only: delthe
    use ARES_Lib_Diffusive, only: Inverse3
    use ARES_Mod_MPI, only: is_local_block
    use ARES_Config_Types_m, only: obj_rans
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    type(centres_t), intent(in) :: xc(:)
    ! Local
    integer, parameter :: gs(3) = [3, 2, 4]                       ! groups of Ms, Msc, Mst
    integer :: b, i, j, k, d, s, m, n(3), c(3), cn(3), f, ii, r, g(3)
    real(R8) :: JT(3,3,3), JTa(3,3), Ac(3,3), xs(3), xg(3,2,3), xa(3,2)
    logical :: rc, turb, turned
    logical, allocatable :: across(:,:,:)

    rc = obj_rans % SpalartShur
    turned = ndir == 3 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 )   ! a sector turned about x
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      n = domain % blk(b) % dim
      turb = any ( btest(xc(b) % em, 2) .neqv. btest(xc(b) % em, 3) )
      if ( rc ) then
        ! ghost cells holding the velocity gradient of the cell across (Ghost_Velocity_Gradient)
        allocate ( across(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc) )
        across = .false.
        do f = 1, 6
          do ii = domain % local_bf_first(b,f), domain % local_bf_first(b,f) + domain % local_bf_count(b,f) - 1
            r = domain % local_bf_idx(ii)
            select case ( domain % bc(r) % type )
              case (101, 201)
              case (200)
                if ( .not. turned ) cycle
              case default
                cycle
            end select
            g = [ domain % bc(r) % i, domain % bc(r) % j, domain % bc(r) % k ] - guide(f,:)
            across(g(1),g(2),g(3)) = .true.
          enddo
        enddo
      end if
      if ( .not. allocated(domain % blk(b) % Ms) ) then
        allocate ( domain % blk(b) % Ms(lbound(domain % blk(b) % M,1):ubound(domain % blk(b) % M,1), &
                                        lbound(domain % blk(b) % M,2):ubound(domain % blk(b) % M,2), &
                                        lbound(domain % blk(b) % M,3):ubound(domain % blk(b) % M,3)) )
        allocate ( domain % blk(b) % Msc, mold = domain % blk(b) % Ms )
        if ( rc ) allocate ( domain % blk(b) % Msa, mold = domain % blk(b) % Ms )
        if ( turb ) allocate ( domain % blk(b) % Mst, mold = domain % blk(b) % Ms )
      end if
      domain % blk(b) % Ms = domain % blk(b) % M                          ! the ghost layers keep the cell metric
      domain % blk(b) % Msc = domain % blk(b) % M
      if ( rc ) domain % blk(b) % Msa = domain % blk(b) % M
      if ( turb ) domain % blk(b) % Mst = domain % blk(b) % M
      do k = 1, n(3) ; do j = 1, n(2) ; do i = 1, n(1)
        c = [i, j, k]
        xs = xc(b) % x(:,i,j,k)
        Ac = Inverse3 ( domain % blk(b) % M(i,j,k) % c )                   ! columns: the edges of the cell
        do d = 1, 3
          if ( n(d) == 1 ) then
            do m = 1, 3
              JT(:,d,m) = Ac(:,d)
            enddo
            JTa(:,d) = Ac(:,d)
            cycle
          end if
          do s = 1, 2                                                     ! s = 1: c + e_d, s = 2: c - e_d
            cn = c
            cn(d) = c(d) + 3 - 2*s
            do m = 1, 3
              xg(:,s,m) = Pos ( xc(b), cn(1), cn(2), cn(3), gs(m) )
            enddo
            xa(:,s) = xc(b) % x(:,cn(1),cn(2),cn(3))
            if ( rc .and. ( cn(d) < 1 .or. cn(d) > n(d) ) ) then
              if ( .not. across(cn(1),cn(2),cn(3)) ) xa(:,s) = xs
            end if
          enddo
          do m = 1, 3
            JT(:,d,m) = 0.5d0 * ( xg(:,1,m) - xg(:,2,m) )
          enddo
          JTa(:,d) = 0.5d0 * ( xa(:,1) - xa(:,2) )
        enddo
        domain % blk(b) % Ms(i,j,k) % c = Face_Inverse ( JT(:,:,1), b, i, j, k, 0, 'cell, source terms, velocity' )
        domain % blk(b) % Msc(i,j,k) % c = Face_Inverse ( JT(:,:,2), b, i, j, k, 0, 'cell, source terms, density' )
        if ( turb ) domain % blk(b) % Mst(i,j,k) % c = Face_Inverse ( JT(:,:,3), b, i, j, k, 0, &
                                                                      'cell, source terms, turbulence variables' )
        if ( rc ) domain % blk(b) % Msa(i,j,k) % c = Face_Inverse ( JTa, b, i, j, k, 0, 'cell, SA-RC' )
      enddo ; enddo ; enddo
      if ( allocated(across) ) deallocate ( across )
    enddo

  end subroutine Cell_Source_Metrics


  !> The cell metrics of the source terms of block blk: with diffusive-metric = stencil those of Cell_Source_Metrics (Mi
  !> for the velocity, Mc for the density, Mt for the turbulence variables, Ma for the velocity gradients of SA-RC, with
  !> their ghost cells), otherwise the metric of the cell for all of them. The one place that chooses them.
  subroutine Source_Metrics ( blk, Mi, Mc, Ma, Mt )
    use ARES_Advanced_Types_m, only: ARES_block_type
    use ARES_Base_Types_m, only: ARES_tensor_3D_type
    implicit none
    type(ARES_block_type), intent(in), target :: blk
    type(ARES_tensor_3D_type), pointer, intent(out) :: Mi(:,:,:)
    type(ARES_tensor_3D_type), pointer, intent(out), optional :: Mc(:,:,:), Ma(:,:,:), Mt(:,:,:)

    Mi => blk % M
    if ( allocated(blk % Ms) ) Mi => blk % Ms
    if ( present(Mc) ) then
      Mc => blk % M
      if ( allocated(blk % Msc) ) Mc => blk % Msc
    end if
    if ( present(Ma) ) then
      Ma => blk % M
      if ( allocated(blk % Msa) ) Ma => blk % Msa
    end if
    if ( present(Mt) ) then
      Mt => Mi
      if ( allocated(blk % Mst) ) Mt => blk % Mst
    end if

  end subroutine Source_Metrics

end module ARES_Lib_Stencil_Metrics
