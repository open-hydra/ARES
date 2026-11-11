!> Minimal 3D k-d tree for nearest-neighbor queries.
!> Built once from a point cloud, then queried many times (thread-safe for queries).
!> Used by Compute_Yn to replace the O(N_cells x N_walls) brute-force wall-distance
!> search with an O(N_cells x log N_walls) lookup.
!>
!> The query returns the same distance the brute-force scan did, bit for bit:
!> `sqrt` is monotonic and correctly rounded, so taking the minimum of the squared
!> distances and rooting once is exactly the minimum of the rooted distances.
!>
!> The points can also be the centres of objects of finite size (the wall faces
!> of Compute_Yn), each with a radius that bounds the distance from its centre to
!> any of its points: kdtree_nearest_object then returns the smallest distance to
!> the objects themselves, computed by a function of the caller only for the
!> objects that can still beat the best one found.
module ARES_Lib_KDTree
  use iso_fortran_env, only: R8 => real64, I4 => int32

  implicit none
  private

  public :: kdtree_t, kdtree_object_distance
  public :: kdtree_build, kdtree_nearest, kdtree_nearest_object, kdtree_free

  !> Opaque k-d tree. Tree structure is implicit: for the range [lo..hi], the
  !> root of the subtree is at position (lo+hi)/2, elements with smaller coord
  !> along the current axis live in [lo..mid-1], larger in [mid+1..hi]. Axis
  !> cycles x -> y -> z -> x -> ... with depth.
  type :: kdtree_t
    integer               :: n = 0
    real(R8), allocatable :: px(:), py(:), pz(:)
    integer,  allocatable :: perm(:)        ! perm(i) = original index of the point now at position i
    real(R8), allocatable :: pr(:)          ! radius of the object centred at position i (objects only)
    real(R8), allocatable :: rsub(:)        ! largest radius in the subtree rooted at position i (objects only)
  end type kdtree_t

  abstract interface
    !> Distance from (x,y,z) to the object of original index i.
    function kdtree_object_distance(i, x, y, z) result(d)
      import :: R8
      integer,  intent(in) :: i
      real(R8), intent(in) :: x, y, z
      real(R8)             :: d
    end function kdtree_object_distance
  end interface

contains

  !> Build the tree from (fx,fy,fz) of length n. O(N log N) average.
  !> `fr`, if present, is the radius of the object centred at each point: every
  !> point of object i lies within fr(i) of (fx(i),fy(i),fz(i)).
  subroutine kdtree_build(tree, fx, fy, fz, n, fr)
    type(kdtree_t), intent(out) :: tree
    integer,        intent(in)  :: n
    real(R8),       intent(in)  :: fx(n), fy(n), fz(n)
    real(R8), optional, intent(in) :: fr(n)
    integer :: i, j
    real(R8) :: r

    tree%n = n
    if (n <= 0) return

    allocate(tree%px(n), tree%py(n), tree%pz(n), tree%perm(n))
    do i = 1, n
      tree%px(i)   = fx(i)
      tree%py(i)   = fy(i)
      tree%pz(i)   = fz(i)
      tree%perm(i) = i
    end do

    ! Randomize input order so quickselect with last-element pivot avoids O(N^2)
    ! adversarial inputs (e.g. already-sorted wall centroids). The shuffle only
    ! changes the tree layout, never the distance a query returns, so ranks that
    ! happen to draw different sequences still agree exactly.
    do i = n, 2, -1
      call random_number(r)
      j = 1 + int(r * real(i, R8))
      if (j > i) j = i
      call swap_point(tree%px, tree%py, tree%pz, tree%perm, i, j)
    end do

    call build_rec(tree%px, tree%py, tree%pz, tree%perm, 1, n, 1)

    if (present(fr)) then
      allocate(tree%pr(n), tree%rsub(n))
      do i = 1, n
        tree%pr(i) = fr(tree%perm(i))
      end do
      call subtree_radius(tree%pr, tree%rsub, 1, n)
    end if
  end subroutine kdtree_build


  !> Free the tree storage.
  subroutine kdtree_free(tree)
    type(kdtree_t), intent(inout) :: tree
    if (allocated(tree%px))   deallocate(tree%px)
    if (allocated(tree%py))   deallocate(tree%py)
    if (allocated(tree%pz))   deallocate(tree%pz)
    if (allocated(tree%perm)) deallocate(tree%perm)
    if (allocated(tree%pr))   deallocate(tree%pr)
    if (allocated(tree%rsub)) deallocate(tree%rsub)
    tree%n = 0
  end subroutine kdtree_free


  !> Query: return the min Euclidean distance from (x,y,z) to any point in the tree.
  !> Read-only on the tree; safe to call concurrently from OpenMP threads.
  !> `idx`, if present, is the original index of the nearest point. Among
  !> equidistant points the lowest original index wins, so the answer does not
  !> depend on the shuffle and every rank picks the same point.
  subroutine kdtree_nearest(tree, x, y, z, d, idx)
    type(kdtree_t), intent(in)     :: tree
    real(R8),       intent(in)     :: x, y, z
    real(R8),       intent(out)    :: d
    integer, optional, intent(out) :: idx
    real(R8) :: best2
    integer  :: best_i

    if (tree%n <= 0) then
      d = huge(1.0_R8)
      if (present(idx)) idx = 0
      return
    end if
    best2  = huge(1.0_R8)
    best_i = huge(1)
    call nearest_rec(tree%px, tree%py, tree%pz, tree%perm, 1, tree%n, 1, x, y, z, best2, best_i)
    d = sqrt(best2)
    if (present(idx)) idx = best_i
  end subroutine kdtree_nearest


  !> Query on a tree built with object radii: return the smallest dist_fn(i,x,y,z)
  !> over the objects, and in `idx` the original index of the object that gives it
  !> (the lowest index among equal distances, so the answer does not depend on the
  !> shuffle or on the visiting order). dist_fn is called only for the objects whose
  !> lower bound |q - c_i| - r_i does not exceed the best distance found so far, and
  !> a subtree is skipped when its splitting plane is farther than that distance plus
  !> the largest radius in it. Read-only on the tree; safe to call concurrently from
  !> OpenMP threads if dist_fn is.
  !> With `resume` true the search starts from the d and idx passed in (a distance
  !> and object already found, for instance for another image of the query point)
  !> instead of from no object, and returns the better of the two with the same rule.
  subroutine kdtree_nearest_object(tree, x, y, z, dist_fn, d, idx, resume)
    type(kdtree_t), intent(in)    :: tree
    real(R8),       intent(in)    :: x, y, z
    procedure(kdtree_object_distance) :: dist_fn
    real(R8),       intent(inout) :: d
    integer,        intent(inout) :: idx
    logical, optional, intent(in) :: resume
    logical :: cont

    cont = .false.
    if (present(resume)) cont = resume
    if (.not. cont) then
      d   = huge(1.0_R8)
      idx = huge(1)
    end if
    if (tree%n <= 0 .or. .not. allocated(tree%pr)) then
      if (.not. cont) idx = 0
      return
    end if
    call nearest_object_rec(tree, 1, tree%n, 1, x, y, z, dist_fn, d, idx)
  end subroutine kdtree_nearest_object


  ! -------------------------------------------------------------------------
  ! Internal helpers
  ! -------------------------------------------------------------------------

  recursive subroutine build_rec(px, py, pz, perm, lo, hi, axis)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: lo, hi, axis
    integer :: mid, next_axis

    if (lo >= hi) return
    mid = (lo + hi) / 2
    call nth_element(px, py, pz, perm, lo, hi, mid, axis)
    next_axis = mod(axis, 3) + 1
    call build_rec(px, py, pz, perm, lo,    mid - 1, next_axis)
    call build_rec(px, py, pz, perm, mid+1, hi,      next_axis)
  end subroutine build_rec


  recursive subroutine nearest_rec(px, py, pz, perm, lo, hi, axis, qx, qy, qz, best2, best_i)
    real(R8), intent(in)    :: px(:), py(:), pz(:)
    integer,  intent(in)    :: perm(:)
    real(R8), intent(in)    :: qx, qy, qz
    real(R8), intent(inout) :: best2
    integer,  intent(inout) :: best_i
    integer,  intent(in)    :: lo, hi, axis
    integer :: mid, next_axis
    real(R8) :: d2, dx, dy, dz, diff

    if (lo > hi) return
    mid = (lo + hi) / 2

    dx = px(mid) - qx
    dy = py(mid) - qy
    dz = pz(mid) - qz
    d2 = dx*dx + dy*dy + dz*dz
    if (d2 < best2 .or. (d2 == best2 .and. perm(mid) < best_i)) then
      best2  = d2
      best_i = perm(mid)
    end if

    select case (axis)
      case (1) ; diff = qx - px(mid)
      case (2) ; diff = qy - py(mid)
      case (3) ; diff = qz - pz(mid)
    end select
    next_axis = mod(axis, 3) + 1

    if (diff <= 0.0_R8) then
      ! query on left side: descend left first
      call nearest_rec(px, py, pz, perm, lo,      mid - 1, next_axis, qx, qy, qz, best2, best_i)
      ! `<=`, not `<`: a point on the splitting plane can tie the best distance
      ! with a lower original index, and the tie-break must still see it.
      if (diff*diff <= best2) &
        call nearest_rec(px, py, pz, perm, mid+1, hi,      next_axis, qx, qy, qz, best2, best_i)
    else
      call nearest_rec(px, py, pz, perm, mid+1,   hi,      next_axis, qx, qy, qz, best2, best_i)
      if (diff*diff <= best2) &
        call nearest_rec(px, py, pz, perm, lo,    mid - 1, next_axis, qx, qy, qz, best2, best_i)
    end if
  end subroutine nearest_rec


  !> rsub(mid) = largest radius in the subtree rooted at mid = (lo+hi)/2.
  recursive subroutine subtree_radius(pr, rsub, lo, hi)
    real(R8), intent(in)    :: pr(:)
    real(R8), intent(inout) :: rsub(:)
    integer,  intent(in)    :: lo, hi
    integer :: mid

    if (lo > hi) return
    mid = (lo + hi) / 2
    rsub(mid) = pr(mid)
    call subtree_radius(pr, rsub, lo, mid - 1)
    call subtree_radius(pr, rsub, mid + 1, hi)
    if (lo <= mid - 1) rsub(mid) = max(rsub(mid), rsub((lo + mid - 1) / 2))
    if (mid + 1 <= hi) rsub(mid) = max(rsub(mid), rsub((mid + 1 + hi) / 2))
  end subroutine subtree_radius


  recursive subroutine nearest_object_rec(tree, lo, hi, axis, qx, qy, qz, dist_fn, best, best_i)
    type(kdtree_t), intent(in)    :: tree
    integer,        intent(in)    :: lo, hi, axis
    real(R8),       intent(in)    :: qx, qy, qz
    procedure(kdtree_object_distance) :: dist_fn
    real(R8),       intent(inout) :: best
    integer,        intent(inout) :: best_i
    integer :: mid, next_axis
    real(R8) :: d, diff

    if (lo > hi) return
    mid = (lo + hi) / 2

    ! The object at mid can reach the query no closer than its centre less its radius
    if (sqrt((tree%px(mid) - qx)**2 + (tree%py(mid) - qy)**2 + (tree%pz(mid) - qz)**2) - tree%pr(mid) <= best) then
      d = dist_fn(tree%perm(mid), qx, qy, qz)
      if (d < best .or. (d == best .and. tree%perm(mid) < best_i)) then
        best   = d
        best_i = tree%perm(mid)
      end if
    end if

    select case (axis)
      case (1) ; diff = qx - tree%px(mid)
      case (2) ; diff = qy - tree%py(mid)
      case (3) ; diff = qz - tree%pz(mid)
    end select
    next_axis = mod(axis, 3) + 1

    ! Every centre of the far subtree lies at least |diff| from the query along the
    ! axis, so its objects lie at least |diff| - rsub(its root) from it. `<=`, as in
    ! nearest_rec, lets an equal distance with a lower index be seen.
    if (diff <= 0.0_R8) then
      call nearest_object_rec(tree, lo, mid - 1, next_axis, qx, qy, qz, dist_fn, best, best_i)
      if (mid + 1 <= hi) then
        if (-diff - tree%rsub((mid + 1 + hi) / 2) <= best) &
          call nearest_object_rec(tree, mid + 1, hi, next_axis, qx, qy, qz, dist_fn, best, best_i)
      end if
    else
      call nearest_object_rec(tree, mid + 1, hi, next_axis, qx, qy, qz, dist_fn, best, best_i)
      if (lo <= mid - 1) then
        if (diff - tree%rsub((lo + mid - 1) / 2) <= best) &
          call nearest_object_rec(tree, lo, mid - 1, next_axis, qx, qy, qz, dist_fn, best, best_i)
      end if
    end if
  end subroutine nearest_object_rec


  !> Rearrange (px, py, pz, perm) in [lo..hi] so that the element at position k
  !> is the k-th smallest along the current axis, with smaller-or-equal to the
  !> left and greater-or-equal to the right. Quickselect, average O(hi-lo).
  recursive subroutine nth_element(px, py, pz, perm, lo, hi, k, axis)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: lo, hi, k, axis
    integer :: p

    if (lo >= hi) return
    p = partition(px, py, pz, perm, lo, hi, axis)
    if (k < p) then
      call nth_element(px, py, pz, perm, lo,  p - 1, k, axis)
    else if (k > p) then
      call nth_element(px, py, pz, perm, p+1, hi,    k, axis)
    end if
  end subroutine nth_element


  !> Lomuto partition around the pivot at position hi. Returns final pivot index.
  function partition(px, py, pz, perm, lo, hi, axis) result(p)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: lo, hi, axis
    integer :: p, i
    real(R8) :: pivot

    select case (axis)
      case (1) ; pivot = px(hi)
      case (2) ; pivot = py(hi)
      case (3) ; pivot = pz(hi)
    end select

    p = lo
    do i = lo, hi - 1
      select case (axis)
        case (1) ; if (px(i) < pivot) then
                     call swap_point(px, py, pz, perm, p, i); p = p + 1
                   end if
        case (2) ; if (py(i) < pivot) then
                     call swap_point(px, py, pz, perm, p, i); p = p + 1
                   end if
        case (3) ; if (pz(i) < pivot) then
                     call swap_point(px, py, pz, perm, p, i); p = p + 1
                   end if
      end select
    end do
    call swap_point(px, py, pz, perm, p, hi)
  end function partition


  subroutine swap_point(px, py, pz, perm, i, j)
    real(R8), intent(inout) :: px(:), py(:), pz(:)
    integer,  intent(inout) :: perm(:)
    integer,  intent(in)    :: i, j
    real(R8) :: t
    integer  :: it
    if (i == j) return
    t = px(i); px(i) = px(j); px(j) = t
    t = py(i); py(i) = py(j); py(j) = t
    t = pz(i); pz(i) = pz(j); pz(j) = t
    it = perm(i); perm(i) = perm(j); perm(j) = it
  end subroutine swap_point

end module ARES_Lib_KDTree
