module ARES_Lib_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64, IK8 => int64
  use ARES_Base_Types_m
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m

  implicit none
  real(R8), public :: delthe = 0d0  ! grid axisymmetric angle
  real(R8), public :: cos_delthe = 1d0, sin_delthe = 0d0  ! its cosine and sine, set with it (Check_Mesh_Type)

  ! Wall faces of Compute_Yn, read by the distance functions while the k-d tree is queried: for face n, wall_corner(:,1:4,n)
  ! are its corners in cyclic order and wall_corner(:,0,n) their mean; on an axisymmetric wedge (wall_meridional) the
  ! corners are in the meridional plane (x, r, 0), where corners 1 and 2 are the ends of the trace of the face. Across
  ! periodic boundaries each face is moved, whole, by a periodic translation (Wall_Images). wall_fast(n): the face takes
  ! the four triangles of Wall_Face_Distance, which are the face to the resolution of its coordinates (Face_Fast_Test);
  ! otherwise the certified minimum over the bilinear face (Bilinear_Face_Distance). wall_ghat(n): an upper bound of the gap
  ! between those triangles and the face; wall_fan_ok(n): the triangles are conditioned enough for the margin of the filter
  ! of the second pass (Wall_Distance_Blk); wall_radius(n): every point of the face lies within it of wall_corner(:,0,n).
  real(R8), allocatable :: wall_corner(:,:,:)
  logical               :: wall_meridional = .false.
  logical,  allocatable :: wall_fast(:), wall_fan_ok(:)
  real(R8), allocatable :: wall_ghat(:), wall_radius(:)
  ! Upper bound of the wall distance of the point being searched, from the first pass of Wall_Distance_Blk, and the least
  ! distance computed so far in the second, read by the filter of the second pass (Wall_Exact_Distance): one per thread
  real(R8)              :: wall_upper = 0d0, wall_least = 0d0
  ! The four triangles in local coordinates (Local_Fan_Distance) of the faces met in the first pass of a point, with the
  ! image of the wall they were met in (wall_image): the second pass of the same point takes them from this table instead of
  ! computing them again. A slot holds face, image, the point it belongs to (fanc_tag against fanc_now, which grows by one
  ! for each point, so that the table needs no clearing) and the value; a full table only means more computing. One table
  ! per thread
  integer, parameter    :: FANC = 1024
  integer               :: fanc_face(0:FANC-1) = 0, fanc_img(0:FANC-1) = 0, wall_image = 0
  integer(IK8)          :: fanc_tag(0:FANC-1) = 0, fanc_now = 0
  real(R8)              :: fanc_val(0:FANC-1) = 0d0
  !$omp threadprivate(wall_upper, wall_least, fanc_face, fanc_img, fanc_tag, fanc_now, fanc_val, wall_image)
  private :: wall_corner, wall_meridional, wall_fast, wall_fan_ok, wall_ghat, wall_radius, wall_upper, wall_least
  private :: FANC, fanc_face, fanc_img, fanc_tag, fanc_now, fanc_val, wall_image, Fan_Cached
  private :: Cell_Centre, Meridional, Wall_Face_Distance, Point_Segment_Distance, Point_Triangle_Distance, Cross3
  private :: Wall_Upper_Distance, Wall_Exact_Distance, Wall_Face_Exact, Local_Fan_Distance, Filter_Margin, Tri_Dist_Lb, Orient
  private :: Vnorm
  ! A degenerate or inverted cell found by Compute_Metric_Tensor, inside the OpenMP threads of Setup_Metrics: the run stops
  ! after the parallel region, on every rank (First_Bad_Metric, Stop_On_Bad_Metric)
  logical, public, protected :: bad_metric = .false.
  private :: First_Bad_Metric

contains

  !> Mesh type, from the nodes of every block.
  !>  - ndir is the number of index directions with more than one cell: 3 if a block has more than one cell in k,
  !>    otherwise 2 if a block has more than one cell in j, otherwise 1.
  !>  - delthe is the angle of the axisymmetric wedge of a 2-D mesh, or of the sectors closed by the rotational
  !>    periodicity (type 200 on faces 5 and 6) of a 3-D mesh, and 0 otherwise. A 2-D mesh is a wedge when, in every
  !>    block, the k = 1 node plane is the k = 0 plane turned about the x axis by the same angle; it is plane when no
  !>    block is so turned. In 3-D the k = dim(3) node plane of every block with a type-200 face is the k = 0 plane turned
  !>    about the x axis by the same angle, or translated (delthe = 0); a block closed on itself, a full turn about the x
  !>    axis with its last node plane on the first, is translated by zero. Every node is checked (Node_Planes); blocks
  !>    that disagree stop the run.
  !>  - The angle is the azimuth turned by the node (0,1) of the block, or by the node farthest from
  !>    the axis when that one lies on it. The angles of two blocks agree when they differ by at most ETA of the angle.
  !>  - The checks are relative to the move of the nodes per cell in k, measured along the k lines (K_Move_Per_Cell; the
  !>    thickness of a 2-D block), not to the size of the coordinates: the
  !>    coordinates of a mesh are written with a finite number of digits, and their round-off, relative to that move, is
  !>    what a turn or a translation leaves over. A check relative to the size of the coordinates, at 1e-6 of it, fails a
  !>    wedge of 5 degrees written with 5 significant digits, which is then taken for a plane mesh, and passes a plane slab
  !>    far from the x axis, which is then taken for a wedge. A 2-D block nearer a turn than a translation, or nearer a
  !>    translation than a turn, but not to within ETA of its thickness, stops the run: its nodes are written with too
  !>    few digits to tell its type, or its k planes are neither parallel nor turned about the x axis.
  !>  - cos_delthe and sin_delthe: its cosine and sine, computed here once for the fluxes and the ghost cells of the
  !>    type-200 faces, which turn vectors by delthe at every stage.
  subroutine Check_Mesh_Type( domain )
    use ARES_Mod_MPI, only: mpi_stop_all
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    real(R8), parameter :: TOL = 1d-6, ETA = 1d-2, PI = 3.141592653589793d0
    integer  :: b, l, nb, nturn, nmove
    real(R8) :: sr, ang, ang_b, fit, fit_move
    logical  :: turned, moved, nearer_turn, nearer_move
    logical, allocatable :: sector(:)
    character(len=512) :: msg   ! room for the longest message: 255 characters and the block number

    nb = domain % nb

    ndir = 1
    do b = 1, nb
      if ( domain % Blk(b) % dim(2) > 1 ) ndir = max ( ndir, 2 )
      if ( domain % Blk(b) % dim(3) > 1 ) ndir = 3
    enddo

    ! Size of the distance from the x axis, to tell a probe node on the axis
    sr = 0d0
    do b = 1, nb
      associate ( node => domain % Blk(b) % node )
        sr = max ( sr, maxval ( sqrt ( node(:,:,:) % c(2)**2 + node(:,:,:) % c(3)**2 ) ) )
      end associate
    enddo

    delthe = 0d0
    select case ( ndir )

      case ( 2 )
        ! Wedge if every block is turned by the same angle; plane if none is
        nturn = 0
        do b = 1, nb
          call Node_Planes ( domain % Blk(b), 1, sr, TOL, ETA, turned, ang_b, moved, fit, nearer_turn, fit_move, &
                             nearer_move )
          if ( nearer_turn .and. .not. turned ) then
            write(msg,'(A,I0,A,ES10.3,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' is nearer a wedge of angle ', ang_b, &
              ' rad than a plane block, but its nodes match the turn only to ', fit, ' of its thickness (1e-2 needed):'// &
              ' write the node coordinates with more significant digits'
            call mpi_stop_all ( trim(msg) )
          end if
          if ( nearer_move .and. .not. moved ) then
            write(msg,'(A,I0,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' is nearer a plane block than a wedge, but its'// &
              ' k = 1 node plane is the k = 0 plane translated only to ', fit_move, ' of its thickness (1e-2 needed):'// &
              ' its k planes must be parallel, or the coordinates written with more significant digits'
            call mpi_stop_all ( trim(msg) )
          end if
          if ( .not. turned ) cycle
          nturn = nturn + 1
          if ( nturn == 1 ) then
            ang = ang_b
          elseif ( abs ( ang_b - ang ) > ETA * abs ( ang ) ) then
            write(msg,'(A,I0,A,ES12.5,A,ES12.5,A)') '[ERROR] Mesh type: block ', b, ' is a wedge of angle ', ang_b, &
              ' rad, another block one of ', ang, ' rad'
            call mpi_stop_all ( trim(msg) )
          end if
        enddo
        if ( nturn == nb ) then
          delthe = ang
        elseif ( nturn > 0 ) then
          write(msg,'(A,I0,A,I0,A)') '[ERROR] Mesh type: ', nturn, ' of the ', nb, &
            ' blocks of the 2-D mesh are axisymmetric wedges (k = 1 plane turned about the x axis), the others are not'
          call mpi_stop_all ( trim(msg) )
        end if

      case ( 3 )
        ! Sectors closed by the rotational periodicity: the blocks with a type-200 face on faces 5-6
        allocate ( sector(nb) )
        sector = .false.
        do l = 1, domain % nbound
          if ( domain % bc(l) % type == 200 .and. domain % bc(l) % f >= 5 ) sector(domain % bc(l) % b) = .true.
        enddo
        nturn = 0
        nmove = 0
        do b = 1, nb
          if ( .not. sector(b) ) cycle
          call Node_Planes ( domain % Blk(b), domain % Blk(b) % dim(3), sr, TOL, ETA, turned, ang_b, moved, fit, &
                             nearer_turn, fit_move, nearer_move )
          if ( turned ) then
            nturn = nturn + 1
            if ( nturn == 1 ) then
              ang = ang_b
            elseif ( abs ( ang_b - ang ) > ETA * abs ( ang ) ) then
              write(msg,'(A,I0,A,ES12.5,A,ES12.5,A)') '[ERROR] Mesh type: block ', b, ' is a sector of angle ', ang_b, &
                ' rad, another block one of ', ang, ' rad'
              call mpi_stop_all ( trim(msg) )
            end if
          elseif ( moved ) then
            nmove = nmove + 1
          else
            write(msg,'(A,I0,A,ES9.2,A,ES9.2,A)') '[ERROR] Mesh type: block ', b, ' has type 200 on its k faces, but its'// &
              ' last k node plane is neither the first one turned about the x axis nor translated, to 1e-2 of the move of'// &
              ' its nodes per cell in k (the turn fits to ', fit, ', the translation to ', fit_move, ')'
            call mpi_stop_all ( trim(msg) )
          end if
        enddo
        if ( nturn > 0 .and. nmove > 0 ) then
          call mpi_stop_all ( '[ERROR] Mesh type: type 200 closes sectors turned about the x axis in some blocks and'// &
                               ' translated ones in others' )
        end if
        if ( nturn > 0 ) delthe = ang
        deallocate ( sector )

    end select

    cos_delthe = cos ( delthe )
    sin_delthe = sin ( delthe )

  contains

    !> Compare the node planes k = 0 and k = ktop of a block. cell is the move of the nodes per cell in k measured along
    !> the k lines (K_Move_Per_Cell; in 2-D, ktop = 1, the thickness of the block). turned: every node of the plane ktop
    !> is the node of the plane 0 turned about the x axis by ang, with |ang| > tol, to within eta of cell, and closer
    !> than as a translation. moved: every node is the node of the plane 0 translated by the move of the node (0,0), to
    !> within eta of cell, and closer than as a turn, or with |ang| <= tol, which is no turn: the turn by a null angle is
    !> the translation by zero, so a block closed on itself (a full turn about the x axis, its plane ktop on the plane 0)
    !> is translated by zero, the periodicity by the identity. Comparing the two keeps a wedge from passing for a
    !> translation and a slab far from the axis, which a small turn nearly fits, from passing for a wedge; taking cell,
    !> not the thickness, keeps a periodic plane of a sector of many cells that does not match by a fraction of a cell
    !> from passing for a turn. The move along the k lines, not the chord from the plane 0 to the plane ktop, makes
    !> cell the move per cell of a sector past half a turn too (its arcs in the largest-component norm, between 0.71 and
    !> 1 of the arc; the chord understates it, 4.2 times at 270 degrees), and of a block closed on itself, where the chord
    !> is zero or round-off. In 2-D a block whose probe node keeps its azimuth while its plane k = 1 is not the plane
    !> k = 0 translated to eta, and which the translation does not fit better than the turn by the null angle (planes
    !> that touch or cross, or a plane k = 1 scaled radially about the axis), stops here: without the null angle it would
    !> pass for plane (with crossed planes it would stop at its folded cells, Check_Folded_Cells). fit, fit_move: the
    !> misfits of the turn and of the translation over cell; nearer_turn: |ang| > tol and the turn fits better than the
    !> translation; nearer_move: the translation fits better than the turn, or |ang| <= tol.
    subroutine Node_Planes ( blk, ktop, sr, tol, eta, turned, ang, moved, fit, nearer_turn, fit_move, nearer_move )
      type(ARES_block_type), intent(in) :: blk
      integer, intent(in)               :: ktop
      real(R8), intent(in)              :: sr, tol, eta
      logical, intent(out)              :: turned, moved, nearer_turn, nearer_move
      real(R8), intent(out)             :: ang, fit, fit_move
      ! Local
      integer  :: i, j, i0, j0
      real(R8) :: theta1, theta2, ca, sa, rmax, r, a(3), c(3), shift(3), turn_off, move_off, cell

      associate ( node => blk % node )
        ! Probe node: (0,1), or the node of the plane 0 farthest from the axis when (0,1) lies on it
        i0 = 0
        j0 = 1
        if ( sqrt ( node(0,1,0)%c(2)**2 + node(0,1,0)%c(3)**2 ) <= tol * sr ) then
          rmax = -1d0
          do j = 0, blk % dim(2)
            do i = 0, blk % dim(1)
              r = sqrt ( node(i,j,0)%c(2)**2 + node(i,j,0)%c(3)**2 )
              if ( r > rmax ) then
                rmax = r
                i0 = i
                j0 = j
              end if
            enddo
          enddo
        end if
        theta1 = atan2 ( node(i0,j0,0)%c(3), node(i0,j0,0)%c(2) )
        theta2 = atan2 ( node(i0,j0,ktop)%c(3), node(i0,j0,ktop)%c(2) )
        ang = theta2 - theta1
        if ( ang >  pi ) ang = ang - 2d0 * pi
        if ( ang < -pi ) ang = ang + 2d0 * pi

        ca = cos ( ang )
        sa = sin ( ang )
        shift = node(0,0,ktop)%c - node(0,0,0)%c
        turn_off = 0d0
        move_off = 0d0
        do j = 0, blk % dim(2)
          do i = 0, blk % dim(1)
            a = node(i,j,0)%c
            c = node(i,j,ktop)%c
            turn_off = max ( turn_off, abs ( c(1) - a(1) ), abs ( c(2) - ( a(2) * ca - a(3) * sa ) ), &
                                       abs ( c(3) - ( a(2) * sa + a(3) * ca ) ) )
            move_off = max ( move_off, maxval ( abs ( c - a - shift ) ) )
          enddo
        enddo
        cell = K_Move_Per_Cell ( blk, ktop )
        nearer_turn = abs ( ang ) > tol .and. turn_off < move_off
        nearer_move = move_off < turn_off .or. abs ( ang ) <= tol
        turned = nearer_turn .and. turn_off <= eta * cell
        moved  = nearer_move .and. move_off <= eta * cell
        fit      = huge ( 1d0 )
        fit_move = huge ( 1d0 )
        if ( cell > 0d0 ) fit = turn_off / cell
        if ( cell > 0d0 ) fit_move = move_off / cell
      end associate

    end subroutine Node_Planes

  end subroutine Check_Mesh_Type


  !> The move of the nodes of a block per cell in k, from the node plane 0 to the plane ktop, measured along the k
  !> lines: for every node of the plane 0 the sum, over the cells of its k line, of its move across each cell (the
  !> largest component), and the largest of these sums over ktop. Along straight k lines it is the move from the plane 0
  !> to the plane ktop per cell (to the round-off of the sum); with one cell in k (ktop = 1) it is that move, bit for
  !> bit; along the arcs of a sector it is the length of the arc per cell in the largest-component norm, between 0.71
  !> and 1 of the arc (0.90 of it over a full turn), which the chord from plane to plane understates past half a turn
  !> (4.2 times at 270 degrees); for a block closed on itself (a full turn about the x axis) it is that length per cell,
  !> where the chord is zero or round-off. The scale of the checks of the periodic planes (Check_Mesh_Type) and of the
  !> periodic translations of the blocks closed by type 200 (Wall_Images).
  pure function K_Move_Per_Cell ( blk, ktop ) result ( cell )
    type(ARES_block_type), intent(in) :: blk
    integer, intent(in)               :: ktop
    real(R8)                          :: cell
    ! Local
    integer  :: i, j, k
    real(R8) :: s

    cell = 0d0
    do j = 0, blk % dim(2)
      do i = 0, blk % dim(1)
        s = 0d0
        do k = 1, ktop
          s = s + maxval ( abs ( blk % node(i,j,k) % c - blk % node(i,j,k-1) % c ) )
        enddo
        cell = max ( cell, s )
      enddo
    enddo
    cell = cell / ktop

  end function K_Move_Per_Cell


  subroutine Compute_Norm_Area( b )
    implicit none
    class(ARES_block_type) :: b  
    ! Local
    real(R8) :: Ai, snix, sniy, sniz, Aj, snjx, snjy, snjz, Ak, snkx, snky, snkz
    integer      :: i, j, k, im, jm, km
    real(R8) :: d1(3), d2(3), d3(3)
    real(R8) :: snixx, sniyy, snizz, snjxx, snjyy, snjzz, snkxx, snkyy, snkzz
    real(R8) :: scal, signi, signj, signk
  
    ! Peliminary operations
    im = b%dim(1)
    jm = b%dim(2)
    km = b%dim(3)
  
    ! compute sign of normal vectors to the intefaces
  
    ! i direction
    d1 = b%node(1,1,1)%c - b%node(1,0,0)%c
    d2 = b%node(1,0,1)%c - b%node(1,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snixx=d3(1)
    sniyy=d3(2)
    snizz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(1,0,1)%c + b%node(1,0,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(0,1,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snixx+d3(2)*sniyy+d3(3)*snizz
  
    signi=sign(1.d0,scal)

    ! j direction
    d1 = b%node(0,1,1)%c - b%node(1,1,0)%c
    d2 = b%node(1,1,1)%c - b%node(0,1,0)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snjxx=d3(1)
    snjyy=d3(2)
    snjzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,1,0)%c + b%node(1,1,0)%c )
    d2 = 0.25d0*( b%node(1,0,1)%c + b%node(0,0,1)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snjxx+d3(2)*snjyy+d3(3)*snjzz
  
    signj=sign(1.d0,scal)

    ! k direction
    d1 = b%node(0,1,1)%c - b%node(1,0,1)%c
    d2 = b%node(0,0,1)%c - b%node(1,1,1)%c
  
    d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))
    d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))
    d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))
  
    snkxx=d3(1)
    snkyy=d3(2)
    snkzz=d3(3)
  
    d1 = 0.25d0*( b%node(1,1,1)%c + b%node(0,1,1)%c + b%node(0,0,1)%c + b%node(1,0,1)%c )
    d2 = 0.25d0*( b%node(1,1,0)%c + b%node(0,1,0)%c + b%node(0,0,0)%c + b%node(1,0,0)%c )
  
    d3(1)=d1(1)-d2(1)
    d3(2)=d1(2)-d2(2)
    d3(3)=d1(3)-d2(3)
  
    scal=d3(1)*snkxx+d3(2)*snkyy+d3(3)*snkzz
  
    signk=sign(1.d0,scal)
    
    ! compute metrics: n, A
    
    !$omp parallel private (d1,d2,d3,i,j,k,snix,sniy,sniz,Ai,Aj,snjx,snjy,snjz,Ak,snkx,snky,snkz)
    
    ! i direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 1, jm ; do i = 0, im
  
      d1 = b%node(i,j,k)%c - b%node(i,j-1,k-1)%c
      d2 = b%node(i,j-1,k)%c - b%node(i,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ai = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Ai == 0d0) then
        snix = 0d0
        sniy = 0d0
        sniz = 0d0
      else
        snix = d3(1)/Ai*signi
        sniy = d3(2)/Ai*signi
        sniz = d3(3)/Ai*signi
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(1)%f(i,j,k)%a = Ai
      b%dir(1)%f(i,j,k)%n = [ snix, sniy, sniz ]
      
    enddo; enddo; enddo
  
    ! j direction
    !$omp do collapse(3)
    do k = 1, km ; do j = 0, jm ; do  i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j,k-1)%c
      d2 = b%node(i,j,k)%c - b%node(i-1,j,k-1)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Aj = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Aj == 0d0) then
        snjx = 0d0
        snjy = 0d0
        snjz = 0d0
      else
        snjx = d3(1)/Aj*signj
        snjy = d3(2)/Aj*signj
        snjz = d3(3)/Aj*signj
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(2)%f(i,j,k)%A = Aj
      b%dir(2)%f(i,j,k)%n = [ snjx, snjy, snjz ]
  
    enddo; enddo; enddo
  
    ! k direction
    !$omp do collapse(3)
    do k = 0, km ; do j = 1, jm ; do i = 1, im
  
      d1 = b%node(i-1,j,k)%c - b%node(i,j-1,k)%c
      d2 = b%node(i-1,j-1,k)%c - b%node(i,j,k)%c
  
      d3(1)=(d1(2)*d2(3)-d1(3)*d2(2))*.5d0
      d3(2)=(d1(3)*d2(1)-d1(1)*d2(3))*.5d0
      d3(3)=(d1(1)*d2(2)-d1(2)*d2(1))*.5d0
  
      Ak = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
  
      if (Ak == 0d0) then
        snkx = 0d0
        snky = 0d0
        snkz = 0d0
      else
        snkx = d3(1)/Ak*signk
        snky = d3(2)/Ak*signk
        snkz = d3(3)/Ak*signk
      end if
  
      !% Assign computed normal and area to metrics object
      b%dir(3)%f(i,j,k)%A = Ak
      b%dir(3)%f(i,j,k)%n = [ snkx, snky, snkz ]
  
    enddo; enddo; enddo
    !$omp end parallel
  
  end subroutine Compute_Norm_Area


  !> True for the first degenerate or inverted cell of this rank, in any OpenMP thread: the flag is shared and set in a
  !> critical section, so that one thread writes the cell. A STOP there would end the process while the other threads go
  !> on (with two threads the exit code can be 0, the message of the cell lost or written twice).
  logical function First_Bad_Metric ( ) result ( first )
    !$omp critical (ares_bad_metric)
    first = .not. bad_metric
    bad_metric = .true.
    !$omp end critical (ares_bad_metric)
  end function First_Bad_Metric


  !> Every rank stops, exit code 1, when one of them found a degenerate or inverted cell (bad_metric). Each rank computes
  !> the cells of its own blocks, so the ranks take the decision together (the largest of a flag) before mpi_stop_all, as
  !> Check_IC does. Setup_Metrics calls it on every rank, outside its parallel regions.
  subroutine Stop_On_Bad_Metric ( )
    use ARES_Mod_MPI, only: mpi_allreduce_max_r8, mpi_stop_all
    implicit none
    real(R8) :: flag, flag_all

    flag = merge ( 1d0, 0d0, bad_metric )
    call mpi_allreduce_max_r8 ( flag, flag_all )
    if ( flag_all > 0d0 ) call mpi_stop_all ( 'degenerate or inverted cell (see above)' )
  end subroutine Stop_On_Bad_Metric


  subroutine Compute_Metric_Tensor( N1, N2, N3, N4, N5, N6, N7, N8, M, dl, vol, det1, loc, ratio )
    implicit none
    type(ARES_vector_3D_type), intent(in)  :: N1, N2, N3, N4, N5, N6, N7, N8
    type(ARES_tensor_3D_type), intent(out) :: M
    type(ARES_vector_3D_type), intent(out) :: dl
    real(R8), intent(out)             :: vol
    real(R8), optional, intent(in)    :: det1     ! ghost cell: det of the first cell, whose sign and size it must keep
                                                  ! (0: no such check: the second layer, the first behind a face that
                                                  ! is not a chimera face, and behind a face of zero area)
    integer,  optional, intent(in)    :: loc(5)   ! block, cell i j k, face of a ghost cell (0 otherwise): stop message
    real(R8), optional, intent(in)    :: ratio    ! ghost cell: r = dl_2/dl_1 of the first two cells (stop message)
    ! Local
    integer  :: h
    real(R8) :: det, A(3,3), cofactor(3,3), scale
    logical  :: bad, first, check_vol
    real(R8) :: r(3,8), ea(3,3), v
    real(R8), volatile :: rv(3,8)   ! the nodes relative to N1, rounded once as written (the volume, below)

    ! Compute average cell dimension. A is M^-1, inverse of metric tensor.
    ! M1 = [ xcs, ycs, zcs ; xet, yet, zet ; xzi, yzi, zzi ].
    
    A(1,:) = N5 % c - N1 % c + &
             N6 % c - N2 % c + &
             N7 % c - N3 % c + &
             N8 % c - N4 % c
  
    A(2,:) = N3 % c - N1 % c + &
             N4 % c - N2 % c + &
             N7 % c - N5 % c + &
             N8 % c - N6 % c
             
    A(3,:) = N2 % c - N1 % c + &
             N4 % c - N3 % c + &
             N6 % c - N5 % c + &
             N8 % c - N7 % c

    A = 0.25d0 * A

    ! Determinant of M^-1.
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    ! Volume of the cell: the integral of det J of the trilinear map of its nodes, whose faces are the bilinear faces of
    ! Compute_Norm_Area, shared with the cells next door (docs/theory/numerics.md, Cell volume):
    !   V = 1/6 sum_d e_d . (S_d- + S_d+),
    ! e_d the mean edge of direction d and S_d-, S_d+ the area vectors, (diagonal x diagonal)/2, of the two faces of
    ! direction d, both along +d. It is exact on every cell, twisted faces included, and the volumes of the cells sum to the
    ! volume enclosed by the boundary faces. Signed: negative on a left-handed cell, stored in absolute value below.
    ! Only differences of nodes enter, so the rounding is that of the size of the cell, not of its coordinates, as long as
    ! the compiler rounds the differences as written. With a fast floating-point model (gfortran -fassociative-math, ifx
    ! -fp-model fast) it may rewrite a difference of differences across statements, (q - c1) - (c2 - c1) as q - c2, which
    ! rounds at the size of the coordinates. A VOLATILE object is defined and referenced each time
    ! as the program says (Fortran 2018, 8.5.19), so the nodes relative to N1 are stored in rv and read back into r, and every
    ! later expression has terms of the size of the cell only. A pure procedure cannot hold VOLATILE locals (C1589).
    rv(:,1) = 0d0
    rv(:,2) = N2 % c - N1 % c
    rv(:,3) = N3 % c - N1 % c
    rv(:,4) = N4 % c - N1 % c
    rv(:,5) = N5 % c - N1 % c
    rv(:,6) = N6 % c - N1 % c
    rv(:,7) = N7 % c - N1 % c
    rv(:,8) = N8 % c - N1 % c
    r = rv
    ea(:,1) = 0.25d0 * ( r(:,5) + ( r(:,6) - r(:,2) ) + ( r(:,7) - r(:,3) ) + ( r(:,8) - r(:,4) ) )
    ea(:,2) = 0.25d0 * ( r(:,3) + ( r(:,4) - r(:,2) ) + ( r(:,7) - r(:,5) ) + ( r(:,8) - r(:,6) ) )
    ea(:,3) = 0.25d0 * ( r(:,2) + ( r(:,4) - r(:,3) ) + ( r(:,6) - r(:,5) ) + ( r(:,8) - r(:,7) ) )
    ! The faces, corners in the cyclic order that turns about +d: i N1 N3 N4 N2 and N5 N7 N8 N6, j N1 N2 N6 N5 and
    ! N3 N4 N8 N7, k N1 N5 N7 N3 and N2 N6 N8 N4; twice their area vectors are (3rd - 1st) x (4th - 2nd)
    v =     dot_product ( ea(:,1), Cross3 ( r(:,4), r(:,2) - r(:,3) ) + Cross3 ( r(:,8) - r(:,5), r(:,6) - r(:,7) ) )
    v = v + dot_product ( ea(:,2), Cross3 ( r(:,6), r(:,5) - r(:,2) ) + Cross3 ( r(:,8) - r(:,3), r(:,7) - r(:,4) ) )
    v = v + dot_product ( ea(:,3), Cross3 ( r(:,7), r(:,3) - r(:,5) ) + Cross3 ( r(:,8) - r(:,2), r(:,4) - r(:,6) ) )
    v = v / 12d0

    ! Degenerate cell: |det| against the product of the lengths of the rows of A, its bound (Hadamard's inequality:
    ! equal for an orthogonal cell, whatever its stretching). A ghost cell of the
    ! first layer must also keep the sign and the size of det of the first cell, to 1e-12 (a block may be left-handed:
    ! Compute_Norm_Area takes the orientation of the normals from its first cell): a row of A that rounding leaves tiny
    ! but still in its direction passes the first test, not this one. A mesh error, or ghost cells extrapolated from a
    ! wall-thin first cell whose coordinates were read back in single precision.
    ! This runs inside the OpenMP threads of Setup_Metrics, where a STOP ends the process while the other threads go on:
    ! a bad cell raises bad_metric, the first thread of the rank that finds one writes it, M is the identity (finite for
    ! what the region still computes), and Setup_Metrics stops every rank after the region (Stop_On_Bad_Metric).
    scale = norm2 ( A(1,:) ) * norm2 ( A(2,:) ) * norm2 ( A(3,:) )
    bad = .not. ( abs(det) > 1d-12 * scale )
    if ( present(det1) ) then
      if ( det1 /= 0d0 ) bad = bad .or. .not. ( det * det1 > 1d-12 * det1**2 )
    end if
    ! The volume has the sign of det: det J has the sign of the block in every point of a cell that is not turned over,
    ! and so has its integral. A cell turned over away from its centre can keep det at the centre and give a volume of
    ! the other sign or zero. Not checked where det1 = 0 (the second ghost layer, the first behind a face that is not a
    ! chimera face, and behind a face of zero area), whose size alone is checked.
    check_vol = .true.
    if ( present(det1) ) check_vol = det1 /= 0d0
    if ( check_vol ) bad = bad .or. .not. ( v * det > 0d0 )
    first = .false.
    if ( bad ) then
      first = First_Bad_Metric ( )
      M % c = 0d0
      do h = 1, 3
        M % c(h,h) = 1d0
      enddo
    end if
    if ( first ) then
      if ( present(loc) ) then
        if ( loc(5) > 0 ) then
          write(*,'(A,I0,A,3(1X,I0),A,I0,A)') ' [ERROR] Metric tensor: block ', loc(1), ', cell', loc(2:4), &
                                               ', ghost cell behind face ', loc(5), ': degenerate or inverted'
          if ( present(ratio) ) write(*,'(A,F0.4)') ' [ERROR]   r = dl_2/dl_1 of the first two cells across the face: ', ratio
        else
          write(*,'(A,I0,A,3(1X,I0),A)') ' [ERROR] Metric tensor: block ', loc(1), ', cell', loc(2:4), &
                                          ': degenerate or inverted'
        end if
      else
        write(*,'(A)') ' [ERROR] Metric tensor: degenerate or inverted cell'
      end if
      write(*,'(A,ES12.4,A,ES12.4,A)') ' [ERROR]   det = ', det, ', |A1| |A2| |A3| = ', scale, ' (equal for an orthogonal cell)'
      write(*,'(A,ES12.4,A)') ' [ERROR]   volume = ', v, ' (the integral of det J, of the sign of det on a cell not turned over)'
      if ( present(det1) ) then
        if ( det1 /= 0d0 ) write(*,'(A,ES12.4)') ' [ERROR]   det of the first cell next to the face: ', det1
        write(*,'(A)') ' [ERROR]   The ghost cells take the edge across the face of the first cell, and their edges along the'
        write(*,'(A)') ' [ERROR]   face are extrapolated: behind a chimera face (102), whose viscous flux takes its metric, the'
        write(*,'(A)') ' [ERROR]   first layer turns over where those edges change sharply across the first cell (on a plane'
        write(*,'(A)') ' [ERROR]   2-D mesh, growing 3 times or more); behind any face a layer of zero size stops too.'
        write(*,'(A)') ' [ERROR]   Remedy: index lines that change smoothly across the first two node layers at that face.'
      else
        write(*,'(A)') ' [ERROR]   Check the mesh, or the precision of the restart file.'
      end if
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N1:', N1 % c
      write(*,'(A,3(1X,ES14.6))') ' [ERROR]   corner node N8:', N8 % c

    else if ( .not. bad ) then

      cofactor(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))
      cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      cofactor(1,3) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))
      cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      cofactor(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))
      cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      cofactor(3,1) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))
      cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      cofactor(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))

      M % c = cofactor / det

    endif

    ! Cell length in i,j,k
    do h = 1, 3
      dl % c(h) = sqrt( A(h,1)**2 + A(h,2)**2 + A(h,3)**2 )
    enddo

    vol = abs ( v )

  end subroutine Compute_Metric_Tensor

  
  subroutine Compute_Yn ( blk, bc, nb, nbound )
    use ARES_Lib_KDTree,     only: kdtree_t, kdtree_build, kdtree_free
    use ARES_Mod_MPI,        only: is_local_block, mpi_is_root, mpi_stop_all
    use ARES_Config_Types_m, only: obj_rans
    implicit none
    integer, intent(in)                                  :: nb, nbound
    type(ARES_block_type), dimension(nb), intent(inout)  :: blk
    type(ARES_bc_type), dimension(nbound), intent(in)    :: bc
    ! Local
    real(R8), allocatable :: facex(:), facey(:), facez(:), facekr(:), facer(:)
    real(R8) :: corner(3,4)
    integer :: Fm, Bm, Im, Jm, Km, l, b, nwall, n
    logical, dimension(nb) :: needs_yn
    logical :: all_fast
    type(kdtree_t) :: wall_tree
    ! Periodic images of the wall: image n is the wall turned about the x axis by img_angle(n), then translated by
    ! img_shift(:,n); wall_box(:,1:2) are the low and high corners of the box around the wall faces. Lattice of the
    ! periodic translations: basis(:,1:nv) and its dual basis dual(:,1:nv) (dual(:,i) . basis(:,j) = 1 if i = j, 0
    ! otherwise); the wall faces, and each cell centre before its search, are moved into the cell of the lattice
    ! around anchor (Wall_Images, Reduced_Point)
    integer               :: nimg, nv
    real(R8), allocatable :: img_shift(:,:), img_angle(:)
    real(R8)              :: wall_box(3,2), basis(3,3), dual(3,3), anchor(3)

    ! Validation switch: checks every yn against a brute-force scan, roughly
    ! doubling setup cost.  Keep .false. in production.
    logical, parameter  :: DEBUG_KDTREE = .false.
    real(R8), parameter :: KDTREE_TOL   = 1.0e-10_R8

    ! The distance is the minimum distance from the cell centre to the wall surface
    ! (NASA Turbulence Modeling Resource, Spalart-Allmaras page: not the distance to
    ! the nearest wall node or face centre; next to a convex corner, the distance to
    ! the corner), the union of the wall faces, each the bilinear face through its four
    ! nodes. On a 2-D axisymmetric wedge it is measured in the meridional plane
    ! (x, r), from the wall curve: that is the distance to the wall surface of
    ! revolution, whatever the angle of the wedge. The wedge: a 2-D mesh with delthe
    ! not zero, of either sign.
    wall_meridional = ndir == 2 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 )

    nwall = 0
    do l = 1, nbound
      if ( Is_Wall(bc(l) % type) ) nwall = nwall + 1
    enddo
    allocate ( wall_corner(3, 0:4, max(nwall,1)) )
    allocate ( facex(max(nwall,1)), facey(max(nwall,1)), facez(max(nwall,1)), facekr(max(nwall,1)), facer(max(nwall,1)) )

    nwall=0

    ! bc vector processing.
    do l = 1, nbound

      if ( Is_Wall(bc(l) % type) ) then

        Bm = bc(l) % b
        Im = bc(l) % i
        Jm = bc(l) % j
        Km = bc(l) % k
        Fm = bc(l) % f
        nwall = nwall+1

        ! Boundary face nodes, in cyclic order around the face.
        corner = Face_Corners ( blk(Bm), Im, Jm, Km, Fm )

        if ( wall_meridional ) then
          ! Faces 1-4 of a wedge cell cross the meridional plane along a segment,
          ! whose ends are nodes 1 and 2 (the same x and r as nodes 4 and 3).
          if ( Fm > 4 ) error stop 'Compute_Yn: a wall on a face of the axisymmetric wedge (faces 5-6) is not possible'
          wall_corner(:,1,nwall) = Meridional ( corner(:,1) )
          wall_corner(:,2,nwall) = Meridional ( corner(:,2) )
          wall_corner(:,3,nwall) = wall_corner(:,2,nwall)
          wall_corner(:,4,nwall) = wall_corner(:,1,nwall)
        else
          wall_corner(:,1:4,nwall) = corner
        end if

        ! Face centre (mean of the corners) and radius: every point of the face lies
        ! within the radius of the centre, with a margin for round-off.
        wall_corner(:,0,nwall) = 0.25d0 * ( wall_corner(:,1,nwall) + wall_corner(:,2,nwall) + &
                                            wall_corner(:,3,nwall) + wall_corner(:,4,nwall) )
        facex(nwall) = wall_corner(1,0,nwall)
        facey(nwall) = wall_corner(2,0,nwall)
        facez(nwall) = wall_corner(3,0,nwall)
        facer(nwall) = 0d0
        do n = 1, 4
          facer(nwall) = max ( facer(nwall), norm2 ( wall_corner(:,n,nwall) - wall_corner(:,0,nwall) ) )
        enddo
        facer(nwall) = facer(nwall) * ( 1d0 + 1d-10 )

        ! Face sand-grain roughness, inherited by the cells it is nearest to (with SA-rough only, Wall_Distance_Blk).
        facekr(nwall) = bc(l) % k_rough

      endif

    enddo

    ! After bc processing nwall is the viscous face counter.
    !
    ! Blocks needing yn here: the locally-owned ones, the remote sources of a
    ! local block-connection BC (Yn_Connection reads blk(Bs)%yn), and on root
    ! all of them, as root keeps the arrays of every block.
    ! The rest keep the zero Allocate_Block wrote and are read by no one.
    needs_yn = .false.
    do b = 1, nb
      if ( is_local_block(b) ) needs_yn(b) = .true.
    enddo
    do l = 1, nbound
      if ( (bc(l)%type == 101 .or. bc(l)%type == 201) .and. is_local_block(bc(l)%b) ) &
        needs_yn(bc(l)%bs) = .true.
    enddo
    if ( mpi_is_root ) needs_yn = .true.

    ! No wall anywhere in the domain: leave behind the same sentinel the
    ! brute-force scan did, rather than the k-d tree's huge().
    if ( nwall == 0 ) then
      do b = 1, nb
        if (.not. needs_yn(b)) cycle
        blk(b) % yn ( 1:blk(b)%dim(1), 1:blk(b)%dim(2), 1:blk(b)%dim(3) ) = 1d8
      enddo
      deallocate ( wall_corner )
      return
    end if

    ! Periodic images of the wall. Across a periodic connection (201), and across type
    ! 200 on faces 5-6 of a 3-D mesh without sectors, the mesh goes on translated; across
    ! the rotational periodicity of a 3-D sector (200 on faces 5-6) it goes on turned by
    ! delthe about the x axis: the wall nearest to a cell can be an image of the wall
    ! across those boundaries, and the distance is the minimum over the images too.
    ! The wall faces are moved into one cell of the lattice of the translations, and
    ! so is each cell centre before its search: a few images then suffice, whatever
    ! the size of the mesh against its periods.
    call Wall_Images ( )

    ! Fast path or certified minimum, face by face, on the faces as the search takes them
    ! (moved by Wall_Images). On an axisymmetric wedge the faces are segments of the
    ! meridional plane, whose distance is exact.
    allocate ( wall_fast(nwall), wall_fan_ok(nwall), wall_ghat(nwall), wall_radius(nwall) )
    do n = 1, nwall
      wall_radius(n) = facer(n)
      if ( wall_meridional ) then
        wall_fast(n) = .true.
        wall_fan_ok(n) = .true.
        wall_ghat(n) = 0d0
      else
        call Face_Fast_Test ( wall_corner(:,1:4,n), wall_fast(n), wall_ghat(n), wall_fan_ok(n) )
      end if
    enddo
    all_fast = all ( wall_fast )

    ! One k-d tree over the wall-face centres, with the face radii, shared by every
    ! block and thread: it hands the point-face distance only the faces that can
    ! still be nearer than the nearest found so far (centre within that distance
    ! plus the face radius).
    call kdtree_build ( wall_tree, facex(1:nwall), facey(1:nwall), facez(1:nwall), nwall, facer(1:nwall) )

    ! Compute distance in block b by querying the shared k-d tree.
    do b = 1, nb
      if (.not. needs_yn(b)) cycle
      call Wall_Distance_Blk ( blk(b), wall_tree, facekr(1:nwall) )
    enddo

    ! Validation, off unless DEBUG_KDTREE: prints max |tree - brute| per block
    ! and aborts past KDTREE_TOL.
    if ( DEBUG_KDTREE ) then
      do b = 1, nb
        if (.not. needs_yn(b)) cycle
        call Validate_Tree_Yn ( blk(b), nwall, b )
      enddo
    end if

    call kdtree_free ( wall_tree )
    deallocate ( wall_corner, img_shift, img_angle, wall_fast, wall_fan_ok, wall_ghat, wall_radius )

    contains

      !> Corners of face Fm of cell (Im, Jm, Km) of a block, in cyclic order around the face.
      function Face_Corners ( blk, Im, Jm, Km, Fm ) result ( corner )
        implicit none
        type(ARES_block_type), intent(in) :: blk
        integer, intent(in)               :: Im, Jm, Km, Fm
        real(R8)                          :: corner(3,4)
        ! Local
        integer :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4

        select case(Fm)
          case(1:2)
            i1 = Im - mod(Fm,2)
            j1 = Jm - 1
            k1 = Km - 1
            i2 = Im - mod(Fm,2)
            j2 = Jm
            k2 = Km - 1
            i3 = Im - mod(Fm,2)
            j3 = Jm - 1
            k3 = Km
            i4 = Im - mod(Fm,2)
            j4 = Jm
            k4 = Km
          case(3:4)
            i1 = Im - 1
            j1 = Jm - mod(Fm,2)
            k1 = Km - 1
            i2 = Im
            j2 = Jm - mod(Fm,2)
            k2 = Km - 1
            i3 = Im - 1
            j3 = Jm - mod(Fm,2)
            k3 = Km
            i4 = Im
            j4 = Jm - mod(Fm,2)
            k4 = Km
          case(5:6)
            i1 = Im - 1
            j1 = Jm - 1
            k1 = Km - mod(Fm,2)
            i2 = Im
            j2 = Jm - 1
            k2 = Km - mod(Fm,2)
            i3 = Im - 1
            j3 = Jm
            k3 = Km - mod(Fm,2)
            i4 = Im
            j4 = Jm
            k4 = Km - mod(Fm,2)
        end select

        corner(:,1) = blk % node(i1,j1,k1) % c
        corner(:,2) = blk % node(i2,j2,k2) % c
        corner(:,3) = blk % node(i4,j4,k4) % c
        corner(:,4) = blk % node(i3,j3,k3) % c

      end function Face_Corners


      !> Thickness of cell (Im, Jm, Km) of a block across its face Fm: the distance from the centre of the face to the
      !> centre of the opposite face of the cell, along the normal of the face (the cross product of its diagonals), or
      !> the distance between the two centres when the face has no normal.
      function Cell_Across ( blk, Im, Jm, Km, Fm ) result ( h )
        implicit none
        type(ARES_block_type), intent(in) :: blk
        integer, intent(in)               :: Im, Jm, Km, Fm
        real(R8)                          :: h
        ! Local
        real(R8) :: c(3,4), o(3,4), d(3), nrm(3)

        c = Face_Corners ( blk, Im, Jm, Km, Fm )
        o = Face_Corners ( blk, Im, Jm, Km, merge ( Fm + 1, Fm - 1, mod ( Fm, 2 ) == 1 ) )
        d = 0.25d0 * ( sum ( o, dim=2 ) - sum ( c, dim=2 ) )
        nrm = Cross3 ( c(:,3) - c(:,1), c(:,4) - c(:,2) )
        if ( norm2 ( nrm ) > 0d0 ) then
          h = abs ( dot_product ( d, nrm ) ) / norm2 ( nrm )
        else
          h = norm2 ( d )
        end if

      end function Cell_Across


      !> Images of the wall across the periodic boundaries (img_shift, img_angle, nimg), the lattice of the periodic
      !> translations (nv, basis, dual), the wall faces moved into one cell of it (around anchor), and their box.
      !>  - Translations: each periodic connection (201) gives the vector from its source face to its own face (the
      !>    difference of the face centres), in the meridional plane on a wedge, where it must lie along x. On a 3-D mesh
      !>    without sectors (delthe = 0) so does each block whose k faces are closed by type 200: its last k node plane
      !>    is the first one translated (Check_Mesh_Type), and the cells next to face 5 go on in those next to face 6
      !>    (Ghost_Axysimmetry). The two sides of a periodic boundary match only to the round-off of the coordinates in
      !>    the mesh file, so vectors that differ by at most ETA of the thinnest cell across the periodic boundaries are
      !>    the same: the rule of Check_Mesh_Type, relative to the cells and not to the size of the coordinates, whose
      !>    round-off grows with their distance from the origin (at 1e-6 of the size of the coordinates, a periodic
      !>    channel turned by 30 degrees and written with 5 significant digits stops the run). The thickness of a cell
      !>    across a periodic face is the distance from the centre of the face to the centre of the opposite face of
      !>    the cell, along the normal of the face (Cell_Across); for a block closed in k by type 200, the move of its
      !>    nodes per cell in k (Check_Mesh_Type). The distinct vectors, shortest first, are reduced to a basis of the
      !>    lattice they generate; a vector that is not an integer combination of the basis, within that tolerance
      !>    times the number of vectors combined, stops the run. The basis is then reduced pairwise (Lagrange) and,
      !>    with three translations, by threes (Reduce_By_Two): the same lattice, with shorter vectors.
      !>  - Every wall face is moved, whole, by the vector of the lattice that brings its centre into the cell of the
      !>    lattice around anchor, the centre of the first wall face, where |B*_i . (x - anchor)| <= 1/2 (B* the dual
      !>    basis); Reduced_Point moves each cell centre the same way. The periodic wall is the same set.
      !>  - Images: the integer combinations sum k_i B_i with |k_i| <= floor(1 + rho_i + beta_i), rho_i the largest
      !>    |B*_i . (corner - centre)| of a face and beta_i = (1/2) sum_j |(G^-1)_ij| G_jj, G the Gram matrix of the
      !>    basis. The image that holds the point p of a face nearest to a moved cell centre c is a lattice vector
      !>    closest to v = c - p: v less that image is nearer to 0 than to any +-B_j, so its coordinates along B* are
      !>    at most beta_i in size, and those of v at most 1 + rho_i (docs/theory/turbulence.md). The bound does not
      !>    depend on the size of the mesh against its periods.
      !>  - Rotations: on a 3-D mesh whose sectors are closed by type 200 on faces 5-6 (delthe /= 0), the turns by
      !>    m delthe, |m delthe| <= pi; then every translation must lie along x, the axis of the turns, which keep the
      !>    coordinates along x.
      !>  - A number of images that a default integer cannot hold stops the run.
      subroutine Wall_Images ( )
        implicit none
        ! Local
        integer, parameter    :: MAXV = 3
        real(R8), parameter   :: ETA = 1d-2, PI = 3.141592653589793d0
        real(R8), allocatable :: raw(:,:), cand(:,:), clen(:)
        integer,  allocatable :: order(:)
        logical,  allocatable :: closed(:)
        real(R8)              :: ginv(MAXV,MAXV), coef(MAXV), rhs(MAXV), rho(MAXV), beta(MAXV)
        real(R8)              :: t(3), cface(3,4), hmin, tol, lim, mr, images
        integer               :: nraw, nc, kmax(MAXV), mrot, i, j, m, n, l, b, k1, k2, k3, kk(MAXV)
        logical               :: sector, span, changed
        character(len=512)    :: msg

        ! Translations: of the periodic connections, and of the blocks closed in k by type 200 on a 3-D mesh without
        ! sectors; hmin, the thinnest cell across them
        span = ndir == 3 .and. delthe >= 0d0 .and. delthe <= 0d0
        nraw = 0
        do l = 1, nbound
          if ( bc(l) % type == 201 ) nraw = nraw + 1
        enddo
        if ( span ) nraw = nraw + nb
        allocate ( raw(3, max(nraw,1)), cand(3, max(nraw,1)), clen(max(nraw,1)), order(max(nraw,1)) )
        nraw = 0
        hmin = huge(1d0)
        do l = 1, nbound
          if ( bc(l) % type /= 201 ) cycle
          cface = Face_Corners ( blk(bc(l)%b), bc(l)%i, bc(l)%j, bc(l)%k, bc(l)%f )
          t = 0.25d0 * sum ( cface, dim=2 )
          cface = Face_Corners ( blk(bc(l)%bs), bc(l)%is, bc(l)%js, bc(l)%ks, bc(l)%fs )
          nraw = nraw + 1
          raw(:,nraw) = t - 0.25d0 * sum ( cface, dim=2 )
          hmin = min ( hmin, Cell_Across ( blk(bc(l)%b), bc(l)%i, bc(l)%j, bc(l)%k, bc(l)%f ), &
                             Cell_Across ( blk(bc(l)%bs), bc(l)%is, bc(l)%js, bc(l)%ks, bc(l)%fs ) )
        enddo
        if ( span ) then
          allocate ( closed(nb) )
          closed = .false.
          do l = 1, nbound
            if ( bc(l) % type == 200 .and. bc(l) % f >= 5 ) closed(bc(l) % b) = .true.
          enddo
          do b = 1, nb
            if ( .not. closed(b) ) cycle
            nraw = nraw + 1
            raw(:,nraw) = blk(b) % node(0,0,blk(b)%dim(3)) % c - blk(b) % node(0,0,0) % c
            ! The move of the nodes per cell in k, as Check_Mesh_Type takes it: along the k lines, so that a block closed
            ! on itself (a full turn about the x axis, translated by zero) gives its cells, not the round-off of its
            ! seam, and its translation, zero or round-off, falls below the tolerance
            hmin = min ( hmin, K_Move_Per_Cell ( blk(b), blk(b) % dim(3) ) )
          enddo
          deallocate ( closed )
        end if
        tol = 0d0
        if ( nraw > 0 ) tol = ETA * hmin
        nc = 0
        do l = 1, nraw
          t = raw(:,l)
          if ( norm2 ( t ) <= tol ) cycle
          if ( wall_meridional ) then
            if ( abs ( t(2) ) > tol .or. abs ( t(3) ) > tol ) &
              call mpi_stop_all ( '[ERROR] Wall distance: a periodic connection (201) of an axisymmetric wedge is not '// &
                                   'a translation along the x axis' )
            t(2:3) = 0d0
          end if
          ! One of the vectors t and -t: the largest component positive
          n = maxloc ( abs ( t ), dim=1 )
          if ( t(n) < 0d0 ) t = -t
          do i = 1, nc
            if ( norm2 ( cand(:,i) - t ) <= tol ) exit
          enddo
          if ( i > nc ) then
            nc = nc + 1
            cand(:,nc) = t
            clen(nc) = norm2 ( t )
          end if
        enddo

        ! Basis of the lattice, shortest vectors first
        do i = 1, nc
          order(i) = i
        enddo
        do i = 2, nc
          j = order(i)
          m = i - 1
          do while ( m >= 1 )
            if ( clen(order(m)) <= clen(j) ) exit
            order(m+1) = order(m)
            m = m - 1
          enddo
          order(m+1) = j
        enddo
        nv = 0
        do i = 1, nc
          t = cand(:,order(i))
          if ( nv > 0 ) then
            ! Least-squares coefficients of t on the basis, and the nearest integer ones
            do m = 1, nv
              rhs(m) = dot_product ( basis(:,m), t )
            enddo
            call Invert_Gram ( basis, nv, ginv )
            coef(1:nv) = matmul ( ginv(1:nv,1:nv), rhs(1:nv) )
            lim = tol * ( 1d0 + sum ( abs ( anint ( coef(1:nv) ) ) ) )
            ! An integer combination of the basis: nothing new
            if ( norm2 ( t - matmul ( basis(:,1:nv), anint ( coef(1:nv) ) ) ) <= lim ) cycle
            ! In the span of the basis, between its integer combinations: the periodic boundaries disagree
            if ( norm2 ( t - matmul ( basis(:,1:nv), coef(1:nv) ) ) <= lim ) then
              write(msg,'(A,3ES12.4,A,ES8.2,A)') '[ERROR] Wall distance: the periodic translation (', t, &
                ') is not an integer combination of the others, to 1e-2 of the thinnest cell across the periodic'// &
                ' boundaries (', hmin, '): the two sides of each periodic boundary must match to that, with the node'// &
                ' coordinates written with enough significant digits'
              call mpi_stop_all ( trim(msg) )
            end if
          end if
          if ( nv == MAXV ) then
            write(msg,'(A,ES8.2,A)') '[ERROR] Wall distance: more than three independent periodic translations, to'// &
              ' 1e-2 of the thinnest cell across the periodic boundaries (', hmin, '): the two sides of each periodic'// &
              ' boundary must match to that, with the node coordinates written with enough significant digits'
            call mpi_stop_all ( trim(msg) )
          end if
          nv = nv + 1
          basis(:,nv) = t
        enddo

        ! Rotational periodicity of 3-D sectors
        sector = .false.
        if ( ndir == 3 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) then
          do l = 1, nbound
            if ( bc(l) % type == 200 .and. bc(l) % f >= 5 ) sector = .true.
          enddo
        end if
        if ( sector ) then
          do m = 1, nv
            if ( abs ( basis(2,m) ) > tol .or. abs ( basis(3,m) ) > tol ) &
              call mpi_stop_all ( '[ERROR] Wall distance: a periodic connection (201) of a mesh with rotational '// &
                                   'periodicity (200) is not a translation along the x axis' )
          enddo
        end if

        ! Pairwise (Lagrange) reduction of the basis: B_j <- B_j - m B_i, m the integer nearest to B_i.B_j / B_i.B_i,
        ! while that shortens B_j by a factor below 1 - 1e-12 (so it ends). With three translations, when no pair
        ! shortens, also by threes: B_l <- B_l - m_i B_i - m_j B_j (Reduce_By_Two), which finds a short combination
        ! of all three that no pair gives, as B_1 + B_2 + B_3 of three nearly dependent translations. Unimodular
        ! steps: the same lattice
        changed = nv > 1
        do while ( changed )
          changed = .false.
          do i = 1, nv
            do j = 1, nv
              if ( i == j ) cycle
              mr = anint ( dot_product ( basis(:,i), basis(:,j) ) / dot_product ( basis(:,i), basis(:,i) ) )
              if ( abs ( mr ) < 0.5d0 ) cycle
              t = basis(:,j) - mr * basis(:,i)
              if ( norm2 ( t ) < ( 1d0 - 1d-12 ) * norm2 ( basis(:,j) ) ) then
                basis(:,j) = t
                changed = .true.
              end if
            enddo
          enddo
          if ( nv == 3 .and. .not. changed ) then
            do l = 1, 3
              t = Reduce_By_Two ( basis(:,l), basis(:,mod(l,3)+1), basis(:,mod(l+1,3)+1) )
              if ( norm2 ( t ) < ( 1d0 - 1d-12 ) * norm2 ( basis(:,l) ) ) then
                basis(:,l) = t
                changed = .true.
              end if
            enddo
          end if
        enddo

        ! Dual basis, B*_i = sum_j (G^-1)_ij B_j, and beta_i = (1/2) sum_j |(G^-1)_ij| G_jj, bound on the coordinates
        ! along B*_i of a vector nearer to 0 than to any +-B_j (1/2 on an orthogonal basis)
        if ( nv > 0 ) call Invert_Gram ( basis, nv, ginv )
        do i = 1, nv
          dual(:,i) = 0d0
          beta(i) = 0d0
          do j = 1, nv
            dual(:,i) = dual(:,i) + ginv(i,j) * basis(:,j)
            beta(i) = beta(i) + 0.5d0 * abs ( ginv(i,j) ) * dot_product ( basis(:,j), basis(:,j) )
          enddo
        enddo

        ! Every wall face moved, whole, by the vector of the lattice that brings its centre into the cell of the lattice
        ! around the anchor, the centre of the first wall face; a moved face takes its radius again from its moved
        ! corners. rho_i: the largest |B*_i . (corner - centre)| of a face
        anchor = wall_corner(:,0,1)
        rho = 0d0
        if ( nv > 0 ) then
          do n = 1, nwall
            t = Lattice_Shift ( wall_corner(:,0,n) )
            if ( any ( abs ( t ) > 0d0 ) ) then
              do m = 0, 4
                wall_corner(:,m,n) = wall_corner(:,m,n) - t
              enddo
              facex(n) = wall_corner(1,0,n)
              facey(n) = wall_corner(2,0,n)
              facez(n) = wall_corner(3,0,n)
              facer(n) = 0d0
              do m = 1, 4
                facer(n) = max ( facer(n), norm2 ( wall_corner(:,m,n) - wall_corner(:,0,n) ) )
              enddo
              facer(n) = facer(n) * ( 1d0 + 1d-10 )
            end if
            do i = 1, nv
              do m = 1, 4
                rho(i) = max ( rho(i), abs ( dot_product ( dual(:,i), wall_corner(:,m,n) - wall_corner(:,0,n) ) ) )
              enddo
            enddo
          enddo
        end if

        ! Box around the wall faces, as moved
        wall_box(:,1) = minval ( reshape ( wall_corner(:,1:4,1:nwall), [3, 4*nwall] ), dim=2 )
        wall_box(:,2) = maxval ( reshape ( wall_corner(:,1:4,1:nwall), [3, 4*nwall] ), dim=2 )

        ! Bounds of the integer coefficients, |k_i| <= floor(1 + rho_i + beta_i), and of the turns, |m delthe| <= pi;
        ! first the number of images in reals (an upper bound), which a default integer must hold
        images = 1d0
        do i = 1, nv
          images = images * ( 2d0 * aint ( 1d0 + rho(i) + beta(i) ) + 1d0 )
        enddo
        if ( sector ) images = images * ( 2d0 * aint ( PI / abs ( delthe ) + 1d0 ) + 1d0 )
        if ( .not. ( images - 1d0 < real ( huge(nimg), R8 ) ) ) then
          write(msg,'(A,ES10.3,A)') '[ERROR] Wall distance: ', images - 1d0, ' periodic images of the wall, more'// &
            ' than an integer holds'
          call mpi_stop_all ( trim(msg) )
        end if
        kmax = 0
        mrot = 0
        do i = 1, nv
          kmax(i) = floor ( 1d0 + rho(i) + beta(i) )
        enddo
        if ( sector ) mrot = ceiling ( PI / abs ( delthe ) )

        ! The images, all but the identity
        nimg = product ( 2*kmax + 1 ) * ( 2*mrot + 1 ) - 1
        allocate ( img_shift(3, max(nimg,1)), img_angle(max(nimg,1)) )
        img_angle = 0d0
        nimg = 0
        do k3 = -kmax(3), kmax(3)
          do k2 = -kmax(2), kmax(2)
            do k1 = -kmax(1), kmax(1)
              kk = [ k1, k2, k3 ]
              t = 0d0
              do m = 1, nv
                t = t + kk(m) * basis(:,m)
              enddo
              do m = -mrot, mrot
                if ( k1 == 0 .and. k2 == 0 .and. k3 == 0 .and. m == 0 ) cycle
                nimg = nimg + 1
                img_shift(:,nimg) = t
                img_angle(nimg) = m * delthe
              enddo
            enddo
          enddo
        enddo
        deallocate ( raw, cand, clen, order )

      end subroutine Wall_Images


      !> ginv = inverse of the Gram matrix of basis(:,1:nv), nv <= 3 (independent vectors).
      subroutine Invert_Gram ( basis, nv, ginv )
        implicit none
        real(R8), intent(in)  :: basis(:,:)
        integer, intent(in)   :: nv
        real(R8), intent(out) :: ginv(:,:)
        ! Local
        integer  :: i, j
        real(R8) :: gram(3,3), det

        do i = 1, nv
          do j = 1, nv
            gram(i,j) = dot_product ( basis(:,i), basis(:,j) )
          enddo
        enddo
        select case ( nv )
          case ( 1 )
            ginv(1,1) = 1d0 / gram(1,1)
          case ( 2 )
            det = gram(1,1) * gram(2,2) - gram(1,2) * gram(2,1)
            ginv(1,1) =  gram(2,2) / det
            ginv(2,2) =  gram(1,1) / det
            ginv(1,2) = -gram(1,2) / det
            ginv(2,1) = -gram(2,1) / det
          case ( 3 )
            ginv(1,1) = gram(2,2)*gram(3,3) - gram(2,3)*gram(3,2)
            ginv(1,2) = gram(1,3)*gram(3,2) - gram(1,2)*gram(3,3)
            ginv(1,3) = gram(1,2)*gram(2,3) - gram(1,3)*gram(2,2)
            ginv(2,1) = gram(2,3)*gram(3,1) - gram(2,1)*gram(3,3)
            ginv(2,2) = gram(1,1)*gram(3,3) - gram(1,3)*gram(3,1)
            ginv(2,3) = gram(1,3)*gram(2,1) - gram(1,1)*gram(2,3)
            ginv(3,1) = gram(2,1)*gram(3,2) - gram(2,2)*gram(3,1)
            ginv(3,2) = gram(1,2)*gram(3,1) - gram(1,1)*gram(3,2)
            ginv(3,3) = gram(1,1)*gram(2,2) - gram(1,2)*gram(2,1)
            det = gram(1,1)*ginv(1,1) + gram(1,2)*ginv(2,1) + gram(1,3)*ginv(3,1)
            ginv(1:3,1:3) = ginv(1:3,1:3) / det
        end select

      end subroutine Invert_Gram


      !> One step of the reduction by threes of Wall_Images: the shortest of bk - m_i bi - m_j bj over the integers
      !> m_i, m_j nearest to the least-squares coefficients of bk on bi, bj and next to them (bk itself when none is
      !> shorter).
      pure function Reduce_By_Two ( bk, bi, bj ) result ( best )
        implicit none
        real(R8), intent(in) :: bk(3), bi(3), bj(3)
        real(R8)             :: best(3)
        ! Local
        real(R8) :: gii, gij, gjj, ri, rj, det, ci, cj, t(3)
        integer  :: a, b

        gii = dot_product ( bi, bi )
        gij = dot_product ( bi, bj )
        gjj = dot_product ( bj, bj )
        ri  = dot_product ( bi, bk )
        rj  = dot_product ( bj, bk )
        det = gii * gjj - gij * gij
        ci  = anint ( ( gjj * ri - gij * rj ) / det )
        cj  = anint ( ( gii * rj - gij * ri ) / det )
        best = bk
        do a = -1, 1
          do b = -1, 1
            t = bk - ( ci + a ) * bi - ( cj + b ) * bj
            if ( norm2 ( t ) < norm2 ( best ) ) best = t
          enddo
        enddo

      end function Reduce_By_Two


      !> Distance from p to the box wall_box (0 inside).
      pure function Box_Distance ( p ) result ( d )
        implicit none
        real(R8), intent(in) :: p(3)
        real(R8)             :: d

        d = norm2 ( max ( wall_box(:,1) - p, 0d0, p - wall_box(:,2) ) )

      end function Box_Distance


      !> The point at which the image n of the wall is seen as the wall itself: p translated by -img_shift(:,n) and
      !> turned by -img_angle(n) about the x axis.
      pure function Image_Point ( p, n ) result ( q )
        implicit none
        real(R8), intent(in) :: p(3)
        integer, intent(in)  :: n
        real(R8)             :: q(3)
        ! Local
        real(R8) :: ca, sa

        q = p - img_shift(:,n)
        if ( .not. ( img_angle(n) >= 0d0 .and. img_angle(n) <= 0d0 ) ) then
          ca = cos ( img_angle(n) )
          sa = sin ( img_angle(n) )
          q(2:3) = [ q(2) * ca + q(3) * sa, - q(2) * sa + q(3) * ca ]
        end if

      end function Image_Point


      !> The vector of the lattice of the periodic translations that brings the point p into the cell of the lattice
      !> around anchor: sum_i k_i B_i, k_i the integer nearest to B*_i . (p - anchor); zero without translations.
      pure function Lattice_Shift ( p ) result ( s )
        implicit none
        real(R8), intent(in) :: p(3)
        real(R8)             :: s(3)
        ! Local
        integer :: i

        s = 0d0
        do i = 1, nv
          s = s + anint ( dot_product ( dual(:,i), p - anchor ) ) * basis(:,i)
        enddo

      end function Lattice_Shift


      !> The point p moved into the cell of the lattice of the periodic translations around anchor, where Wall_Images
      !> moves the wall faces; p itself without translations.
      pure function Reduced_Point ( p ) result ( q )
        implicit none
        real(R8), intent(in) :: p(3)
        real(R8)             :: q(3)

        q = p - Lattice_Shift ( p )

      end function Reduced_Point


      subroutine Wall_Distance_Blk ( blk, tree, face_kr )
        use ARES_Lib_KDTree, only: kdtree_t, kdtree_nearest_object
        implicit none
        type(ARES_block_type), intent(inout) :: blk
        type(kdtree_t), intent(in)           :: tree
        real(R8), intent(in)                 :: face_kr(:)   ! Roughness of each wall face, by original face index
        ! Local
        integer :: i, j, k, iface, n, iup
        real(R8) :: center(3), dist, q(3), upper

        !$omp parallel
        !$omp do collapse(3) private(i, j, k, center, dist, iface, n, q, upper, iup)
        do k = 1, blk % dim(3)
        do j = 1, blk % dim(2)
        do i = 1, blk % dim(1)

          ! The cell centre, moved into the cell of the lattice of the periodic translations that holds the wall faces
          center = Reduced_Point ( Cell_Centre ( blk, i, j, k ) )

          if ( all_fast ) then

            ! Every face is plane to the resolution of its coordinates: its four triangles are the face
            ! (Wall_Face_Distance).
            call kdtree_nearest_object ( tree, center(1), center(2), center(3), Wall_Face_Distance, dist, iface )

            ! Periodic images: an image can be nearer only if the box of the wall, seen
            ! from the moved point, is not farther than the distance found so far. The
            ! query goes on from that distance and face, so ties keep the lowest face.
            do n = 1, nimg
              q = Image_Point ( center, n )
              if ( Box_Distance ( q ) > dist ) cycle
              call kdtree_nearest_object ( tree, q(1), q(2), q(3), Wall_Face_Distance, dist, iface, resume=.true. )
            enddo

          else

            ! First pass, over the faces and their images: an upper bound of the distance, the four triangles of a face
            ! plus the bound of their gap from the face and of the rounding (Wall_Upper_Distance).
            fanc_now = fanc_now + 1
            wall_image = 0
            call kdtree_nearest_object ( tree, center(1), center(2), center(3), Wall_Upper_Distance, upper, iup )
            do n = 1, nimg
              q = Image_Point ( center, n )
              if ( Box_Distance ( q ) > upper ) cycle
              wall_image = n
              call kdtree_nearest_object ( tree, q(1), q(2), q(3), Wall_Upper_Distance, upper, iup, resume=.true. )
            enddo

            ! Second pass: the distance to the bilinear faces (Wall_Exact_Distance), with the search starting from the
            ! bound of the first pass (2^-30 of it added, for the triangles of a face of the fast path computed again)
            ! and no face: the tree prunes by that bound from the start, and the distance is computed only on the faces
            ! whose triangles, less the same bound of their gap and of the rounding, are not farther than it, nor than
            ! the least distance computed so far. The face of the minimum is never left out, and each face computed
            ! gives a pure function of the face, so the distance and its face do not depend on the tree, the threads
            ! or the ranks.
            wall_upper = huge ( 1d0 )
            if ( upper < 0.5d0 * huge ( 1d0 ) ) wall_upper = upper + 2d0**(-30) * upper
            wall_least = huge ( 1d0 )
            dist = wall_upper
            iface = huge ( 1 )
            wall_image = 0
            call kdtree_nearest_object ( tree, center(1), center(2), center(3), Wall_Exact_Distance, dist, iface, &
                                         resume=.true. )
            do n = 1, nimg
              q = Image_Point ( center, n )
              if ( Box_Distance ( q ) > dist ) cycle
              wall_image = n
              call kdtree_nearest_object ( tree, q(1), q(2), q(3), Wall_Exact_Distance, dist, iface, resume=.true. )
            enddo
            if ( iface == huge ( 1 ) ) then
              ! No face below the bound: the four triangles of a face of the fast path rounded differently in the two
              ! passes, by more than the 2^-30 added. Search again with no bound.
              wall_upper = huge ( 1d0 )
              wall_least = huge ( 1d0 )
              wall_image = 0
              call kdtree_nearest_object ( tree, center(1), center(2), center(3), Wall_Exact_Distance, dist, iface )
              do n = 1, nimg
                q = Image_Point ( center, n )
                if ( Box_Distance ( q ) > dist ) cycle
                wall_image = n
                call kdtree_nearest_object ( tree, q(1), q(2), q(3), Wall_Exact_Distance, dist, iface, resume=.true. )
              enddo
            end if

          end if

          blk % yn (i,j,k) = dist
          ! The roughness of the nearest wall face, which only SA-rough models (the other cells keep ks = 0)
          if ( obj_rans % rough ) blk % ks (i,j,k) = face_kr(iface)

        enddo; enddo; enddo
        !$omp end parallel

      end subroutine Wall_Distance_Blk


      !> Brute-force check of the k-d tree result: recompute each cell's
      !> nearest-wall distance against the full face list and its periodic
      !> images, compare with `blk%yn`, and abort past KDTREE_TOL.
      subroutine Validate_Tree_Yn ( blk, nwall_, b_id )
        use ARES_Mod_MPI, only: mpi_rank_
        implicit none
        type(ARES_block_type), intent(in)           :: blk
        integer, intent(in)                         :: nwall_
        integer, intent(in)                         :: b_id
        ! Local
        integer :: i, j, k, ll, n
        real(R8) :: center(3), q(3), dbrute, dummy, err, max_err

        max_err = 0d0

        !$omp parallel
        !$omp do collapse(3) private(i, j, k, ll, n, center, q, dbrute, dummy, err) &
        !$omp             reduction(max:max_err)
        do k = 1, blk % dim(3)
        do j = 1, blk % dim(2)
        do i = 1, blk % dim(1)

          center = Reduced_Point ( Cell_Centre ( blk, i, j, k ) )

          dbrute = 1d8
          do n = 0, nimg
            q = center
            if ( n > 0 ) q = Image_Point ( center, n )
            do ll = 1, nwall_
              dummy = Wall_Face_Exact ( ll, q(1), q(2), q(3) )
              if (dummy < dbrute) dbrute = dummy
            enddo
          enddo

          err = abs( blk % yn(i,j,k) - dbrute )
          if (err > max_err) max_err = err

        enddo; enddo; enddo
        !$omp end parallel

        write(*,'(A,I0,A,I0,A,ES12.5)') &
          ' [KDTREE] rank ', mpi_rank_, '  block ', b_id, '  max |tree - brute| = ', max_err
        if (max_err > KDTREE_TOL) then
          write(*,'(A,ES12.5,A)') ' [KDTREE] ERROR: exceeds tolerance ', KDTREE_TOL, '. Aborting.'
          stop 1
        end if

      end subroutine Validate_Tree_Yn

  end subroutine Compute_Yn


  !> Point where the wall distance of cell (i,j,k) is measured: the mean of its
  !> eight nodes, in the meridional plane (x, r, 0) on an axisymmetric wedge.
  function Cell_Centre ( blk, i, j, k ) result ( center )
    implicit none
    type(ARES_block_type), intent(in) :: blk
    integer, intent(in)               :: i, j, k
    real(R8)                          :: center(3)
    ! Local
    integer :: a, b, c

    center = 0d0
    do c = k-1, k
    do b = j-1, j
    do a = i-1, i
      if ( wall_meridional ) then
        center = center + Meridional ( blk % node(a,b,c) % c )
      else
        center = center + blk % node(a,b,c) % c
      end if
    enddo; enddo; enddo
    center = center / 8d0

  end function Cell_Centre


  !> (x, y, z) -> (x, r, 0), with r = sqrt(y^2 + z^2) the distance from the x axis.
  pure function Meridional ( p ) result ( q )
    implicit none
    real(R8), intent(in) :: p(3)
    real(R8)             :: q(3)

    q = [ p(1), sqrt ( p(2)**2 + p(3)**2 ), 0d0 ]

  end function Meridional


  !> Distance from (x,y,z) to wall face i of Compute_Yn (kdtree_object_distance):
  !> the segment of the face in the meridional plane on an axisymmetric wedge;
  !> otherwise the four triangles between the face centre and the sides, which are
  !> the face itself when it is plane, as on every 2-D and wedge mesh.
  function Wall_Face_Distance ( i, x, y, z ) result ( d )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: d
    ! Local
    real(R8) :: p(3)

    p = [ x, y, z ]
    if ( wall_meridional ) then
      d = Point_Segment_Distance ( p, wall_corner(:,1,i), wall_corner(:,2,i) )
    else
      d = min ( Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,1,i), wall_corner(:,2,i) ), &
                Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,2,i), wall_corner(:,3,i) ), &
                Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,3,i), wall_corner(:,4,i) ), &
                Point_Triangle_Distance ( p, wall_corner(:,0,i), wall_corner(:,4,i), wall_corner(:,1,i) ) )
    end if

  end function Wall_Face_Distance


  !> Distance from p to the segment ab (a point if a = b).
  pure function Point_Segment_Distance ( p, a, b ) result ( d )
    implicit none
    real(R8), intent(in) :: p(3), a(3), b(3)
    real(R8)             :: d
    ! Local
    real(R8) :: ab(3), ap(3), l2, t

    ab = b - a
    ap = p - a
    l2 = dot_product ( ab, ab )
    t  = 0d0
    if ( l2 > 0d0 ) t = min ( max ( dot_product ( ap, ab ) / l2, 0d0 ), 1d0 )
    d = norm2 ( ap - t * ab )

  end function Point_Segment_Distance


  !> Distance from p to the triangle abc. The triangle is convex: if the projection
  !> of p on its plane falls inside it (on the inner side of the three sides), the
  !> distance is the one from the plane; otherwise the nearest point lies on a side.
  !> A triangle with no area next to its sides (a degenerate face) has no plane and
  !> is treated as its sides.
  pure function Point_Triangle_Distance ( p, a, b, c ) result ( d )
    implicit none
    real(R8), intent(in) :: p(3), a(3), b(3), c(3)
    real(R8)             :: d
    ! Local
    real(R8) :: n(3), n2, e2

    n  = Cross3 ( b - a, c - a )
    n2 = dot_product ( n, n )
    e2 = max ( sum ( (b - a)**2 ), sum ( (c - b)**2 ), sum ( (a - c)**2 ) )
    if ( n2 > ( 1d-12 * e2 )**2 ) then
      if ( dot_product ( Cross3 ( b - a, p - a ), n ) >= 0d0 .and. &
           dot_product ( Cross3 ( c - b, p - b ), n ) >= 0d0 .and. &
           dot_product ( Cross3 ( a - c, p - c ), n ) >= 0d0 ) then
        d = abs ( dot_product ( p - a, n ) ) / sqrt ( n2 )
        return
      end if
    end if
    d = min ( Point_Segment_Distance ( p, a, b ), Point_Segment_Distance ( p, b, c ), &
              Point_Segment_Distance ( p, c, a ) )

  end function Point_Triangle_Distance

  !> Fast path of the wall distance on the face with corners c(:,1:4) in cyclic order (docs/theory/turbulence.md, Wall
  !> distance). The four triangles between the mean of the corners and the sides (Wall_Face_Distance) lie within the gap
  !> g* = min(|T|/16, |T.n|/(2|n|)) of the face, T = c1 - c2 + c3 - c4 its twist and n = (c3 - c1) x (c4 - c2) the cross
  !> product of its diagonals, the second bound when the face projected on the plane normal to n is convex (every corner
  !> between sides of non-zero length turns the same way). The two bounds are computed together with eps1, eps2, proven
  !> bounds of their rounding (u = 2^-53, every operation within 2u, sums in any order, fused multiply-adds), and the face
  !> takes the fast path when one of them is within 8 u M plus its eps, M the largest |coordinate| of the corners: the face
  !> is plane to the resolution of its coordinates. Every parallelogram takes it, and every convex plane face but those with
  !> a corner within (49 + 30 kappa_n) u of a straight or a null angle, kappa_n = |c3 - c1| |c4 - c2| / |n|; the distance
  !> of a face so taken is within 8 u M + 2 eps of the distance to the face, plus the rounding of the triangles. ghat is an
  !> upper bound of g*, for the filter of the second pass of Wall_Distance_Blk; fan_ok: in each triangle, in the coordinates
  !> local to c1 of Local_Fan_Distance, 1/sin of the angle at the mean of the corners is below 1e5, under which the margin
  !> of that filter covers the rounding of the triangles. The sides and the diagonals are differences of the corners
  !> rounded once, as written: VOLATILE, as in Compute_Metric_Tensor (and so not pure: C1589).
  subroutine Face_Fast_Test ( c, fast, ghat, fan_ok )
    implicit none
    real(R8), intent(in)  :: c(3,4)
    logical,  intent(out) :: fast, fan_ok
    real(R8), intent(out) :: ghat
    ! Local
    real(R8), parameter :: U = epsilon(1d0) / 2d0, CNU = 4d0 + 4d0 * sqrt(2d0), CW = 10d0 + 4d0 * sqrt(2d0), &
                           CONDMAX = 1d5
    real(R8), volatile  :: ve(3,4), vd(3,2), vl(3,4)
    real(R8) :: e(3,4), d1(3), d2(3), t(3), n(3), ne(4), nt, nn, nu, m, g, eps, w, epsw, cl(3,4), mc(3), u1(3), u2(3)
    integer  :: k, k1
    logical  :: convex

    do k = 1, 4
      ve(:,k) = c(:,mod(k,4)+1) - c(:,k)
    enddo
    vd(:,1) = c(:,3) - c(:,1)
    vd(:,2) = c(:,4) - c(:,2)
    e = ve
    d1 = vd(:,1)
    d2 = vd(:,2)
    t = - ( e(:,1) + e(:,3) )
    n = Cross3 ( d1, d2 )
    nt = Vnorm ( t )
    nn = Vnorm ( n )
    do k = 1, 4
      ne(k) = Vnorm ( e(:,k) )
    enddo
    m = maxval ( abs ( c ) )

    ! |T|/16, on any face: rounding |T~ - T| <= 2u (|e1| + |e3|) + 2u |T|, of the norm 5u
    g = nt / 16d0
    eps = 2d0 * ( U / 8d0 * ( ne(1) + ne(3) ) + 5d0 * U / 8d0 * nt )
    ghat = g + eps
    fast = g <= 8d0 * U * m + eps

    ! |T.n|/(2|n|), when the projection is convex (corners turning within their rounding are not taken as convex) and
    ! the normal is known to 1/8 of its size: |n~ - n| <= nu = (4 + 4 sqrt 2) u |d1| |d2|
    nu = CNU * U * Vnorm ( d1 ) * Vnorm ( d2 )
    if ( nn > 0d0 .and. nu <= nn / 8d0 ) then
      convex = .true.
      do k = 1, 4
        k1 = mod(k,4) + 1
        w = dot_product ( Cross3 ( e(:,k), e(:,k1) ), n )
        epsw = 2d0 * ne(k) * ne(k1) * ( CW * U * nn + nu )
        if ( .not. ( w >= epsw ) ) convex = .false.
      enddo
      if ( convex ) then
        g = abs ( dot_product ( t, n ) ) / ( 2d0 * nn )
        eps = 2d0 * ( U * ( ne(1) + ne(3) ) + 9d0 * U * nt + nt * nu / nn )
        ghat = min ( ghat, g + eps )
        fast = fast .or. g <= 8d0 * U * m + eps
      end if
    end if

    ! Conditioning of the four triangles in the coordinates local to c1
    do k = 1, 4
      vl(:,k) = c(:,k) - c(:,1)
    enddo
    cl = vl
    mc = 0.25d0 * ( cl(:,1) + cl(:,2) + cl(:,3) + cl(:,4) )
    fan_ok = .true.
    do k = 1, 4
      u1 = cl(:,k) - mc
      u2 = cl(:,mod(k,4)+1) - mc
      if ( .not. ( Vnorm ( u1 ) * Vnorm ( u2 ) < CONDMAX * Vnorm ( Cross3 ( u1, u2 ) ) ) ) fan_ok = .false.
    enddo

  end subroutine Face_Fast_Test


  !> Distance from the point p to the bilinear face with corners c(:,1:4) in cyclic order, x(s,t) = c1 + s a + t b + s t T
  !> (a = c2 - c1, b = c4 - c1, T = c1 - c2 + c3 - c4, s and t in [0,1]): its global minimum, certified within
  !> tau = 256 u S, S = |p - c1| + |a| + |b| + |T|, u = 2^-53 (docs/theory/turbulence.md, Wall distance). Branch and bound
  !> over rectangles of (s,t), depth first, with a stack of fixed size:
  !>  - upper bound best: the distance to points of the face, its four sides (closed form), the centre of each piece and
  !>    projected Newton steps (from the best point of a 5 x 5 grid, and at the end from the best point);
  !>  - lower bound of a piece with corners V00, V10, V11, V01: the larger of the distance to its two triangles across
  !>    V00-V11 less |T_P|/4 (T_P the twist of the piece: piece and triangles are within |T_P|/4 of each other), and the
  !>    distance to the convex hull of its four corners, which holds the piece; both less proven bounds of their rounding
  !>    (Tri_Dist_Lb, Orient);
  !>  - a piece is dropped when its lower bound is within tau of best, resolved when the distance to its centre is, and
  !>    otherwise split in two across its longer side in space; at most 128 levels.
  !> d = best, the distance to a point of the face, within tau + 32 u S of the minimum when no piece is cut at the deepest
  !> level; cert = best less the least lower bound of the pieces left (at most tau then), ncut the pieces cut at that level,
  !> npiece the pieces visited. The coordinates local to c1 are differences rounded once, as written (VOLATILE, as in
  !> Compute_Metric_Tensor), and every quantity of the bounds derives from them; not pure (C1589).
  subroutine Bilinear_Face_Distance ( p, c, d, cert, ncut, npiece )
    implicit none
    real(R8), intent(in)            :: p(3), c(3,4)
    real(R8), intent(out)           :: d
    real(R8), intent(out), optional :: cert
    integer,  intent(out), optional :: ncut, npiece
    ! Local
    real(R8), parameter :: U = epsilon(1d0) / 2d0, KTAU = 256d0
    integer,  parameter :: MAXLEV = 128, MAXSTK = MAXLEV + 8
    real(R8), volatile  :: vq(3), va(3), vb(3), ve(3)
    real(R8) :: q(3), a(3), b(3), e(3), t(3), nt, scale, tau, best, bs, bt, zero(3), v(3,4), lb, lbmin, ub, ls, lt, h2
    real(R8) :: pm0(3), pm1(3), kv(3,4,2), klb(2), ks(2), kt(2), khs(2), kht(2), s0, t0, hs, ht, dd, xg(3), si, ti
    real(R8) :: sk_s0(MAXSTK), sk_t0(MAXSTK), sk_hs(MAXSTK), sk_ht(MAXSTK), sk_lb(MAXSTK), sk_v(3,4,MAXSTK)
    integer  :: sk_lev(MAXSTK), nsk, lev, i, j, kk, nc, npc, ncu, o1, o2

    vq = p - c(:,1)
    va = c(:,2) - c(:,1)
    vb = c(:,4) - c(:,1)
    ve = c(:,3) - c(:,1)
    q = vq
    a = va
    b = vb
    e = ve
    t = e - a - b
    nt = Vnorm ( t )
    scale = Vnorm ( q ) + Vnorm ( a ) + Vnorm ( b ) + nt
    tau = KTAU * U * scale
    zero = 0d0

    ! Upper bound: the four sides, then projected Newton from the best point of a 5 x 5 grid
    best = min ( Point_Segment_Distance ( q, zero, a ), Point_Segment_Distance ( q, a, e ), &
                 Point_Segment_Distance ( q, e, b ), Point_Segment_Distance ( q, b, zero ) )
    bs = 0.5d0
    bt = 0.5d0
    dd = huge ( 1d0 )
    s0 = 0d0
    t0 = 0d0
    do j = 0, 4
      do i = 0, 4
        si = 0.25d0 * i
        ti = 0.25d0 * j
        xg = Bil ( si, ti )
        if ( Vnorm ( xg - q ) < dd ) then
          dd = Vnorm ( xg - q )
          s0 = si
          t0 = ti
        end if
      enddo
    enddo
    call Polish ( s0, t0 )

    ! Branch and bound
    v(:,1) = zero
    v(:,2) = a
    v(:,3) = e
    v(:,4) = b
    npc = 1
    ncu = 0
    lbmin = huge ( 1d0 )
    nsk = 1
    sk_s0(1) = 0d0
    sk_t0(1) = 0d0
    sk_hs(1) = 1d0
    sk_ht(1) = 1d0
    sk_lev(1) = 0
    sk_v(:,:,1) = v
    sk_lb(1) = Lower ( 1d0, 1d0, v, best - tau )
    do while ( nsk > 0 )
      s0 = sk_s0(nsk)
      t0 = sk_t0(nsk)
      hs = sk_hs(nsk)
      ht = sk_ht(nsk)
      lev = sk_lev(nsk)
      v = sk_v(:,:,nsk)
      lb = sk_lb(nsk)
      nsk = nsk - 1
      if ( lb >= best - tau ) then
        lbmin = min ( lbmin, lb )
        cycle
      end if
      ub = Vnorm ( Bil ( s0 + 0.5d0 * hs, t0 + 0.5d0 * ht ) - q )
      if ( ub < best ) then
        best = ub
        bs = s0 + 0.5d0 * hs
        bt = t0 + 0.5d0 * ht
      end if
      if ( ub - lb <= tau .or. lev >= MAXLEV ) then
        if ( lev >= MAXLEV .and. ub - lb > tau ) ncu = ncu + 1
        lbmin = min ( lbmin, lb )
        cycle
      end if
      ! split in two across the longer side in space (the longer of its two sides along s against the longer along t)
      ls = max ( Vnorm ( v(:,2) - v(:,1) ), Vnorm ( v(:,3) - v(:,4) ) )
      lt = max ( Vnorm ( v(:,4) - v(:,1) ), Vnorm ( v(:,3) - v(:,2) ) )
      if ( ls >= lt ) then
        h2 = 0.5d0 * hs
        pm0 = Bil ( s0 + h2, t0 )
        pm1 = Bil ( s0 + h2, t0 + ht )
        ks = [ s0, s0 + h2 ]
        kt = [ t0, t0 ]
        khs = h2
        kht = ht
        kv(:,1,1) = v(:,1)
        kv(:,2,1) = pm0
        kv(:,3,1) = pm1
        kv(:,4,1) = v(:,4)
        kv(:,1,2) = pm0
        kv(:,2,2) = v(:,2)
        kv(:,3,2) = v(:,3)
        kv(:,4,2) = pm1
      else
        h2 = 0.5d0 * ht
        pm0 = Bil ( s0, t0 + h2 )
        pm1 = Bil ( s0 + hs, t0 + h2 )
        ks = [ s0, s0 ]
        kt = [ t0, t0 + h2 ]
        khs = hs
        kht = h2
        kv(:,1,1) = v(:,1)
        kv(:,2,1) = v(:,2)
        kv(:,3,1) = pm1
        kv(:,4,1) = pm0
        kv(:,1,2) = pm0
        kv(:,2,2) = pm1
        kv(:,3,2) = v(:,3)
        kv(:,4,2) = v(:,4)
      end if
      nc = 0
      do kk = 1, 2
        klb(kk) = Lower ( khs(kk), kht(kk), kv(:,:,kk), best - tau )
        npc = npc + 1
      enddo
      ! the child with the larger lower bound goes on the stack first, so the smaller one is taken first
      o1 = 1
      o2 = 2
      if ( klb(1) < klb(2) ) then
        o1 = 2
        o2 = 1
      end if
      do kk = 1, 2
        nc = merge ( o1, o2, kk == 1 )
        if ( klb(nc) >= best - tau ) then
          lbmin = min ( lbmin, klb(nc) )
          cycle
        end if
        if ( nsk >= MAXSTK ) then            ! not reached (the stack holds at most one more piece per level)
          ncu = ncu + 1
          lbmin = min ( lbmin, klb(nc) )
          cycle
        end if
        nsk = nsk + 1
        sk_s0(nsk) = ks(nc)
        sk_t0(nsk) = kt(nc)
        sk_hs(nsk) = khs(nc)
        sk_ht(nsk) = kht(nc)
        sk_lev(nsk) = lev + 1
        sk_v(:,:,nsk) = kv(:,:,nc)
        sk_lb(nsk) = klb(nc)
      enddo
    enddo
    call Polish ( bs, bt )

    d = best
    if ( present(cert) ) cert = best - min ( lbmin, best )
    if ( present(ncut) ) ncut = ncu
    if ( present(npiece) ) npiece = npc

  contains

    !> Point (s, t) of the face, in the local coordinates
    function Bil ( s, tt ) result ( x )
      real(R8), intent(in) :: s, tt
      real(R8)             :: x(3)
      x = s * a + tt * b + ( s * tt ) * t
    end function Bil

    !> Projected Newton steps on |x(s,t) - q|^2 from (sa, ta), clamped to the unit square: they only lower best
    subroutine Polish ( sa, ta )
      real(R8), intent(in) :: sa, ta
      real(R8) :: sc, tc, s1, t1, r(3), xs(3), xt(3), g1, g2, h11, h22, h12, det, ds, dt, d1, lam
      integer  :: it
      sc = sa
      tc = ta
      do it = 1, 30
        r = Bil ( sc, tc ) - q
        xs = a + tc * t
        xt = b + sc * t
        g1 = dot_product ( r, xs )
        g2 = dot_product ( r, xt )
        h11 = dot_product ( xs, xs )
        h22 = dot_product ( xt, xt )
        h12 = dot_product ( xs, xt ) + dot_product ( r, t )
        det = h11 * h22 - h12 * h12
        if ( det > 0d0 .and. h11 > 0d0 ) then
          ds = - ( h22 * g1 - h12 * g2 ) / det
          dt = - ( h11 * g2 - h12 * g1 ) / det
        else
          lam = max ( h11, h22, 1d-300 )
          ds = - g1 / lam
          dt = - g2 / lam
        end if
        if ( .not. ( abs ( ds ) <= 2d0 .and. abs ( dt ) <= 2d0 ) ) then
          ds = sign ( min ( abs ( ds ), 2d0 ), ds )
          dt = sign ( min ( abs ( dt ), 2d0 ), dt )
          if ( .not. ( abs ( ds ) <= 2d0 .and. abs ( dt ) <= 2d0 ) ) exit      ! not a number
        end if
        s1 = min ( max ( sc + ds, 0d0 ), 1d0 )
        t1 = min ( max ( tc + dt, 0d0 ), 1d0 )
        d1 = Vnorm ( Bil ( s1, t1 ) - q )
        if ( d1 < best ) then
          best = d1
          bs = s1
          bt = t1
        end if
        if ( abs ( s1 - sc ) + abs ( t1 - tc ) < 4d0 * U ) exit
        sc = s1
        tc = t1
      enddo
    end subroutine Polish

    !> Lower bound of the distance from q to the piece of sides hs, ht and corners vv (V00, V10, V11, V01). When the
    !> bound (a) already reaches thr, the threshold at which the piece is dropped, it is returned alone: the decisions of
    !> the search, and so the distance, are those of the larger of (a) and (b) (only the certificate can be larger, and
    !> stays within tau)
    function Lower ( hs_, ht_, vv, thr ) result ( lbo )
      real(R8), intent(in) :: hs_, ht_, vv(3,4), thr
      real(R8) :: lbo, da1, da2, db1, db2, lbs, lbh, o, eo, oi, ei, amax, thick, sg, w(3,4)
      integer  :: k
      logical  :: outside
      da1 = Tri_Dist_Lb ( q, vv(:,1), vv(:,2), vv(:,3) )
      da2 = Tri_Dist_Lb ( q, vv(:,1), vv(:,3), vv(:,4) )
      ! (a) the two triangles across V00-V11, within |T_P|/4 = |T| hs ht / 4 of the piece
      lbs = min ( da1, da2 ) - 0.25d0 * nt * hs_ * ht_ * ( 1d0 + 16d0 * U )
      if ( lbs >= thr ) then
        lbo = lbs
        return
      end if
      db1 = Tri_Dist_Lb ( q, vv(:,1), vv(:,2), vv(:,4) )
      db2 = Tri_Dist_Lb ( q, vv(:,2), vv(:,3), vv(:,4) )
      ! (b) the convex hull of the four corners: 0 unless q is certainly outside it (a sub-determinant certainly of the
      ! sign opposite to the whole); a hull flat to the rounding lies within its thickness of the four triangles
      call Orient ( vv(:,1), vv(:,2), vv(:,3), vv(:,4), o, eo )
      if ( .not. ( abs ( o ) > eo ) ) then
        amax = 0.5d0 * max ( Vnorm ( Cross3 ( vv(:,2) - vv(:,1), vv(:,3) - vv(:,1) ) ), &
                             Vnorm ( Cross3 ( vv(:,3) - vv(:,1), vv(:,4) - vv(:,1) ) ), &
                             Vnorm ( Cross3 ( vv(:,2) - vv(:,1), vv(:,4) - vv(:,1) ) ), &
                             Vnorm ( Cross3 ( vv(:,3) - vv(:,2), vv(:,4) - vv(:,2) ) ) )
        thick = 0d0
        if ( amax > 0d0 ) thick = ( abs ( o ) + eo ) / ( 2d0 * amax )
        lbh = min ( da1, da2, db1, db2 ) - thick
      else
        sg = sign ( 1d0, o )
        outside = .false.
        do k = 1, 4
          w = vv
          w(:,k) = q
          call Orient ( w(:,1), w(:,2), w(:,3), w(:,4), oi, ei )
          if ( sg * oi < - ei ) outside = .true.
        enddo
        lbh = 0d0
        if ( outside ) lbh = min ( da1, da2, db1, db2 )
      end if
      lbo = max ( lbs, lbh )
    end function Lower

  end subroutine Bilinear_Face_Distance


  !> A lower bound of the distance from p to the triangle abc in binary64: the distance of Point_Triangle_Distance less a
  !> bound of its rounding, every operation within 2u (docs/theory/turbulence.md, Wall distance). Where the triangle has a
  !> plane, n = (b - a) x (c - a): (32 cond + 20) u (|p - a| + e), cond = |b - a| |c - a| / |n|, e the longest side.
  !> Where Point_Triangle_Distance takes its sides (|n| at most 1e-12 e^2): every point of the triangle lies within its
  !> height |n| / e over the longest side, so that height plus 16 u e + 20 u (|p - a| + e); on the border between the two
  !> (|n|^2 within a factor 4 of the threshold, where the two evaluations of n may fall on either side), the larger of the
  !> two; three corners on one point, 20 u |p - a|.
  pure function Tri_Dist_Lb ( p, a, b, c ) result ( d )
    implicit none
    real(R8), intent(in) :: p(3), a(3), b(3), c(3)
    real(R8)             :: d
    ! Local
    real(R8), parameter :: U = epsilon(1d0) / 2d0
    real(R8) :: ab(3), ac(3), n(3), n2, e2, pa, e, mplane, mside

    ab = b - a
    ac = c - a
    n = Cross3 ( ab, ac )
    n2 = dot_product ( n, n )
    e2 = max ( dot_product ( ab, ab ), dot_product ( c - b, c - b ), dot_product ( ac, ac ) )
    pa = Vnorm ( p - a )
    e = sqrt ( e2 )
    d = Point_Triangle_Distance ( p, a, b, c )
    mplane = 0d0
    if ( n2 > 0d0 ) mplane = ( 32d0 * Vnorm ( ab ) * Vnorm ( ac ) / sqrt ( n2 ) + 20d0 ) * U * ( pa + e )
    mside = 20d0 * U * pa
    if ( e2 > 0d0 ) mside = sqrt ( n2 ) / e + 16d0 * U * e + 20d0 * U * ( pa + e )
    if ( n2 > 4d0 * ( 1d-12 * e2 )**2 ) then
      d = d - mplane
    else if ( n2 > 0.25d0 * ( 1d-12 * e2 )**2 ) then
      d = d - max ( mplane, mside )
    else
      d = d - mside
    end if

  end function Tri_Dist_Lb


  !> o = det [b - a, c - a, d - a] and eo, a bound of its rounding, every operation within 2u: 32 u |b - a| |c - a| |d - a|
  pure subroutine Orient ( a, b, c, d, o, eo )
    implicit none
    real(R8), intent(in)  :: a(3), b(3), c(3), d(3)
    real(R8), intent(out) :: o, eo
    ! Local
    real(R8), parameter :: U = epsilon(1d0) / 2d0
    real(R8) :: u1(3), u2(3), u3(3)

    u1 = b - a
    u2 = c - a
    u3 = d - a
    o = dot_product ( u1, Cross3 ( u2, u3 ) )
    eo = 32d0 * U * Vnorm ( u1 ) * Vnorm ( u2 ) * Vnorm ( u3 )

  end subroutine Orient


  !> |v| = sqrt(v . v), a sum of squares and a root: within 5u of the length, every operation within 2u
  pure function Vnorm ( v ) result ( r )
    implicit none
    real(R8), intent(in) :: v(3)
    real(R8)             :: r

    r = sqrt ( dot_product ( v, v ) )

  end function Vnorm


  !> The four triangles of wall face i of Compute_Yn between the mean of its corners and its sides, in the coordinates local
  !> to its first corner (differences rounded once, VOLATILE): the distance from (x,y,z) for the filter of the second pass
  !> of Wall_Distance_Blk, whose margin (Filter_Margin) covers their rounding at that scale (not pure: C1589).
  function Local_Fan_Distance ( i, x, y, z ) result ( d )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: d
    ! Local
    real(R8), volatile :: vq(3), vc(3,4)
    real(R8) :: q(3), cl(3,4), mc(3)
    integer  :: k

    vq = [ x, y, z ] - wall_corner(:,1,i)
    do k = 1, 4
      vc(:,k) = wall_corner(:,k,i) - wall_corner(:,1,i)
    enddo
    q = vq
    cl = vc
    mc = 0.25d0 * ( cl(:,1) + cl(:,2) + cl(:,3) + cl(:,4) )
    d = min ( Point_Triangle_Distance ( q, mc, cl(:,1), cl(:,2) ), Point_Triangle_Distance ( q, mc, cl(:,2), cl(:,3) ), &
              Point_Triangle_Distance ( q, mc, cl(:,3), cl(:,4) ), Point_Triangle_Distance ( q, mc, cl(:,4), cl(:,1) ) )

  end function Local_Fan_Distance


  !> Local_Fan_Distance of wall face i from (x,y,z) through the table of the point (fanc_*): the value met for this face in
  !> the image wall_image of the same point, when the table holds it, else computed and kept. Either way the same function of
  !> the face and the point (the second pass of a point queries the points of the first). Linear probing from a hash of the
  !> face and the image; a slot of another point is free.
  function Fan_Cached ( i, x, y, z ) result ( d )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: d
    ! Local
    integer :: h, k, hk, slot

    h = mod ( mod ( i, FANC ) * 31 + mod ( wall_image, FANC ) * 977, FANC )
    slot = -1
    do k = 0, FANC - 1
      hk = mod ( h + k, FANC )
      if ( fanc_tag(hk) /= fanc_now ) then
        slot = hk
        exit
      end if
      if ( fanc_face(hk) == i .and. fanc_img(hk) == wall_image ) then
        d = fanc_val(hk)
        return
      end if
    enddo
    d = Local_Fan_Distance ( i, x, y, z )
    if ( slot >= 0 ) then
      fanc_tag(slot) = fanc_now
      fanc_face(slot) = i
      fanc_img(slot) = wall_image
      fanc_val(slot) = d
    end if

  end function Fan_Cached


  !> Margin of the filter of the second pass for wall face i seen from (x,y,z): 2^-30 (|q - m| + r), m the mean of the
  !> corners and r the radius of the face. It covers the rounding of Local_Fan_Distance on triangles whose angle at the mean
  !> has 1/sin below 1e5, 3 (25 . 1e5 + 16) u (|q - m| + r), that of the local coordinates, and the 288 u S of the certified
  !> minimum (S <= |q - m| + 9 r); it widens the set of the faces computed, not the result.
  pure function Filter_Margin ( i, x, y, z ) result ( mg )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: mg

    mg = 2d0**(-30) * ( Vnorm ( [ x, y, z ] - wall_corner(:,0,i) ) + wall_radius(i) )

  end function Filter_Margin


  !> First pass of Wall_Distance_Blk: an upper bound of the distance from (x,y,z) to wall face i. On a face of the fast path
  !> its distance (its four triangles); on another, the four triangles in local coordinates plus ghat and the margin; huge
  !> when those triangles are too ill-conditioned for the margin (the face is then computed in the second pass anyway).
  function Wall_Upper_Distance ( i, x, y, z ) result ( d )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: d

    if ( wall_fast(i) ) then
      d = Wall_Face_Distance ( i, x, y, z )
    else if ( wall_fan_ok(i) ) then
      d = Fan_Cached ( i, x, y, z ) + wall_ghat(i) + Filter_Margin ( i, x, y, z )
    else
      d = huge ( 1d0 )
    end if

  end function Wall_Upper_Distance


  !> Second pass of Wall_Distance_Blk: the distance from (x,y,z) to wall face i (Wall_Face_Exact), or huge when the face
  !> cannot be nearer than min(wall_upper, wall_least), the bound of the first pass and the least distance computed so far
  !> in this pass: its four triangles in local coordinates, less ghat and the margin, are farther. Both are at least the
  !> minimum when the second pass ends with a face (wall_least always; wall_upper, else the pass is made again with no
  !> bound), so the face of the minimum is never left out, whatever the order of the search. A face of the fast path is
  !> always computed.
  function Wall_Exact_Distance ( i, x, y, z ) result ( d )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: d

    if ( .not. wall_fast(i) .and. wall_fan_ok(i) ) then
      if ( Fan_Cached ( i, x, y, z ) - wall_ghat(i) - Filter_Margin ( i, x, y, z ) > &
           min ( wall_upper, wall_least ) ) then
        d = huge ( 1d0 )
        return
      end if
    end if
    d = Wall_Face_Exact ( i, x, y, z )
    wall_least = min ( wall_least, d )

  end function Wall_Exact_Distance


  !> Distance from (x,y,z) to wall face i of Compute_Yn: its four triangles on a face of the fast path (Wall_Face_Distance:
  !> the face to the resolution of its coordinates, or the meridional segment on a wedge), the certified minimum over the
  !> bilinear face otherwise (Bilinear_Face_Distance).
  function Wall_Face_Exact ( i, x, y, z ) result ( d )
    implicit none
    integer, intent(in)  :: i
    real(R8), intent(in) :: x, y, z
    real(R8)             :: d

    if ( wall_fast(i) ) then
      d = Wall_Face_Distance ( i, x, y, z )
    else
      call Bilinear_Face_Distance ( [ x, y, z ], wall_corner(:,1:4,i), d )
    end if

  end function Wall_Face_Exact



  pure function Cross3 ( u, v ) result ( w )
    implicit none
    real(R8), intent(in) :: u(3), v(3)
    real(R8)             :: w(3)

    w = [ u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1) ]

  end function Cross3


  subroutine Yn_Connection ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, blks )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm, Is, Js, Ks
    type(ARES_block_type), intent(inout) :: blks
    type(ARES_block_type), intent(inout) :: blkm
    ! Local
    integer :: Ig, Jg, Kg

    ! Ghost cell coordinates
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    blkm % yn(Ig,Jg,Kg) = blks % yn(Is,Js,Ks)
    blkm % ks(Ig,Jg,Kg) = blks % ks(Is,Js,Ks)

  end subroutine Yn_Connection


  subroutine Yn_Chimera ( Im, Jm, Km, Fm, blk )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blk
    ! Local
    integer :: Ig, Jg, Kg

    ! Ghost cell coordinates
    Ig = Im - guide(Fm,1)
    Jg = Jm - guide(Fm,2)
    Kg = Km - guide(Fm,3)

    ! The chimera ghost cell extends the receiver grid (as its metrics do), so it
    ! takes the wall distance and roughness of the adjacent receiver cell.
    ! Zeroth order on purpose: linear extrapolation can reach yn <= 0 near walls.
    blk % yn(Ig,Jg,Kg) = blk % yn(Im,Jm,Km)
    blk % ks(Ig,Jg,Kg) = blk % ks(Im,Jm,Km)

  end subroutine Yn_Chimera


  !> Metrics of the gc ghost cells behind face Fm of cell (Im,Jm,Km) of block b, from nodes extrapolated linearly
  !> along the index lines, N(-g) = 2 N(-g+1) - N(-g+2): every ghost layer has the edge across the face of the first
  !> cell, so its size across the face never vanishes nor turns over; its edges along the face are extrapolated too, and
  !> turn a layer over where they change sharply across the first cell (on a plane 2-D mesh, growing 3 times or more
  !> turns the first layer, 5/3 times the second). Behind a chimera face (chimera true) the first layer must keep the
  !> sign and the size of the first cell (Compute_Metric_Tensor stops otherwise): its metric is the ghost metric of the
  !> viscous flux of chimera faces (BC_Connection_Visc), the one flux that reads the metric of an extrapolated ghost
  !> cell. Behind the other faces no flux reads it: the walls (301, 302) take the metric of the cell, the other faces
  !> carry no viscous flux that needs it, and the wedge (200) rotates the cell,
  !> exactly when its node planes lie at theta = -+delthe/2 (otherwise the ghost cell is flat and the check of its size
  !> stops the run). Of the ghost layers the solver then reads only the length across the face, in the MUSCL reconstruction
  !> (the interior faces next to the boundary that of the first layer, BC_Connection_Eul of chimera faces that of both),
  !> and that is the length of the first cell whatever the edges along the face: there their sign is not
  !> checked, only a zero size stops (no infinite metric). The first layer is centred at
  !> 2 x_f - x_1 (x_f the mean of the nodes of the face, x_1 that of the cell): the position of the ghost values with a
  !> wall value (2 Z_wall - Z_cell, RANS_Extrapolate_Wall), and the MUSCL distance l0 of the first interior face is dl_1.
  !> The 2-D axisymmetric wedge rotates the cell (exact).
  subroutine BC_Extrapolate_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg, b, chimera )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: blk
    type(ARES_tensor_3D_type), intent(out)    :: Mg(2)
    type(ARES_vector_3D_type), intent(out)    :: dlg(2)
    real(R8), intent(out)                :: volg(2)
    integer, intent(in)                  :: b     ! the block, for the stop message of Compute_Metric_Tensor
    logical, intent(in)                  :: chimera   ! a chimera face (102): the sign of the first layer is checked
    ! Local
    integer :: g, loc(5), d, x(3)
    real(R8) :: d1, rr
    type(ARES_vector_3D_type), dimension(-gc:gc):: N1, N2, N3, N4, N5, N6, N7, N8

    ! det of the first cell, 1/det(M), whose sign and size its ghost cells must keep (Compute_Metric_Tensor)
    associate ( Mc => blk % M(Im,Jm,Km) % c )
      d1 = 1d0 / ( Mc(1,1) * ( Mc(2,2)*Mc(3,3) - Mc(2,3)*Mc(3,2) ) - Mc(1,2) * ( Mc(2,1)*Mc(3,3) - Mc(2,3)*Mc(3,1) ) &
                 + Mc(1,3) * ( Mc(2,1)*Mc(3,2) - Mc(2,2)*Mc(3,1) ) )
    end associate
    loc = [b, Im, Jm, Km, Fm]
    ! r = dl_2/dl_1 of the first two cells across the face, for the stop message (0 with one cell)
    d  = (Fm + 1) / 2
    rr = 0d0
    if ( blk % dim(d) >= 2 ) rr = blk % dl(Im+guide(Fm,1),Jm+guide(Fm,2),Km+guide(Fm,3)) % c(d) / blk % dl(Im,Jm,Km) % c(d)
    ! Behind a face of zero area (an axis, as Compute_Norm_Area finds it) the ghost cells are the reflection of the
    ! first cell through the axis: det = -d1 by construction. No check of their sign there (d1 = 0)
    x = [Im, Jm, Km]
    if ( mod(Fm,2) == 1 ) x(d) = x(d) - 1
    if ( blk % dir(d) % f(x(1),x(2),x(3)) % a == 0d0 ) d1 = 0d0
    ! Behind a face that is not a chimera face no flux reads the metric of the first layer: no check of its sign
    ! either, as for the second layer (docs/theory/numerics.md, metrics of the ghost cells)
    if ( .not. chimera ) d1 = 0d0

    select case (Fm)
      case(1)
        do g = 0, 1
          N1(g)%c = blk % node (Im+g-1,Jm-1,Km-1) % c
          N2(g)%c = blk % node (Im+g-1,Jm-1,Km  ) % c
          N3(g)%c = blk % node (Im+g-1,Jm  ,Km-1) % c
          N4(g)%c = blk % node (Im+g-1,Jm  ,Km  ) % c
        enddo
        do g = -1, -gc, -1
          N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
          N2(g)%c = 2d0*N2(g+1)%c - N2(g+2)%c
          N3(g)%c = 2d0*N3(g+1)%c - N3(g+2)%c
          N4(g)%c = 2d0*N4(g+1)%c - N4(g+2)%c
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g+1), N2(g+1), N3(g+1), N4(g+1), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(-g)   = blk % M(Im+g, Jm, Km)
          dlg(-g)  = blk % dl(Im+g, Jm, Km)
          volg(-g) = blk % vol(Im+g, Jm, Km)
        enddo
        
      case(2)
        do g = 0, -1, -1
          N5(g)%c = blk % node (Im+g,Jm-1,Km-1) % c
          N6(g)%c = blk % node (Im+g,Jm-1,Km  ) % c
          N7(g)%c = blk % node (Im+g,Jm  ,Km-1) % c
          N8(g)%c = blk % node (Im+g,Jm  ,Km  ) % c
        enddo
        do g = 1, gc
          N5(g)%c = 2d0*N5(g-1)%c - N5(g-2)%c
          N6(g)%c = 2d0*N6(g-1)%c - N6(g-2)%c
          N7(g)%c = 2d0*N7(g-1)%c - N7(g-2)%c
          N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blk % M(Im+g, Jm, Km), blk % dl(Im+g, Jm, Km), blk % vol(Im+g, Jm, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(g)   = blk % M(Im+g, Jm, Km)
          dlg(g)  = blk % dl(Im+g, Jm, Km)
          volg(g) = blk % vol(Im+g, Jm, Km)
        enddo

      case(3)
        do g = 0, 1
          N1(g)%c = blk % node (Im-1,Jm+g-1,Km-1) % c
          N2(g)%c = blk % node (Im-1,Jm+g-1,Km  ) % c
          N5(g)%c = blk % node (Im  ,Jm+g-1,Km-1) % c
          N6(g)%c = blk % node (Im  ,Jm+g-1,Km  ) % c
        enddo
        do g = -1, -gc, -1
          N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
          N2(g)%c = 2d0*N2(g+1)%c - N2(g+2)%c
          N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
          N6(g)%c = 2d0*N6(g+1)%c - N6(g+2)%c
          call Compute_Metric_Tensor ( N1(g), N2(g), N1(g+1), N2(g+1), N5(g), N6(g), N5(g+1), N6(g+1), &
                                       blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(-g)   = blk % M(Im, Jm+g, Km)
          dlg(-g)  = blk % dl(Im, Jm+g, Km)
          volg(-g) = blk % vol(Im, Jm+g, Km)
        enddo

      case(4)
        do g = 0, -1, -1
          N3(g)%c = blk % node (Im-1,Jm+g,Km-1) % c
          N4(g)%c = blk % node (Im-1,Jm+g,Km  ) % c
          N7(g)%c = blk % node (Im  ,Jm+g,Km-1) % c
          N8(g)%c = blk % node (Im  ,Jm+g,Km  ) % c
        enddo
        do g = 1, gc
          N3(g)%c = 2d0*N3(g-1)%c - N3(g-2)%c
          N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
          N7(g)%c = 2d0*N7(g-1)%c - N7(g-2)%c
          N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
          call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                       blk % M(Im, Jm+g, Km), blk % dl(Im, Jm+g, Km), blk % vol(Im, Jm+g, Km), &
                                       det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
          Mg(g)   = blk % M(Im, Jm+g, Km)
          dlg(g)  = blk % dl(Im, Jm+g, Km)
          volg(g) = blk % vol(Im, Jm+g, Km)
        enddo

      case(5)
        if ( ndir==1 ) then
          do g = 1, gc
            Mg(g)   = blk % M(Im, Jm, Km)
            dlg(g)  = blk % dl(Im, Jm, Km)
            volg(g) = blk % vol(Im, Jm, Km)
          enddo
        elseif ( ndir==2 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c(1) = N1(0) % c(1) 
            N1(g)%c(2) = N1(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N1(g)%c(3) = N1(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N3(g)%c(1) = N3(0) % c(1) 
            N3(g)%c(2) = N3(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N3(g)%c(3) = N3(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N5(g)%c(1) = N5(0) % c(1) 
            N5(g)%c(2) = N5(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N5(g)%c(3) = N5(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            N7(g)%c(1) = N7(0) % c(1) 
            N7(g)%c(2) = N7(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0-g))
            N7(g)%c(3) = N7(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0-g))
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo         
        else
          do g = 0, 1
            N1(g)%c = blk % node (Im-1,Jm-1,Km+g-1) % c
            N3(g)%c = blk % node (Im-1,Jm  ,Km+g-1) % c
            N5(g)%c = blk % node (Im  ,Jm-1,Km+g-1) % c
            N7(g)%c = blk % node (Im  ,Jm  ,Km+g-1) % c
          enddo
          do g = -1, -gc, -1
            N1(g)%c = 2d0*N1(g+1)%c - N1(g+2)%c
            N3(g)%c = 2d0*N3(g+1)%c - N3(g+2)%c
            N5(g)%c = 2d0*N5(g+1)%c - N5(g+2)%c
            N7(g)%c = 2d0*N7(g+1)%c - N7(g+2)%c
            call Compute_Metric_Tensor ( N1(g), N1(g+1), N3(g), N3(g+1), N5(g), N5(g+1), N7(g), N7(g+1), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(-g)   = blk % M(Im, Jm, Km+g)
            dlg(-g)  = blk % dl(Im, Jm, Km+g)
            volg(-g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      case(6)
        if ( ndir==1 ) then
          do g = 1, gc
            Mg(g)   = blk % M(Im, Jm, Km)
            dlg(g)  = blk % dl(Im, Jm, Km)
            volg(g) = blk % vol(Im, Jm, Km)
          enddo
        elseif ( ndir==2 .and. .not. ( delthe >= 0d0 .and. delthe <= 0d0 ) ) then
          ! 2Dax: extrapolation with a rotation angle delthe (exact).
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c(1) = N2(0) % c(1) 
            N2(g)%c(2) = N2(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N2(g)%c(3) = N2(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N4(g)%c(1) = N4(0) % c(1) 
            N4(g)%c(2) = N4(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N4(g)%c(3) = N4(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N6(g)%c(1) = N6(0) % c(1) 
            N6(g)%c(2) = N6(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N6(g)%c(3) = N6(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            N8(g)%c(1) = N8(0) % c(1) 
            N8(g)%c(2) = N8(0) % c(2) / cos(delthe*0.5d0) * cos(delthe*(0.5d0+g))
            N8(g)%c(3) = N8(0) % c(3) / sin(delthe*0.5d0) * sin(delthe*(0.5d0+g))
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo         
        else
          do g = 0, -1, -1
            N2(g)%c = blk % node (Im-1,Jm-1,Km+g) % c
            N4(g)%c = blk % node (Im-1,Jm  ,Km+g) % c
            N6(g)%c = blk % node (Im  ,Jm-1,Km+g) % c
            N8(g)%c = blk % node (Im  ,Jm  ,Km+g) % c
          enddo
          do g = 1, gc
            N2(g)%c = 2d0*N2(g-1)%c - N2(g-2)%c
            N4(g)%c = 2d0*N4(g-1)%c - N4(g-2)%c
            N6(g)%c = 2d0*N6(g-1)%c - N6(g-2)%c
            N8(g)%c = 2d0*N8(g-1)%c - N8(g-2)%c
            call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                         blk % M(Im, Jm, Km+g), blk % dl(Im, Jm, Km+g), blk % vol(Im, Jm, Km+g), &
                                         det1 = merge ( d1, 0d0, abs(g) == 1 ), loc = loc, ratio = rr )
            Mg(g)   = blk % M(Im, Jm, Km+g)
            dlg(g)  = blk % dl(Im, Jm, Km+g)
            volg(g) = blk % vol(Im, Jm, Km+g)
          enddo
        endif

      end select

  end subroutine BC_Extrapolate_Metrics


  !> Layer g of the ghost cells of a symmetry face (300): a copy of the cell g - 1 behind the face, in M, dl and vol of
  !> the block and in the copies of the record, Mg(g), dlg(g), volg(g). In a direction with fewer cells than g that cell
  !> is a ghost cell of the opposite face: Setup_Metrics fills the layers one at a time, every face of a layer before
  !> the next.
  subroutine BC_Symmetry_Metrics( Im, Jm, Km, Fm, blk, Mg, dlg, volg, g )
    implicit none
    integer, intent(in) :: Im, Jm, Km, Fm, g
    type(ARES_block_type), intent(inout)     :: blk
    type(ARES_tensor_3D_type), intent(inout) :: Mg(2)
    type(ARES_vector_3D_type), intent(inout) :: dlg(2)
    real(R8), intent(inout)                  :: volg(2)

    select case (Fm)
      case(1)
        blk % M(Im-g,Jm,Km)   = blk % M(Im+g-1,Jm,Km)
        blk % dl(Im-g,Jm,Km)  = blk % dl(Im+g-1,Jm,Km)
        blk % vol(Im-g,Jm,Km) = blk % vol(Im+g-1,Jm,Km)
        Mg(g)   = blk % M(Im-g,Jm,Km)
        dlg(g)  = blk % dl(Im-g,Jm,Km)
        volg(g) = blk % vol(Im-g,Jm,Km)
      case(2)
        blk % M(Im+g,Jm,Km)   = blk % M(Im-g+1,Jm,Km)
        blk % dl(Im+g,Jm,Km)  = blk % dl(Im-g+1,Jm,Km)
        blk % vol(Im+g,Jm,Km) = blk % vol(Im-g+1,Jm,Km)
        Mg(g)   = blk % M(Im+g,Jm,Km)
        dlg(g)  = blk % dl(Im+g,Jm,Km)
        volg(g) = blk % vol(Im+g,Jm,Km)
      case(3)
        blk % M(Im,Jm-g,Km)   = blk % M(Im,Jm+g-1,Km)
        blk % dl(Im,Jm-g,Km)  = blk % dl(Im,Jm+g-1,Km)
        blk % vol(Im,Jm-g,Km) = blk % vol(Im,Jm+g-1,Km)
        Mg(g)   = blk % M(Im,Jm-g,Km)
        dlg(g)  = blk % dl(Im,Jm-g,Km)
        volg(g) = blk % vol(Im,Jm-g,Km)
      case(4)
        blk % M(Im,Jm+g,Km)   = blk % M(Im,Jm-g+1,Km)
        blk % dl(Im,Jm+g,Km)  = blk % dl(Im,Jm-g+1,Km)
        blk % vol(Im,Jm+g,Km) = blk % vol(Im,Jm-g+1,Km)
        Mg(g)   = blk % M(Im,Jm+g,Km)
        dlg(g)  = blk % dl(Im,Jm+g,Km)
        volg(g) = blk % vol(Im,Jm+g,Km)
      case(5)
        blk % M(Im,Jm,Km-g)   = blk % M(Im,Jm,Km+g-1)
        blk % dl(Im,Jm,Km-g)  = blk % dl(Im,Jm,Km+g-1)
        blk % vol(Im,Jm,Km-g) = blk % vol(Im,Jm,Km+g-1)
        Mg(g)   = blk % M(Im,Jm,Km-g)
        dlg(g)  = blk % dl(Im,Jm,Km-g)
        volg(g) = blk % vol(Im,Jm,Km-g)
      case(6)
        blk % M(Im,Jm,Km+g)   = blk % M(Im,Jm,Km-g+1)
        blk % dl(Im,Jm,Km+g)  = blk % dl(Im,Jm,Km-g+1)
        blk % vol(Im,Jm,Km+g) = blk % vol(Im,Jm,Km-g+1)
        Mg(g)   = blk % M(Im,Jm,Km+g)
        dlg(g)  = blk % dl(Im,Jm,Km+g)
        volg(g) = blk % vol(Im,Jm,Km+g)
    end select

  end subroutine BC_Symmetry_Metrics


  subroutine BC_Connect_Metrics ( Im, Jm, Km, Fm, blkm, Is, Js, Ks, Fs, blks, d11s, d12s, d21s, d22s, Mg, dlg, volg, periodic )
    implicit none
    integer, intent(in)                  :: Im, Jm, Km, Fm, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(ARES_block_type), intent(in)    :: blks
    type(ARES_block_type), intent(inout) :: blkm
    type(ARES_tensor_3D_type), intent(out) :: Mg(2)
    type(ARES_vector_3D_type), intent(out) :: dlg(2)
    real(R8), intent(out)                  :: volg(2)
    logical, intent(in)                    :: periodic
    ! Local
    integer      :: g, Is1, Js1, Ks1
    integer      :: II, JJ, KK, III, JJJ, KKK
    integer      :: guidem(6,3), guides(6,3), guidem2(6,3), guides2(6,3)
    type(ARES_vector_3D_type), dimension(0:gc) :: N1, N2, N3, N4, N5, N6, N7, N8
    type(ARES_vector_3D_type) :: T
    
    ! ---------------------------------------------------------------------------------------------
    ! Preliminary definitions
    guidem(1,1)=-2 ; guidem(1,2)= 0 ; guidem(1,3)= 0 ; guidem(2,1)=1 ; guidem(2,2)=0 ; guidem(2,3)=0
    guidem(3,1)= 0 ; guidem(3,2)=-2 ; guidem(3,3)= 0 ; guidem(4,1)=0 ; guidem(4,2)=1 ; guidem(4,3)=0
    guidem(5,1)= 0 ; guidem(5,2)=0  ; guidem(5,3)=-2 ; guidem(6,1)=0 ; guidem(6,2)=0 ; guidem(6,3)=1

    guidem2(1,1)=-3 ; guidem2(1,2)= 0 ; guidem2(1,3)= 0 ; guidem2(2,1)=2 ; guidem2(2,2)=0 ; guidem2(2,3)=0
    guidem2(3,1)= 0 ; guidem2(3,2)=-3 ; guidem2(3,3)= 0 ; guidem2(4,1)=0 ; guidem2(4,2)=2 ; guidem2(4,3)=0
    guidem2(5,1)= 0 ; guidem2(5,2)= 0 ; guidem2(5,3)=-3 ; guidem2(6,1)=0 ; guidem2(6,2)=0 ; guidem2(6,3)=2

    guides(1,1)= 0 ; guides(1,2)= 0 ; guides(1,3)= 0 ; guides(2,1)=-1 ; guides(2,2)= 0 ; guides(2,3)= 0
    guides(3,1)= 0 ; guides(3,2)= 0 ; guides(3,3)= 0 ; guides(4,1)= 0 ; guides(4,2)=-1 ; guides(4,3)= 0
    guides(5,1)= 0 ; guides(5,2)= 0 ; guides(5,3)= 0 ; guides(6,1)= 0 ; guides(6,2)= 0 ; guides(6,3)=-1

    guides2(1,1)= 1 ; guides2(1,2)= 0 ; guides2(1,3)= 0 ; guides2(2,1)=-2 ; guides2(2,2)= 0 ; guides2(2,3)= 0
    guides2(3,1)= 0 ; guides2(3,2)= 1 ; guides2(3,3)= 0 ; guides2(4,1)= 0 ; guides2(4,2)=-2 ; guides2(4,3)= 0
    guides2(5,1)= 0 ; guides2(5,2)= 0 ; guides2(5,3)= 1 ; guides2(6,1)= 0 ; guides2(6,2)= 0 ; guides2(6,3)=-2
    ! ---------------------------------------------------------------------------------------------

    select case (Fm)
      case(1)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
      case(2)
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(3)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
      case(4)
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
      case(5)
        N1(0)%c = blkm%node(Im-1,Jm-1,Km-1) % c
        N3(0)%c = blkm%node(Im-1,Jm  ,Km-1) % c
        N5(0)%c = blkm%node(Im  ,Jm-1,Km-1) % c
        N7(0)%c = blkm%node(Im  ,Jm  ,Km-1) % c
      case(6)
        N2(0)%c = blkm%node(Im-1,Jm-1,Km  ) % c
        N4(0)%c = blkm%node(Im-1,Jm  ,Km  ) % c
        N6(0)%c = blkm%node(Im  ,Jm-1,Km  ) % c
        N8(0)%c = blkm%node(Im  ,Jm  ,Km  ) % c
    end select

    T%c = 0._R8
    if (periodic) then
      Is1 = Is + gc*guides(Fs,1) - guides2(Fs,1)
      Js1 = Js + gc*guides(Fs,2) - guides2(Fs,2)
      Ks1 = Ks + gc*guides(Fs,3) - guides2(Fs,3)
      select case(Fs)
        case(1:2)
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4)
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6)
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      select case(Fm)
        case(1) ; T%c = N4(0)%c - blks%node(II,JJ,KK)%c
        case(2) ; T%c = N8(0)%c - blks%node(II,JJ,KK)%c
        case(3) ; T%c = N6(0)%c - blks%node(II,JJ,KK)%c
        case(4) ; T%c = N8(0)%c - blks%node(II,JJ,KK)%c
        case(5) ; T%c = N7(0)%c - blks%node(II,JJ,KK)%c
        case(6) ; T%c = N8(0)%c - blks%node(II,JJ,KK)%c
      end select
    end if

    do g = 1, gc

      ! connected cell node index
      Is1 = Is + (gc-g)*guides(Fs,1)+(g-1)*guides2(Fs,1)
      Js1 = Js + (gc-g)*guides(Fs,2)+(g-1)*guides2(Fs,2)
      Ks1 = Ks + (gc-g)*guides(Fs,3)+(g-1)*guides2(Fs,3)

      ! compute the indexes in the connected block corresponding to the ghost nodes to update in the local block 
      select case(Fs)
        case(1:2) ! i-faces
          II=Is1
          JJ=(2*Js1+d11s+d21s)/2
          KK=(2*Ks1+d12s+d22s)/2
        case(3:4) ! j-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=Js1
          KK=(2*Ks1+d12s+d22s)/2
        case(5:6) ! k-faces
          II=(2*Is1+d11s+d21s)/2
          JJ=(2*Js1+d12s+d22s)/2
          KK=Ks1
      end select
      if (Fm == 1) then
          N4(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 2) then
          N8(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 3) then
          N6(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 4) then
          N8(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 5) then
          N7(g)%c = blks%node(ii,jj,kk)%c + T%c
      elseif (Fm == 6) then
          N8(g)%c = blks%node(ii,jj,kk)%c + T%c
      endif

      ! update the remaining 3 indexes (only necessary in the corners, but for good measure)
      ! node i1-1,i2 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s
          KKK=KK-d12s
        case(3:4)
          III=II-d11s
          JJJ=JJ
          KKK=KK-d12s
        case(5:6)
          III=II-d11s
          JJJ=JJ-d12s
          KKK=KK
      end select
      if (Fm == 1) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 2) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 3) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 4) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 5) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 6) then
        N4(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      endif

      ! node i1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d21s
          KKK=KK-d22s
        case(3:4)
          III=II-d21s
          JJJ=JJ
          KKK=KK-d22s
        case(5:6)
          III=II-d21s
          JJJ=JJ-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 2) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 3) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 4) then
        N7(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 5) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 6) then
        N6(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      endif

      ! node i1-1,i2-1 update
      select case(Fs)
        case(1:2)
          III=II
          JJJ=JJ-d11s-d21s
          KKK=KK-d12s-d22s
        case(3:4)
          III=II-d11s-d21s
          JJJ=JJ
          KKK=KK-d12s-d22s
        case(5:6)
          III=II-d11s-d21s
          JJJ=JJ-d12s-d22s
          KKK=KK
      end select
      if (Fm == 1) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 2) then
        N5(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 3) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 4) then
        N3(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 5) then
        N1(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      elseif (Fm == 6) then
        N2(g)%c = blks%node(iii,jjj,kkk)%c + T%c
      endif

      ! compute metric variables
      select case(Fm)
        case(1)
          call Compute_Metric_Tensor ( N1(g), N2(g), N3(g), N4(g), N1(g-1), N2(g-1), N3(g-1), N4(g-1), &
                                       blkm%M(Im-g, Jm, Km), blkm%dl(Im-g, Jm, Km), blkm%vol(Im-g, Jm, Km))
          Mg(g)   = blkm % M(Im-g, Jm, Km)
          dlg(g)  = blkm % dl(Im-g, Jm, Km)
          volg(g) = blkm % vol(Im-g, Jm, Km)
        case(2)
          call Compute_Metric_Tensor ( N5(g-1), N6(g-1), N7(g-1), N8(g-1), N5(g), N6(g), N7(g), N8(g), &
                                       blkm%M(Im+g, Jm, Km), blkm%dl(Im+g, Jm, Km), blkm%vol(Im+g, Jm, Km))
          Mg(g)   = blkm % M(Im+g, Jm, Km)
          dlg(g)  = blkm % dl(Im+g, Jm, Km)
          volg(g) = blkm % vol(Im+g, Jm, Km)
        case(3)
          call Compute_Metric_Tensor ( N1(g), N2(g), N1(g-1), N2(g-1), N5(g), N6(g), N5(g-1), N6(g-1), &
                                       blkm%M(Im, Jm-g, Km), blkm%dl(Im, Jm-g, Km), blkm%vol(Im, Jm-g, Km))
          Mg(g)   = blkm % M(Im, Jm-g, Km)
          dlg(g)  = blkm % dl(Im, Jm-g, Km)
          volg(g) = blkm % vol(Im, Jm-g, Km)
        case(4)
          call Compute_Metric_Tensor ( N3(g-1), N4(g-1), N3(g), N4(g), N7(g-1), N8(g-1), N7(g), N8(g), &
                                       blkm%M(Im, Jm+g, Km), blkm%dl(Im, Jm+g, Km), blkm%vol(Im, Jm+g, Km))
          Mg(g)   = blkm % M(Im, Jm+g, Km)
          dlg(g)  = blkm % dl(Im, Jm+g, Km)
          volg(g) = blkm % vol(Im, Jm+g, Km)
        case(5)
          call Compute_Metric_Tensor ( N1(g), N1(g-1), N3(g), N3(g-1), N5(g), N5(g-1), N7(g), N7(g-1), &
                                       blkm%M(Im, Jm, Km-g), blkm%dl(Im, Jm, Km-g), blkm%vol(Im, Jm, Km-g))
          Mg(g)   = blkm % M(Im, Jm, Km-g)
          dlg(g)  = blkm % dl(Im, Jm, Km-g)
          volg(g) = blkm % vol(Im, Jm, Km-g)
        case(6)
          call Compute_Metric_Tensor ( N2(g-1), N2(g), N4(g-1), N4(g), N6(g-1), N6(g), N8(g-1), N8(g), &
                                       blkm%M(Im, Jm, Km+g), blkm%dl(Im, Jm, Km+g), blkm%vol(Im, Jm, Km+g))
          Mg(g)   = blkm % M(Im, Jm, Km+g)
          dlg(g)  = blkm % dl(Im, Jm, Km+g)
          volg(g) = blkm % vol(Im, Jm, Km+g)
        end select

    enddo ! End loop over ghost cell           

  end subroutine BC_Connect_Metrics


  !> The walls of the wall distance: coupled to another solver (103, 104), heat flux (301), temperature (302).
  function Is_Wall(type) result (ans)
      implicit none
      integer, intent(in) :: type
      logical :: ans

      select case (type)
        case (103, 104, 301, 302)  ! wall BCs
          ans = .true.
        case default
          ans = .false.
      end select
    end function Is_Wall

end module ARES_Lib_Metrics
