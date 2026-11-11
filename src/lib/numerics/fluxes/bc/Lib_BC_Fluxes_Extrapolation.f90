module ARES_Lib_BC_Fluxes_Extrapolation
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use ARES_Advanced_Types_m
  use ARES_Global_m
  use ARES_Parameters_m
  use ARES_Lib_BC_Fluxes, only: Face_Index

  implicit none
  public

  !> Flux evaluations of 400 faces where the gas of the boundary cell moves into the domain, and the largest -w1/a1 among them
  !> (w1 the velocity of that gas along the normal out of the domain, a1 its speed of sound): reported at the end of a run in
  !> which the gas entered (Report_400_Entering)
  integer(I8) :: n_400_entering = 0
  real(R8)    :: w_400_entering = 0d0
  !> The lines of the gas entering through each 400 boundary (the 400 faces of face f of block b, on each grid level), written by
  !> Report_400_Entering_Lines: x_400_line(f,b,level) is the largest -w1/a1 of the last line of the boundary, zero before its
  !> first one; n_400_line(f,b) counts the flux evaluations of a stage on the boundary with -w1/a1 above 10 times x_400_line,
  !> w_400_line(f,b) is their largest -w1/a1 and count_400_line their sum over the boundaries (Setup_400_Entering allocates them)
  integer,  allocatable :: n_400_line(:,:)
  real(R8), allocatable :: w_400_line(:,:)
  real(R8), allocatable :: x_400_line(:,:,:)
  integer :: count_400_line = 0

contains

  subroutine BC_Extrapolation ( Bm, Im, Jm, Km, Fm, Blk )
    use FLINT_Lib_Thermodynamic

    implicit none
    integer, intent(in) :: Bm, Im, Jm, Km, Fm
    type(ARES_block_type), intent(inout) :: Blk
    ! Local
    integer :: Dir, Face_i, Face_j, Face_k, FluxSign
    real(R8) :: Normal(3), Area, Prim(nprim), Un, Fmass, Flux(nprim), rho
    real(R8) :: Sound_c


    ! Boundary face coordinates
    call Face_Index ( Fm, dir, Im, Jm, Km, Face_i, Face_j, Face_k )

    ! Metric stuff
    Normal = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % n
    Area = Blk % dir(Dir) % f(Face_i,Face_j,Face_k) % a

    ! boundary cell primitive variables
    Prim = Blk % P(:,Im,Jm,Km)
    rho = ph2vars( Prim(np), Prim(nh), rho_tab )
    ! boundary cell: normal velocity to the boundary face
    Un = dot_product ( Prim(nu:nw), Normal )

    Flux(np) = rho * Un * Area
    Fmass = rho * Un * Area
    Flux(nu:nw) =  Fmass * Prim(nu:nw) + Prim(np) * Area * Normal
    Flux(nh) = Fmass * ( Prim(nh) + 0.5_R8* norm2 ( Prim(nu:nw) )**2 )
    if (model==2) Flux(nt:nprim) = Prim(nt:nprim) * Un * Area

    ! Residual update
    FluxSign = 1 - 2 * mod(Fm,2)
    Blk % r(:,Im,Jm,Km) = Blk % r(:,Im,Jm,Km) + FluxSign * Flux

    ! Where the gas of the boundary cell moves into the domain, w1 = FluxSign Un < 0 (its velocity along the normal out of the
    ! domain), it enters with the state of the cell: the evaluation is counted, with -w1/a1, for the lines of the run
    if ( FluxSign * Un < 0d0 ) then
      Sound_c = ph2vars( Prim(np), Prim(nh), sound_tab )
      call Note_400_Entering ( Bm, Fm, -FluxSign * Un / Sound_c )
    end if

  end subroutine BC_Extrapolation


  !> One flux evaluation of a 400 face of face Fm of block Bm where the gas of the boundary cell moves into the domain, x = -w1/a1:
  !> counted for the end line (Report_400_Entering) and, when x is more than 10 times the largest -w1/a1 of the last line of the
  !> boundary on the current grid level (zero before its first line), for its next line (Report_400_Entering_Lines). Called
  !> inside the worksharing loops of the boundary fluxes: the updates are atomic, and x_400_line changes only between two
  !> barriers that no flux evaluation crosses
  subroutine Note_400_Entering ( Bm, Fm, x )
    use ARES_Config_Types_m, only: obj_multigrid
    implicit none
    integer,  intent(in) :: Bm, Fm
    real(R8), intent(in) :: x

    !$omp atomic
    n_400_entering = n_400_entering + 1
    !$omp atomic
    w_400_entering = max ( w_400_entering, x )
    if ( x > 10d0 * x_400_line(Fm,Bm,obj_multigrid%MG_level) ) then
      !$omp atomic
      n_400_line(Fm,Bm) = n_400_line(Fm,Bm) + 1
      !$omp atomic
      w_400_line(Fm,Bm) = max ( w_400_line(Fm,Bm), x )
      !$omp atomic
      count_400_line = count_400_line + 1
    end if

  end subroutine Note_400_Entering


  !> One line for each 400 boundary through which the gas enters at this stage with -w1/a1 above 10 times the largest of its last
  !> line, on each grid level (the first line of a boundary comes at its first entry): the block, the face (the grid level), the
  !> iteration (the Iter of the lines on the screen), the faces counted at this stage (all those with gas entering in the first
  !> line, those above 10 times the last line in the next ones) and their largest -w1/a1, which is the largest of the boundary at
  !> this stage and becomes that of its last line. Every thread of the team calls it after the boundary fluxes of the stage, past
  !> the barrier that ends the last worksharing loop of their evaluations, where count_400_line is the same for every thread: in
  !> BC_Fluxes_core the call stays after the last of its !$omp do loops over the faces of the blocks, which have no nowait
  !> (Connection_Fluxes runs before them); with a nowait there, or the call before that barrier, a thread could read
  !> count_400_line while another one is adding to it, and the threads would not meet the same barriers. The rank that owns the
  !> block writes the line and flushes it before any thread goes on to the update of the cells, so that a run that then stops on
  !> an error has it in its output (such a run prints no end line). No MPI call
  subroutine Report_400_Entering_Lines ()
    use iso_fortran_env,     only: output_unit
    use ARES_Config_Types_m, only: obj_multigrid, obj_sim_param
    implicit none
    integer :: b, f, lev
    character(len=:), allocatable :: faces

    if ( count_400_line == 0 ) return
    ! Every thread has read count_400_line before the single region below sets it to zero
    !$omp barrier
    !$omp single
    lev = obj_multigrid%MG_level
    do b = 1, size ( n_400_line, 2 )
      do f = 1, 6
        if ( n_400_line(f,b) == 0 ) cycle
        if ( x_400_line(f,b,lev) > 0d0 ) then
          faces = ', faces with |w|/a above 10 times its last line: '
        else
          faces = ', faces: '
        end if
        if ( obj_multigrid%MGL > 1 ) then
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES12.5)') ' 400 boundary of block ', b, ', face ', f, ', grid level ', lev, &
            ', gas entering at iteration ', obj_sim_param%iter_from_call, faces, n_400_line(f,b), ', largest |w|/a: ', &
            w_400_line(f,b)
        else
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES12.5)') ' 400 boundary of block ', b, ', face ', f, &
            ', gas entering at iteration ', obj_sim_param%iter_from_call, faces, n_400_line(f,b), ', largest |w|/a: ', &
            w_400_line(f,b)
        end if
        x_400_line(f,b,lev) = w_400_line(f,b)
        n_400_line(f,b) = 0
        w_400_line(f,b) = 0d0
      end do
    end do
    count_400_line = 0
    flush ( output_unit )
    !$omp end single

  end subroutine Report_400_Entering_Lines


  !> The record of the lines of the gas entering (Report_400_Entering_Lines) for the blocks of the domain on its grid levels, and
  !> the counts of the end line (Report_400_Entering) set to zero: every rank, at set-up (the levels have the same blocks)
  subroutine Setup_400_Entering ( domain, levels )
    implicit none
    integer, intent(in) :: levels
    type(ARES_domain_type), intent(in) :: domain(levels)
    ! Local
    integer :: nb

    nb = domain(1) % nb
    if ( allocated ( n_400_line ) ) deallocate ( n_400_line, w_400_line, x_400_line )
    allocate ( n_400_line(6,nb), w_400_line(6,nb), x_400_line(6,nb,levels) )
    n_400_line = 0
    w_400_line = 0d0
    x_400_line = 0d0
    count_400_line = 0
    n_400_entering = 0
    w_400_entering = 0d0

  end subroutine Setup_400_Entering


  !> One line at the end of a run in which the gas entered through 400 faces: the flux evaluations of 400 faces where the gas of
  !> the boundary cell moved into the domain, over every stage, grid level and rank, and the largest -w1/a1 among them: about 1e-9
  !> or less where only round-off moves the gas towards the domain, 0.01 and above where the gas really enters,
  !> docs/user/boundary-conditions.md. A run in which no gas entered prints no line. Collective: every rank calls it, and both
  !> reductions, in every run, with or without 400 faces, so that no rank takes another MPI path; rank 0 writes the line when the
  !> total of the evaluations is not zero (only BC_Extrapolation counts them, for the faces of type 400)
  subroutine Report_400_Entering ()
    use ARES_Mod_MPI, only: mpi_is_root, mpi_allreduce_sum_r8, mpi_allreduce_max_r8
    implicit none
    real(R8) :: n, x

    call mpi_allreduce_sum_r8 ( real(n_400_entering, R8), n )
    call mpi_allreduce_max_r8 ( w_400_entering, x )
    if ( mpi_is_root .and. n > 0d0 ) write(*,'(A,I0,A,ES12.5)') ' 400 faces with gas entering, flux evaluations: ', &
                                                               nint(n, I8), ', largest |w|/a: ', x

  end subroutine Report_400_Entering

end module ARES_Lib_BC_Fluxes_Extrapolation
