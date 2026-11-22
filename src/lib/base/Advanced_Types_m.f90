module ARES_Advanced_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Base_Types_m
  use ARES_Parameters_m
  use Lib_ORION_Data
  
  implicit none
  
  ! FOUNDAMENTAL TYPES
  type :: block_type
    integer                                :: dim(3)           ! Number of cells in i-j-k (ghost not included)
    real(R8), allocatable                  :: vol(:,:,:)       ! Cell volume
    type(ARES_vector_3D_type), allocatable :: node(:,:,:)      ! Mesh grid points (including ghost)
    type(ARES_tensor_3D_type), allocatable :: M(:,:,:)         ! Metric transformation tensor
    type(ARES_tensor_3D_type), allocatable :: Ms(:,:,:)        ! Cell metric of the RANS source terms, diffusive-metric = stencil (velocity)
    type(ARES_tensor_3D_type), allocatable :: Msc(:,:,:)       ! the same for the density
    type(ARES_tensor_3D_type), allocatable :: Msa(:,:,:)       ! the same for the velocity gradients of SA-RC, with their ghost cells
    type(ARES_tensor_3D_type), allocatable :: Mst(:,:,:)       ! the same for the turbulence variables, where it is not Ms
    type(ARES_vector_3D_type), allocatable :: dl(:,:,:)        ! Average cell length (in i/j/k direction). eg: dl%c(1) is sqrt(dx**2+dy**2+dz**2) of the cell in the i direction
    type(d_metrics_type)                   :: dir(3)           ! Direction object. Contains: i-faces, j-faces, k-faces; eg: dir(1)%face(i,j,k)%n
    real(R8), allocatable                  :: yn(:,:,:)        ! Nearest wall distance
    real(R8), allocatable                  :: ks(:,:,:)        ! Roughness of nearest wall element
  end type block_type

  ! BOUNDARY CONDITIONS TYPES
  type :: bc_type
    integer                    :: i, j, k, b, f                             ! ijk coordinates, block and face in which the boundary element is located
    integer                    :: type                                      ! BC type (1,9)
    integer                    :: bs, is, js, ks, fs, d11, d12, d21, d22    ! BC 1 (connection) specifications
    type(ARES_tensor_3D_type)  :: Mg(2)                                     ! Ghost cell metric tensor
    type(ARES_vector_3D_type)  :: dlg(2)                                    ! Ghost cell average cell length
    real(R8)                   :: volg(2)                                   ! Ghost cell volume
    real(R8), allocatable      :: Pg(:,:)                                   ! Ghost cell primitive stencil
    ! BC 101/201: record of the same face seen from the other block
    integer                    :: partner = 0
    logical                    :: owner = .false.                           ! BC 101/201: this record computes the flux of the face
    ! BC 101/201: residual of the face for this record's cell, (:,1) inviscid, (:,2) viscous
    real(R8), allocatable      :: Fc(:,:)
    ! BC 200 (ARES_Lib_BC_Fluxes_Rotational): face metric of the viscous flux, set once (Setup_Metrics)
    real(R8)                   :: rM(3,3) = 0.0_R8                          ! Rotational_Periodic_Face_Metric
    integer                    :: ni(2)                                     ! BC chimera
    integer, allocatable       :: donorID(:,:)                              ! BC chimera
    real(R8), allocatable      :: volume_fraction(:)                        ! BC chimera
    real(R8), allocatable      :: ext_flux(:)                               ! Multi-Solver Coupling
  end type bc_type


  ! REAL GAS EXTENSION 
  type, extends(block_type) :: ARES_block_type
    real(R8), dimension(:,:,:,:), allocatable :: P, PO                 ! Primitive variables at time n and n-1: { p vel h }
    real(R8), dimension(:,:,:,:), allocatable :: R                     ! Residuals
    real(R8), dimension(:,:,:,:), allocatable :: Visc                  ! Viscous variables of the cells, once per stage
    real(R8), dimension(:,:,:,:), allocatable :: TE                    ! Truncation error
    real(R8), dimension(:,:,:,:), allocatable :: RC                    ! Conservation form residuals (temporary storage)
    real(R8), dimension(:,:,:,:), allocatable :: RS1, RS2              ! Implicit smoothing residuals (temporary storage)
    real(R8), dimension(:,:,:), allocatable   :: dtlocal               ! Local time step
    type(ARES_tensor_3D_type), allocatable    :: vel_gradient(:,:,:)   ! Gradient of velocity
    real(R8), allocatable                     :: rc_term1(:,:,:)       ! Spalart-Shur rotation/curvature correction terms
    real(R8), allocatable                     :: rc_term2(:,:,:)       ! Spalart-Shur rotation/curvature correction terms
  end type ARES_block_type


  type, extends(bc_type) :: ARES_bc_type
    real(R8)                  :: qw, Tw, Taw, hg, qrad             ! BC wall specifications 
    real(R8)                  :: eps_wall = 0.0_R8, k_rough = 0.0_R8 ! BC viscous wall specifications   
    real(R8)                  :: sx = 0.0_R8                       ! Walls and 400: l0/l1 of the MUSCL distances (Extrapolation_Ratio)
    ! Wall faces with diffusive-metric = stencil (ARES_Lib_Stencil_Metrics: Wall_Face_Stencils)
    logical                   :: wstencil = .false.                ! wM, wdT set for this wall face
    real(R8)                  :: wM(3,3) = 0.0_R8                  ! face metric of the wall gradient
    real(R8)                  :: wdT(2) = 0.0_R8                   ! tangential differences of T_w (302)
    real(R8)                  :: T0, p0, alpha, beta, mach, pamb   ! BC 4 (inflow/outflow) specifications
    real(R8)                  :: mdot                              ! BC 4 (inflow/outflow) specifications
    real(R8)                  :: rel_fac                           ! BC 4 (inflow/outflow) specifications
    real(R8), allocatable     :: ci(:)                             ! BC 4 (inflow/outflow) specifications
    ! BC 501 (manifold): the group of the record, its role (1 outlet, 2 inlet), its index in domain % man, and the
    ! parameters of the group (rel_fac above): the loss K on the dynamic pressure of the outlet faces, the mass flow,
    ! total pressure and total enthalpy of the inlet faces against those of the outlet faces (Lib_BC_Fluxes_Manifold)
    integer                   :: man_gid = 0, man_role = 0, man_idx = 0
    real(R8)                  :: man_K = 0.0_R8, man_mratio = 1.0_R8, man_p0ratio = 1.0_R8, man_dh0 = 0.0_R8
  end type ARES_bc_type


  !> A manifold (BC 501): one group of faces joined by a lumped plenum, built from the records of a grid level
  !> (Manifold_Setup) and given the state its faces take at each stage (Manifold_Gather, Lib_BC_Fluxes_Manifold)
  type :: ARES_manifold_type
    integer  :: gid = 0                                  ! the group id of the records
    integer  :: n_out = 0, n_in = 0                      ! records of each role in the whole domain
    integer  :: first = 0                                ! the first record of the group (its parameters)
    real(R8) :: K = 0.0_R8, rf = 1.0_R8                  ! loss on the dynamic pressure, relaxation factor: outlet faces
    real(R8) :: mratio = 1.0_R8, p0ratio = 1.0_R8        ! mass flow and total pressure of the inlet faces over the outlet's
    real(R8) :: dh0 = 0.0_R8                             ! total enthalpy taken out between them [J/kg]
    logical  :: active = .false.                         ! gas leaves through the outlet faces at this stage
    real(R8) :: g_in = 0.0_R8, h0_in = 0.0_R8            ! mass flux and total enthalpy the inlet faces take
    real(R8) :: p_out = 0.0_R8                           ! static pressure the outlet faces take
    real(R8), allocatable :: ci_in(:)                    ! turbulence (1:nrans, specific) of the inlet faces
    real(R8) :: mdot_out = 0.0_R8, mdot_in = 0.0_R8      ! the line of the log: mass flow, total enthalpy and total pressure
    real(R8) :: h0_out = 0.0_R8                          !   of the two sides (h0_in above)
    real(R8) :: p0_out = 0.0_R8, p0_in = 0.0_R8
    integer  :: last_line = -1                           ! the iteration of the last line written
  end type ARES_manifold_type


  ! COMPOUND TYPES
  type :: ARES_domain_type
    real(R8)                                         :: time             ! Solution time
    real(R8)                                         :: dtglobal         ! Global dt for time accurate simulation
    integer                                          :: iter, itermax    ! Iteration number
    integer                                          :: nb, nbound       ! Number of blocks, number of boundary faces
    integer, dimension(:,:), allocatable             :: n_bf             ! Number of bc elements per faces per block
    type(ARES_block_type), dimension(:), allocatable :: blk              ! Allocatable block type
    type(ARES_bc_type), dimension(:), allocatable    :: bc               ! Allocatable bc object
    ! MPI local BC indices
    integer                                          :: n_local_bc = 0
    integer, dimension(:), allocatable               :: local_bc_idx
    integer                                          :: n_local_bs = 0
    integer, dimension(:), allocatable               :: local_bs_idx
    ! The same local BC indices grouped by block and face, for BC_Fluxes: the records of face f
    ! of block b are local_bf_idx(local_bf_first(b,f) : local_bf_first(b,f) + local_bf_count(b,f) - 1)
    integer, dimension(:,:), allocatable             :: local_bf_first
    integer, dimension(:,:), allocatable             :: local_bf_count   ! 0 on blocks of other ranks
    integer, dimension(:), allocatable               :: local_bf_idx
    integer                                          :: mg_level = 1     ! Multigrid level of this domain (selects exchange schedules)
    ! The manifolds (501) of this grid level (Manifold_Setup)
    integer                                          :: nman = 0
    type(ARES_manifold_type), dimension(:), allocatable :: man
  end type ARES_domain_type


  type :: ARES_simulation_type
    type(ARES_domain_type), dimension(:), allocatable  :: domain
    type(ORION_data), dimension(:), allocatable        :: IOfield
  end type ARES_simulation_type


  abstract interface
    !> The procedure that the caller of ARES%solve passes to it (the program of ARES and the couplings of hydra pass an
    !> empty one): called once per stage by every thread of the team, inside the parallel region of Explicit_Step, after
    !> the fluxes and before the turbulence sources. It takes no argument; a loop over cells in it needs a worksharing
    !> construct, anything done once a single or master construct
    subroutine ARES_external_function ()
    end subroutine ARES_external_function
  end interface


end module ARES_Advanced_Types_m