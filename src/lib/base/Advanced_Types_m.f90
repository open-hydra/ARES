module ARES_Advanced_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ARES_Base_Types_m
  use ARES_Parameters_m
  use ARES_Series_Data_m
  use Lib_ORION_Data
  
  implicit none
  
  ! FOUNDAMENTAL TYPES
  type :: block_type
    integer                                :: dim(3)           ! Number of cells in i-j-k (ghost not included)
    real(R8), allocatable                  :: vol(:,:,:)       ! Cell volume
    type(ARES_vector_3D_type), allocatable :: node(:,:,:)      ! Mesh grid points (including ghost)
    type(ARES_tensor_3D_type), allocatable :: M(:,:,:)         ! Metric transformation tensor
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
    !integer                    :: ni(2)                                    ! BC chimera
    !integer, allocatable       :: donorID(:,:)                             ! BC chimera
    !real(R8), allocatable      :: volume_fraction(:)                       ! BC chimera
    real(R8), allocatable      :: ext_flux(:)                               ! Multi-Solver Coupling
  end type bc_type


  ! REAL GAS EXTENSION 
  type, extends(block_type) :: ARES_block_type
    real(R8), dimension(:,:,:,:), allocatable :: P, PO                 ! Primitive variables at time n and n-1: { p vel h }
    real(R8), dimension(:,:,:,:), allocatable :: R                     ! Residuals
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
    real(R8)                  :: eps_wall, k_rough = 0.0_R8        ! BC viscous wall specifications   
    real(R8)                  :: T0, p0, alpha, beta, mach, pamb   ! BC 4 (inflow/outflow) specifications
    real(R8)                  :: mdot                              ! BC 4 (inflow/outflow) specifications
    real(R8)                  :: rel_fac                           ! BC 4 (inflow/outflow) specifications
    real(R8), allocatable     :: ci(:)                             ! BC 4 (inflow/outflow) specifications
    type(time_series_type)    :: p0time         
    logical                   :: q2d_periodic = .false.            ! BC 667 (Q2D mapped) periodic time signal flag
  end type ARES_bc_type


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
  end type ARES_domain_type


  type :: ARES_simulation_type
    type(ARES_domain_type), dimension(:), allocatable  :: domain
    type(ORION_data), dimension(:), allocatable        :: IOfield
  end type ARES_simulation_type


end module ARES_Advanced_Types_m