module ARES_Mod_Multigrid
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Setup_Multigrid, Restriction, Prolongation

contains

  subroutine Setup_Multigrid ( simulation )
    use ARES_Advanced_Types_m
    use ARES_Config_Types_m, only: obj_multigrid
    use ARES_Lib_Multigrid, only: Check_Multigrid, Coarse_Grid, Coarse_IOfield
    use ARES_Mod_Allocate_Data, only: Allocate_Block
    implicit none
    type(ARES_simulation_type), intent(inout) :: simulation
    ! Local
    integer :: b, m, i, j, k, nb, rap

    call Check_Multigrid ( simulation%domain(1) )

    nb = simulation%domain(1) % nb

    do m = 2, obj_multigrid%MGL
      
      allocate ( simulation%domain(m) % Blk(nb) )
      simulation%domain(m) % nb = nb
      simulation%domain(m) % time = simulation%domain(1) % time

      rap = 2
      do b = 1, nb
        simulation%domain(m) % Blk(b) % Dim = simulation%domain(m-1) % Blk(b) % Dim / rap
        simulation%domain(m) % Blk(b) % Dim(3) = Max ( 1, simulation%domain(m) % Blk(b) % Dim(3) ) ! 2D case
        call Allocate_Block ( simulation%domain(m) % Blk(b), simulation%domain(m) % Blk(b) % Dim )
      enddo

      call Coarse_Grid ( simulation%domain(m-1), simulation%domain(m) )

      call Coarse_IOfield ( simulation%IOfield(m-1), simulation%IOfield(m) )

      do b = 1, nb
        do k = 0, simulation%domain(m)%Blk(b)%Dim(3)
        do j = 0, simulation%domain(m)%Blk(b)%Dim(2)
        do i = 0, simulation%domain(m)%Blk(b)%Dim(1)
          simulation%IOfield(m)%block(b)%mesh(:,i,j,k) = simulation%domain(m) % Blk(b) % Node(i,j,k) % c
        enddo; enddo; enddo
      enddo

    enddo

    ! Read number of iterations for each domain
    ! (HYDRA coupling provides them through the HYDRA-Multigrid section)
    do m = 1, obj_multigrid%MGL
      simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
    enddo

  end subroutine Setup_Multigrid


  subroutine Restriction ( Fine, Coarse )
    use ARES_Advanced_Types_m
    use ARES_Lib_Multigrid, only: fine2coarse_prim
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b

    do b = 1, Coarse % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      call fine2coarse_prim ( Fine % Blk(b) % P,   Coarse % Blk(b) % P  , &
                              Fine % Blk(b) % vol, Coarse % Blk(b) % vol, &
                              Fine % Blk(b) % dim, Coarse % Blk(b) % dim )
    enddo

  end subroutine Restriction


  subroutine Prolongation ( Fine, Coarse )
    use ARES_Advanced_Types_m
    use ARES_Lib_Multigrid, only: coarse2fine_prim
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: Fine, Coarse
    ! Local
    integer :: b

    do b = 1, Coarse % nb ! Loop over blocks
      if (.not. is_local_block(b)) cycle
      call coarse2fine_prim ( Fine % Blk(b) % P, Coarse % Blk(b) % P, &
                              Fine % Blk(b) % dim, Coarse % Blk(b) % dim )
    enddo

  end subroutine Prolongation

end module ARES_Mod_Multigrid