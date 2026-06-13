module ARES_Mod_dt
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_dt, Set_Global_dt

contains

  subroutine Compute_dt ( domain, cfl, vnn, rampa_iter )
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Lib_RANS
    use ARES_Mod_MPI, only: is_local_block, mpi_allreduce_min_r8
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: cfl, vnn
    integer, intent(in)  :: rampa_iter
    ! Local
    integer  :: i, j, k, b
    real(R8) :: dtcell, dtglobal, dtglobal_mpi

    dtglobal = domain % dtglobal

    do b = 1, domain % nb ! loop over blocks
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3) private ( dtcell ), reduction ( min : dtglobal )
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)

        ! Compute local cell dt according to CFL and VNN numbers
        call compute ( pressure = domain%blk(b)%p(np,i,j,k), &
                       vel = domain%blk(b)%p(nu:nw,i,j,k), &
                       entalpy = domain%blk(b)%p(nh,i,j,k), &
                       rans_ = domain%blk(b)%P(nt:nprim,i,j,k), &
                       met = domain%blk(b)%m(i,j,k)%c, &
                       dl = domain%blk(b)%dl(i,j,k)%c, &
                       dtmin = dtcell, &
                       cfl = cfl, &
                       vnn = vnn )

        ! Apply CFL reduction if required
        if ( domain%iter < rampa_iter ) dtcell = dtcell * domain%iter / rampa_iter

        ! Update local cell dt and global minimum dt
        domain % blk(b) % dtlocal(i,j,k) = dtcell
        dtglobal = min ( dtcell, dtglobal )

      enddo; enddo; enddo
      !$omp end parallel
    enddo ! end of loop over blocks

    ! MPI: global minimum across all ranks
    call mpi_allreduce_min_r8(dtglobal, dtglobal_mpi)
    domain % dtglobal = dtglobal_mpi

    contains
      
      subroutine compute ( pressure, vel, entalpy, rans_, met, dl, dtmin, cfl, vnn )
        use ARES_Global_m
        use FLINT_Lib_Thermodynamic
        use ARES_Lib_Preconditioning
        use ARES_Config_Types_m, only: obj_prec

        implicit none
        real(R8), intent(in)  :: pressure, vel(3), entalpy, rans_(:)
        real(R8), intent(in)  :: met(3,3), dl(3), cfl, vnn
        real(R8), intent(out) :: dtmin
        ! Local
        integer :: d
        real(R8) :: rho, Rgas, Sound, dt, versor(3), lambda, mie, mil, mi, dx
        real(R8) :: Alpha, Beta, vel_, U_r, dummy(3,3)=0d0

        Sound = ph2vars( pressure, entalpy, sound_tab )
        rho   = ph2vars( pressure, entalpy, rho_tab   )
        if ( obj_prec%enabled ) then
          dx = min( dl(1), dl(2), dl(3))
          U_r = comp_Ur ( vel=Norm2( Vel ), sound=Sound ) 
          Beta = 1.d0/Sound**2
          if ( Uref > 0.d0 ) Beta = 1.d0/Uref**2 ! in case of reference sound velocity
          Alpha = 5.d-1 * ( 1.d0 - Beta * U_r**2 )
        endif

        dtmin = 1d8

        do d = 1, ndir
        
          ! CFL condition along d-direction
          versor = met(d,:) / norm2 ( met(d,:) )
          vel_ = abs ( dot_product ( vel, versor ) )

          if ( obj_prec%enabled ) then
            Sound = Sqrt ( Alpha**2 * Vel_**2 + U_r**2 )
            vel_ = vel_ * ( 1.d0 - Alpha )
          endif

          lambda  = vel_ + Sound
          dt = dl(d) / lambda * cfl
          dtmin = min ( dt, dtmin )
        
        enddo
      
        if ( model==0 ) return
        
        mie = 0.d0
        mil = ph2vars( pressure, entalpy, mi_tab2D )
        if (model==2) then 
          ! Note: approximate mit for 2 equation models. 
          ! velocity gradient assumed 0 and small distance from wall
          call Eddy_Viscosity ( mut=mie, rans_variables=rans_, mul=mil, rho=rho, &
                                vel_gradient=dummy, walldist=1d-6, ks=0d0 )
        endif
        mi = mie + mil

        do d = 1, ndir

          ! VNN condition along d-direction
          dt = ( rho * dl(d)**2 * vnn ) / mi
          dtmin = min ( dt, dtmin )

        end do

      end subroutine compute

  end subroutine Compute_dt


  subroutine Set_Global_dt ( domain )
    use ARES_Advanced_Types_m
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k

    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      !$omp parallel
      !$omp do collapse(3)
      do k = 1, domain % blk(b) % dim(3)
      do j = 1, domain % blk(b) % dim(2)
      do i = 1, domain % blk(b) % dim(1)
        domain % blk(b) % dtlocal(i,j,k) = domain % dtglobal
      enddo; enddo; enddo
      !$omp end parallel
    enddo
    
  end subroutine Set_Global_dt

end module ARES_Mod_dt