module ARES_Mod_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Fluxes

contains

    subroutine Fluxes ( domain )
      use ARES_Advanced_Types_m
      use ARES_Mod_MPI, only: is_local_block
      use ARES_Config_Types_m, only: obj_riemann, obj_rans

      implicit none
      type(ARES_domain_type), intent(inout) :: domain
      ! Local
      integer :: b
      real(R8) :: Sc, Sct, Prt
      logical :: Prt_corr

      Sc  = obj_rans%Sc
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