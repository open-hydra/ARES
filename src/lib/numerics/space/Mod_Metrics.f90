module ARES_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

contains

  subroutine Setup_Metrics ( domain )
    use ARES_Base_Types_m
    use ARES_Advanced_Types_m
    use ARES_Global_m
    use ARES_Lib_Metrics
    use ARES_Mod_MPI, only: is_local_block
    implicit none
    type(ARES_domain_type), intent(inout) :: domain
    ! Local
    integer :: b, i, j, k
    integer :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s
    type(ARES_vector_3D_type) :: N1, N2, N3, N4, N5, N6, N7, N8

    call Check_Mesh_Type ( domain )

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
        call Compute_Metric_Tensor ( N1, N2, N3, N4, N5, N6, N7, N8, domain % blk(b) % M(i,j,k), domain % blk(b) % dl(i,j,k), domain % blk(b) % vol(i,j,k) )
      enddo; enddo; enddo
      !$omp end parallel

      ! Compute Normal & Area
      call Compute_Norm_Area ( domain % blk(b) )
    enddo

    if (model>1) then

      ! Compute wall minimum distance => block%yn
      call Compute_Yn ( domain % blk(:), domain % bc, domain % nb, domain % nbound ) 
      
      !$omp parallel
      ! Set wall-distance in connected cells
      !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks)
      do i = 1, domain % nbound
        select case ( domain % bc(i) % type )
          case(101) ! block connection
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
        end select
      enddo
      !$omp end parallel
    endif
    
    !$omp parallel
    ! Create the nodes for gc layers of ghost cell
    !$omp do schedule (dynamic) private(i, Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs, d11s, d12s, d21s, d22s)
    do i = 1, domain % nbound
      Bm = domain % bc(i) % b
      Im = domain % bc(i) % i 
      Jm = domain % bc(i) % j 
      Km = domain % bc(i) % k 
      Fm = domain % bc(i) % f
      select case ( domain % bc(i) % type )
        case(101) ! block connection
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
                                    domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg)
        case(300) ! symmetry
          call BC_Symmetry_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                     domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg )
        case default  ! Rientra qui anche chimera
          call BC_Extrapolate_Metrics ( Im, Jm, Km, Fm, domain % blk(Bm), &
                                        domain % bc(i) % Mg, domain % bc(i) % dlg, domain % bc(i) % volg )
      end select
    enddo
    !$omp end parallel

  end subroutine Setup_Metrics

end module ARES_Mod_Metrics