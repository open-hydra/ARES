!>@brief Wall-clock instrumentation of the solver loop.
!> Every rank accumulates the time spent in one `Explicit_Step`, the time
!> blocked on the ghost-cell exchange and the time blocked on collectives.
!> The accumulators are reduced over the ranks and reported every `timer-diter`
!> iterations (0 = off) and once at the end of the run.
!>
!> Times are reported as the **maximum** over the ranks, which is the critical
!> path and therefore what sets the time to solution; the spread between the
!> slowest and the mean rank is the load imbalance.
module ARES_Mod_Timers
#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env, only: R8 => real64, I8 => int64

  implicit none
  private

  !> Reduced timings for one window of iterations.
  type :: stats_type
    real(R8) :: tmax = 0.0_R8       !< iteration time on the slowest rank
    real(R8) :: tmin = 0.0_R8       !< iteration time on the fastest rank
    real(R8) :: tavg = 0.0_R8       !< iteration time, mean over ranks
    real(R8) :: imbalance = 0.0_R8  !< (tmax-tavg)/tavg, %
    real(R8) :: commfrac = 0.0_R8   !< time blocked on the halo exchange, %
    real(R8) :: syncfrac = 0.0_R8   !< time blocked on collectives, %
    real(R8) :: wmax = 0.0_R8       !< compute time (iteration minus waits), max
    real(R8) :: wmin = 0.0_R8       !< compute time, min
    real(R8) :: wavg = 0.0_R8       !< compute time, mean
    real(R8) :: wspread = 0.0_R8    !< (wmax-wavg)/wavg, %
    integer  :: rmax = 0            !< rank holding wmax
    integer  :: rmin = 0            !< rank holding wmin
  end type stats_type

  real(R8) :: t_iter_beg = 0.0_R8   !< start of the current iteration
  real(R8) :: t_comm_beg = 0.0_R8   !< start of the current exchange wait
  real(R8) :: t_sync_beg = 0.0_R8   !< start of the current collective
  real(R8) :: t_iter_acc = 0.0_R8   !< iteration time since the last report
  real(R8) :: t_comm_acc = 0.0_R8   !< exchange wait since the last report
  real(R8) :: t_sync_acc = 0.0_R8   !< collective wait since the last report
  integer  :: n_iter_acc = 0        !< iterations since the last report

  real(R8) :: t_run_beg = 0.0_R8    !< start of the solver loop
  real(R8) :: t_run_iter = 0.0_R8   !< iteration time over the whole run
  real(R8) :: t_run_comm = 0.0_R8   !< exchange wait over the whole run
  real(R8) :: t_run_sync = 0.0_R8   !< collective wait over the whole run
  integer  :: n_run = 0             !< iterations over the whole run

  public :: timer_wtime
  public :: timer_run_begin, timer_summary
  public :: timer_iter_begin, timer_iter_end
  public :: timer_comm_begin, timer_comm_end
  public :: timer_sync_begin, timer_sync_end
  public :: timer_report

contains


  !> Wall-clock reading, in seconds. MPI_WTIME when available (monotonic and
  !> consistent across the ranks of a job), SYSTEM_CLOCK otherwise.
  real(R8) function timer_wtime()
#ifdef USE_MPI
    timer_wtime = MPI_WTIME()
#else
    integer(I8) :: count, rate
    call system_clock(count, rate)
    timer_wtime = real(count, R8) / real(rate, R8)
#endif
  end function timer_wtime


  !> Mark the start of the solver loop, i.e. the end of set-up.
  !> Clears the window accumulators as well as the run totals.
  subroutine timer_run_begin()
    t_run_beg  = timer_wtime()
    t_run_iter = 0.0_R8
    t_run_comm = 0.0_R8
    t_run_sync = 0.0_R8
    n_run      = 0

    t_iter_acc = 0.0_R8
    t_comm_acc = 0.0_R8
    t_sync_acc = 0.0_R8
    n_iter_acc = 0
  end subroutine timer_run_begin


  !> Open the iteration timer. Call outside any OpenMP parallel region.
  subroutine timer_iter_begin()
    t_iter_beg = timer_wtime()
  end subroutine timer_iter_begin


  !> Close the iteration timer and accumulate.
  subroutine timer_iter_end()
    real(R8) :: dt
    dt = timer_wtime() - t_iter_beg
    t_iter_acc = t_iter_acc + dt
    t_run_iter = t_run_iter + dt
    n_iter_acc = n_iter_acc + 1
    n_run      = n_run + 1
  end subroutine timer_iter_end


  !> Open the exchange-wait timer. Call from a single-threaded region only.
  subroutine timer_comm_begin()
    t_comm_beg = timer_wtime()
  end subroutine timer_comm_begin


  !> Close the exchange-wait timer and accumulate.
  subroutine timer_comm_end()
    real(R8) :: dt
    dt = timer_wtime() - t_comm_beg
    t_comm_acc = t_comm_acc + dt
    t_run_comm = t_run_comm + dt
  end subroutine timer_comm_end


  !> Open the collective timer: the residual reduction and the control
  !> broadcasts. These are synchronisation points, so a rank that finished its
  !> own work early blocks here — which is where load imbalance shows up when
  !> the halo exchange does not reveal it.
  subroutine timer_sync_begin()
    t_sync_beg = timer_wtime()
  end subroutine timer_sync_begin


  !> Close the collective timer and accumulate.
  subroutine timer_sync_end()
    real(R8) :: dt
    dt = timer_wtime() - t_sync_beg
    t_sync_acc = t_sync_acc + dt
    t_run_sync = t_run_sync + dt
  end subroutine timer_sync_end


  !> Reduce the accumulators over the ranks. Collective: every rank must call
  !> it on the same iteration.
  function reduce_times(t_iter, t_comm, t_sync, n) result(s)
    use ARES_Mod_MPI, only: mpi_size_, mpi_allreduce_sum_r8, &
                            mpi_allreduce_min_r8, mpi_allreduce_max_r8, &
                            mpi_gather_r8
    real(R8), intent(in) :: t_iter  !< iteration time on this rank
    real(R8), intent(in) :: t_comm  !< exchange wait on this rank
    real(R8), intent(in) :: t_sync  !< collective wait on this rank
    integer,  intent(in) :: n       !< iterations the times were accumulated over
    type(stats_type)     :: s
    ! Local
    real(R8) :: tsum, csum, ssum, t_work
    real(R8), allocatable :: work(:)
    integer  :: r, ni

    ni = max(n, 1)

    ! Time this rank spent on its own cells: everything not blocked waiting for
    ! a neighbour's halo or for a collective.
    t_work = t_iter - t_comm - t_sync

    call mpi_allreduce_max_r8(t_iter, s%tmax)
    call mpi_allreduce_min_r8(t_iter, s%tmin)
    call mpi_allreduce_sum_r8(t_iter, tsum)
    call mpi_allreduce_sum_r8(t_comm, csum)
    call mpi_allreduce_sum_r8(t_sync, ssum)

    allocate(work(max(mpi_size_, 1)))
    work = 0.0_R8
    call mpi_gather_r8(t_work, work)

    s%tavg = tsum / real(mpi_size_, R8)
    if (s%tavg > 0.0_R8) s%imbalance = (s%tmax - s%tavg) / s%tavg * 100.0_R8
    if (tsum   > 0.0_R8) s%commfrac  = csum / tsum * 100.0_R8
    if (tsum   > 0.0_R8) s%syncfrac  = ssum / tsum * 100.0_R8

    s%wmin = work(1); s%wmax = work(1)
    do r = 1, mpi_size_
      if (work(r) > s%wmax) then; s%wmax = work(r); s%rmax = r - 1; end if
      if (work(r) < s%wmin) then; s%wmin = work(r); s%rmin = r - 1; end if
      s%wavg = s%wavg + work(r)
    end do
    s%wavg = s%wavg / real(mpi_size_, R8)
    if (s%wavg > 0.0_R8) s%wspread = (s%wmax - s%wavg) / s%wavg * 100.0_R8

    deallocate(work)

    ! Per-iteration figures
    s%tmax = s%tmax / real(ni, R8); s%tmin = s%tmin / real(ni, R8)
    s%tavg = s%tavg / real(ni, R8)
    s%wmax = s%wmax / real(ni, R8); s%wmin = s%wmin / real(ni, R8)
    s%wavg = s%wavg / real(ni, R8)

  end function reduce_times


  !> Report the current window and reset it. Collective.
  subroutine timer_report(level, iter)
    use ARES_Mod_MPI, only: mpi_is_root, mpi_size_
    integer, intent(in) :: level   !< multigrid level the iterations were run on
    integer, intent(in) :: iter    !< iteration counter on that level
    ! Local
    type(stats_type) :: s

    s = reduce_times(t_iter_acc, t_comm_acc, t_sync_acc, n_iter_acc)

    if (mpi_is_root) then
      write(*,'(A,I0,A,I0,A,I0,A,ES11.4,A,ES11.4,A,ES11.4,A,F7.1,A,F7.1,A,F7.1,A)') &
        ' ARES Timing | Level ', level, ' | Iter ', iter, ' | ', &
        max(n_iter_acc, 1), ' iters | wall/iter ', s%tmax, &
        ' s | rank min ', s%tmin, ' avg ', s%tavg, &
        ' | imbalance ', s%imbalance, &
        ' % | exchange wait ', s%commfrac, &
        ' % | collective wait ', s%syncfrac, ' %'

      if (mpi_size_ > 1) &
        write(*,'(A,ES11.4,A,I0,A,ES11.4,A,I0,A,ES11.4,A,F7.1,A)') &
          '  ARES Ranks | compute/iter max ', s%wmax, &
          ' s (rank ', s%rmax, ') | min ', s%wmin, &
          ' s (rank ', s%rmin, ') | mean ', s%wavg, &
          ' s | spread ', s%wspread, ' %'
    end if

    t_iter_acc = 0.0_R8
    t_comm_acc = 0.0_R8
    t_sync_acc = 0.0_R8
    n_iter_acc = 0

  end subroutine timer_report


  !> End-of-run timing summary. Collective, so every rank must call it.
  !> `Solver` is the time inside the iteration loop and is the figure to quote
  !> when timing the code; `Elapsed` adds what is spent between iterations,
  !> which is dominated by solution output.
  subroutine timer_summary()
    use ARES_Mod_MPI, only: mpi_is_root, mpi_size_, mpi_allreduce_max_r8
    ! Local
    type(stats_type) :: s
    real(R8) :: elapsed

    ! Both figures are taken on the critical path, i.e. the slowest rank
    call mpi_allreduce_max_r8(timer_wtime() - t_run_beg, elapsed)
    s = reduce_times(t_run_iter, t_run_comm, t_run_sync, n_run)

    if (mpi_is_root) then
      write(*,*)
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Timing'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A,T35,I0)')          '   Iterations', n_run
      write(*,'(A,T35,ES12.5,A)')    '   Solver', s%tmax * real(max(n_run,1), R8), ' s'
      write(*,'(A,T35,ES12.5,A)')    '   Solver per iteration', s%tmax, ' s'
      write(*,'(A,T35,ES12.5,A)')    '   Elapsed', elapsed, ' s'
      if (mpi_size_ > 1) then
        write(*,'(A,T35,F7.1,A)')    '   Load imbalance', s%imbalance, ' %'
        write(*,'(A,T35,F7.1,A)')    '   Exchange wait', s%commfrac, ' %'
        write(*,'(A,T35,F7.1,A)')    '   Collective wait', s%syncfrac, ' %'
        write(*,'(A,T35,F7.1,A)')    '   Compute spread over ranks', s%wspread, ' %'
      end if
      write(*,'(A)') ' ========================================================================================='
    end if

  end subroutine timer_summary

end module ARES_Mod_Timers
