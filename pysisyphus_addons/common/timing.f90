module mod_pa_timing

   use mod_pa_constants, only: i4, dp

   implicit none

   private

   public :: wtimer

contains

  !! Helper procedure for measuring wall times in seconds with millisecond resolution.
   subroutine wtimer(wtime, init)
      real(dp), intent(in out):: wtime
      !! Variable that will hold the elapsed wall time in seconds after the
      !! procedure call.
      logical, optional :: init
      !! Whether we wan't to start measuring the wall time in wtime. For the
      !! first call init should be .true. . For the second call it should be
      !! omitted or set to .false. .

      ! Intermediate variables
      real(dp) :: factor
      !! Multiplicative factor that will be -1 when init is .true.
      logical :: act_init
      integer(i4) :: count_, count_rate, count_max

      if (present(init)) then
         act_init = init
      else
         act_init = .false.
      end if

      ! Initialize wall time and determine factor
      if (act_init) then
         wtime = 0d0
         factor = -1d0
      else
         factor = 1d0
      end if

      call system_clock(count_, count_rate, count_max)
      wtime = wtime + factor*count_/real(count_rate, dp)
   end subroutine wtimer
end module mod_pa_timing
