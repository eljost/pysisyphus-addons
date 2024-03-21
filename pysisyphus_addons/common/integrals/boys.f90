module mod_pa_boys
   use iso_fortran_env, only: int32, int64, real32

   use mod_pa_constants, only: dp, PI
   use mod_boys_data, only: boys_table

   implicit none

   private

   public :: boys_init, boys

   !real(dp), parameter :: PI = 4*atan(1.0d0)
   integer(kind=int32), parameter :: BOYS_N_MAX = 15
   ! Store double factorials as reals, as they are later divided by an integer
   ! and we don't want to do integer arithmetic.
   real(dp) :: double_facts(0:2*BOYS_N_MAX)
   !real(dp) :: boys_table(0:64, 2711)
   integer(kind=int32) :: i
   ! Constants, associated with the precaclulated Boys table.
   ! One row per N, one column per x-step.
   ! Step between two successive x_values
   real(dp), parameter :: table_step = 0.01d0
   integer(kind=int32), parameter :: table_factor = 1/table_step
   ! Beyond this x-value the boys_factorial is used.
   real(dp), parameter :: table_x_asym = 27.0d0
   integer(kind=int32), parameter :: table_add_points = 10
   ! Number of x-points in the table
   integer(kind=int32), parameter :: table_points = table_x_asym/table_step + table_add_points
   ! The x-values in the table
   real(dp), parameter :: xs_table(table_points + 1) = [(0.0 + i*table_step, i=0, table_points)]

contains

   ! https://en.wikipedia.org/wiki/Double_factorial
   real(kind=int64) function double_factorial(n)
      integer(kind=int32) :: N, i

      double_factorial = 1
      ! TODO: check single precision kind
      do i = 0, (ceiling(real(N, kind=real32)/2, kind=int32) - 1)
         double_factorial = double_factorial*(N - 2*i)
      end do
   end function double_factorial

   subroutine boys_init()
      integer(kind=int32) :: i
      ! Precompute double factorials. Maybe this is not needed with Fortran,
      ! but in the Python version this was a useful idea.
      do i = 0, 2*BOYS_N_MAX
         double_facts(i) = double_factorial(i)
      end do
   end subroutine boys_init

   ! Factorial boys for x >= 27.0
   pure real(dp) function boys_factorial(N, x)
      integer(kind=int32), intent(in) :: N
      real(dp), intent(in) :: x

      boys_factorial = double_facts(max(0, 2*N - 1))/(2**(N + 1))*sqrt(PI/(x**(2*N + 1)))
   end function boys_factorial

   ! 5-point Neville interpolation of pretabulated Boys-values
   ! Adapted from the code-generated version in pysisyphus.
   pure real(dp) function boys_neville(N, x)
      integer(kind=int32), intent(in) :: N
      real(dp), intent(in) :: x
      real(dp) :: xs(5), ys(5)
      integer(kind=int32) :: x_closest
      integer(kind=int32) :: indices(5)
      real(dp) :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11

      x_closest = 1 + int(x*table_factor, kind=int32)
      indices = x_closest + (/0, 1, 2, 3, 4/)
      xs = xs_table(indices)
      ys = boys_table(N, indices)
      x0 = -xs(5)
      x1 = x - xs(1)
      x2 = -xs(2)
      x3 = x + x2
      x4 = -xs(3)
      x5 = x + x4
      x6 = -xs(4)
      x7 = x + x6
      x8 = x + x0
      x9 = (x5*ys(4) - x7*ys(3))/(-x6 - xs(3))
      x10 = (x3*ys(3) - x5*ys(2))/(-x4 - xs(2))
      x11 = (-x10*x7 + x3*x9)/(-x6 - xs(2))

      !&< The following lines are not improved w/ fprettify ...
      boys_neville = (x1* &
          (-x11*x8 + x3*(x5*(x7*ys(5) - x8*ys(4))/(-x0 - xs(4)) - x8*x9) /(-x0 - xs(3))) &
          /(-x0 - xs(2)) - x8 &
          *(x1*x11 - x7 *(x1*x10 - x5*(x1*ys(2) - x3*ys(1))/(-x2 - xs(1)))  /(-x4 - xs(1)) ) &
          /(-x6 - xs(1)))/(-x0 - xs(1))
      !&>
   end function boys_neville

   ! Wrapper for Boys function that calls the appropriate implementation
   pure real(dp) function boys(N, x)
      integer(kind=int32), intent(in) :: N
      real(dp), intent(in) :: x

      if (x >= 27.0) then
         boys = boys_factorial(N, x)
      else
         boys = boys_neville(N, x)
      end if
   end function boys
end module mod_pa_boys
