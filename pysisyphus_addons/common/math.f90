module mod_pa_math

   use mod_pa_constants, only: i4, i8

   implicit none

contains

   ! https://en.wikipedia.org/wiki/Double_factorial
   integer(i8) function double_factorial(n)
      integer(i4) :: i, N

      double_factorial = 1
      do i = 0, (ceiling(real(N)/2, kind=i4) - 1)
         double_factorial = double_factorial*(N - 2*i)
      end do
   end function double_factorial

end module mod_pa_math
