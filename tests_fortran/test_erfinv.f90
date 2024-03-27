module test_erfinv
   use testdrive, only: new_unittest, unittest_type, error_type, check

   use mod_pa_erfinv, only: erfinv

   implicit none

   private

   public :: collect_erfinv

contains

   subroutine collect_erfinv(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("test_erfinv_eval", test_erfinv_eval) &
      ]

   end subroutine collect_erfinv

   subroutine test_erfinv_eval(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, erfinv(0d0) == 0d0)
      call check(error, erfinv(0.9d0) == 1.1630871536766727d0)
      call check(error, erfinv(-0.9d0) == -1.1630871536766727d0)
      call check(error, erfinv(0.5d0), 0.4769362762044698d0)
      call check(error, erfinv(-0.5d0), -0.4769362762044698d0)
   end subroutine test_erfinv_eval

end module test_erfinv
