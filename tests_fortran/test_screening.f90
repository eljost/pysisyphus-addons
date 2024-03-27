 module test_screening
  use testdrive, only : new_unittest, unittest_type, error_type, check

  use mod_pa_constants, only: dp
  use mod_pa_screening, only: contracted_center

  implicit none

  private

  public :: collect_screening

contains

!> Collect all exported unit tests
subroutine collect_screening(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("contracted_center_prim", test_contracted_center_prim), &
    new_unittest("contracted_center_contr_sym", test_contracted_center_contr_sym), &
    new_unittest("contracted_center_contr_asym", test_contracted_center_contr_asym) &
  ]

end subroutine collect_screening

subroutine test_contracted_center_prim(error)
  type(error_type), allocatable, intent(out) :: error

  real(dp), parameter :: A(3) = (/ 0, 0, 0 /)
  real(dp), parameter :: axs(1) = (/ 1.0 /)
  real(dp) :: R(3)

  R = contracted_center(axs, axs, A, axs, axs, A)
  call check(error, abs(norm2(R-A)) <= 1e-14)
end subroutine test_contracted_center_prim

subroutine test_contracted_center_contr_sym(error)
  type(error_type), allocatable, intent(out) :: error

  real(dp), parameter :: A(3) = (/ 0, 0, 0 /)
  real(dp), parameter :: B(3) = (/ 1, 1, 1 /)
  real(dp), parameter :: axs(2) = (/ 1.0, 0.5 /)
  real(dp) :: R(3)

  R = contracted_center(axs, axs, A, axs, axs, B)
  call check(error, abs(norm2(R-(A + B) / 2)) <= 1e-14)
end subroutine test_contracted_center_contr_sym

subroutine test_contracted_center_contr_asym(error)
  type(error_type), allocatable, intent(out) :: error

  real(dp), parameter :: A(3) = (/ 0, 0, 0 /)
  real(dp), parameter :: B(3) = (/ 1, 1, 1 /)
  real(dp), parameter :: axs(3) = (/ 1.0, 0.5, 0.25 /)
  real(dp), parameter :: bxs(2) = (/ 1.0, 0.5 /)
  real(dp) :: R(3)
  real(dp) :: ri = 0.53650793650794d0

  ! Reuse exponents as contraction coefficients
  R = contracted_center(axs, axs, A, bxs, bxs, B)
  call check(error, abs(norm2(R-(/ ri, ri, ri /))) <= 1e-14)
end subroutine test_contracted_center_contr_asym

end module test_screening
