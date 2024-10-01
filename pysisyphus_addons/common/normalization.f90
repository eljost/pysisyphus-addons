module mod_pa_normalization
  use iso_fortran_env, only: stderr=>error_unit

  use mod_pa_constants, only: dp, PI

  implicit none

contains

  subroutine norm_cgto_coeffs(coeffs, exps, L)
    real(dp), intent(in out) :: coeffs(:)
    real(dp), intent(in) :: exps(:)
    integer, intent(in) :: L
    real(dp) :: N
    integer :: i, j

    if (L < 0 ) then
      write(stderr, *) "L must be >= 0!"
      error stop
    end if

    if (size(exps) /= size(coeffs)) then
      write(stderr, *) "coeffs and exps must have the same size!"
      error stop
    end if

    N = 0.0
    do i = 1, size(exps)
        do j = 1, size(exps)
            N = N + coeffs(i) * coeffs(j) / (exps(i) + exps(j))**(L + 1.5) &
                * sqrt(exps(i) * exps(j))**(L + 1.5)
        end do
    end do

    do i = 1, size(exps)
        coeffs(i) = coeffs(i) * sqrt(exps(i)**(L + 1.5) / (PI**1.5 / 2**L * N))
    end do
  end subroutine norm_cgto_coeffs

end module mod_pa_normalization
