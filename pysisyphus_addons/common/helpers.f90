module mod_pa_helpers

  use mod_pa_constants, only: dp, i4

  implicit none

contains

  subroutine print_array (arr)
    ! Helper procedure to pretty-print an double array row-by-row for debugging.
    real(dp), intent(in) :: arr(:, :)
    integer(i4) :: i, nrows

    nrows = size(arr, 1)
    do i = 1, nrows
      print "(*(F8.4,:',',1x))", arr(i, :)
    end do

  end subroutine print_array
end module mod_pa_helpers
