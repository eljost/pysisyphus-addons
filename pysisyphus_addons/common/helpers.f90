module mod_pa_helpers

  use mod_pa_constants, only: dp, i4

  implicit none

  interface print_array
    module procedure print_array2d
    module procedure print_array3d
  end interface print_array

contains

  subroutine print_array2d (arr)
    ! Helper procedure to pretty-print an double array row-by-row for debugging.
    real(dp), intent(in) :: arr(:, :)
    integer(i4) :: i, nrows

    nrows = size(arr, 1)
    do i = 1, nrows
      print "(*(F8.4,:',',1x))", arr(i, :)
    end do
  end subroutine print_array2d

  subroutine print_array3d (arr)
    real(dp), intent(in) :: arr(:, :, :)
    integer(i4) :: i

    do i = 1, size(arr, 3)
      call print_array(arr(:, :, i))
    end do
  end subroutine print_array3d
end module mod_pa_helpers
