module mod_pa_linalg
  use mod_pa_constants, only: dp, i4

  implicit none

contains
  subroutine eigh(mat, eigvals, eigvecs)
    real(dp), intent(in) :: mat(:, :)
    real(dp), allocatable, intent(out) :: eigvals(:)
    real(dp), allocatable, intent(out) :: eigvecs(:, :)
    real(dp), allocatable :: work(:)
    integer(i4), allocatable :: iwork(:)
    integer(i4) :: n, lwork, liwork, info

    n = size(mat, 1)

    allocate(eigvals(n))
    allocate(eigvecs(n, n))
    allocate(work(1))
    allocate(iwork(1))

    ! Workspace query to determine optimal lwork, liwork values
    lwork = -1
    liwork = -1
    call dsyevd("V", "U", n, eigvecs, n, eigvals, work, lwork, iwork, liwork, info)
    lwork = int(work(1), i4)
    liwork = iwork(1)
    deallocate(work, iwork)
    allocate(work(lwork), iwork(liwork))

    eigvecs = mat
    call dsyevd("V", "U", n, eigvecs, n, eigvals, work, lwork, iwork, liwork, info)
    
  end subroutine eigh

  subroutine matrix_powerh(mat, power)
    real(dp), intent(in out) :: mat(:, :)
    real(dp), intent(in) :: power
    real(dp), allocatable :: eigvals(:), eigvecs(:, :)
    real(dp), allocatable :: work(:, :)
    integer :: n, i

    allocate(work, mold=mat)
    n = size(mat, 1)

    ! Calculate eigenvalues and eigenvectors
    call eigh(mat, eigvals, eigvecs)
    eigvals = eigvals**power
    
    ! Turn 'mat' into diagonal matrix of eigenvalues
    mat = 0d0
    do i = 1, n
      mat(i, i) = eigvals(i)
    end do
    ! eigvecs @ diag(eigvals) @ eigvecs^T
    call dgemm("N", "N", n, n, n, 1d0, eigvecs, n, mat, n, 0d0, work, n)
    call dgemm("N", "T", n, n, n, 1d0, work, n, eigvecs, n, 0d0, mat, n)
  end subroutine matrix_powerh

end module mod_pa_linalg
