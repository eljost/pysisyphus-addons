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

  subroutine inv_chol(mat)
    ! Cholesky decomposition of symmetric positive-definite matrix mat into lower
    ! triangular matrix L, followed by inversion of L. The input matrix mat is
    ! overwritten in the process and will finally contain the inverse of L, L⁻¹.
    !
    ! This capability is useful for calculating density fitting metric. L⁻¹^T can be
    ! absorbed into the DF-tensor (rs|P).

    ! Symmetric, positive-definite matrix
    real(dp), intent(in out) :: mat(:, :)
    ! Number of rows/columns in mat
    integer(i4) :: N
    integer(i4) :: info

    N = size(mat, 1)

    ! Cholesky decomposition of mat into L L^T
    call dpotrf("L", N, mat, N, info)

    ! Invert lower triangular matrix
    call dtrtri("L", "N", N, mat, N, info)
  end subroutine inv_chol

end module mod_pa_linalg
