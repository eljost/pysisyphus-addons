! [1]   https://doi.org/10.1063/1.1567253
!       Fast evaluation of the Coulomb potential for electron densities
!       using multipole accelerated resolution of identity approximation
!       Sierka, Hogekamp, Ahlrichs, 2003
! [2]   https://doi.org/10.1063/1.3693908
!       Distance-dependent Schwarz-based integral estimates for
!       two-electron integrals: Reliable tightness vs. rigorous upper bounds
!       Maurer, Lambrecht, Flaig, Ochsenfeld, 2012
! [3]   https://doi.org/10.1063/1.4983393
!       Efficient evaluation of three-center Coulomb integrals
!       Samu, Kállay, 2017
! [4]   A tight distance-dependent estimator for screening three-center Coulomb integrals
!       over Gaussian basis functions
!       Hollman, Schaefer, Valeev, 2015

module mod_pa_screening
   use mod_pa_constants, only: dp, i4, PI
   use mod_pa_erfinv, only: erfinv
   use mod_pa_math, only: double_factorial
   use mod_pa_shells, only: t_shell, t_shells
   use mod_pa_intor, only: intor_schwarz_aux, intor_schwarz_bra

   implicit none

   type t_df_screener
      type(t_shells) :: shells
      type(t_shells) :: shells_aux
      real(dp) :: thresh_ws
      ! Number of shell pairs in the bra and number of auxiliary shells
      integer :: nbra
      integer :: naux
      ! Extents
      real(dp), allocatable :: extents_bra(:)
      real(dp), allocatable :: extents_aux(:)
      ! Centers
      real(dp), allocatable :: centers_bra(:, :)
      real(dp), allocatable :: centers_aux(:, :)
      ! Schwarz integral norms
      real(dp), allocatable :: schwarz_bra(:)
      real(dp), allocatable :: schwarz_aux(:)
      ! Effective exponents
      real(dp), allocatable :: effexps_bra(:)
      real(dp), allocatable :: effexps_aux(:)
      ! Common factor beta, eq. (28) in [4]
      real(dp), allocatable :: betas_aux(:)
      ! Factor in the denominator of the QVL-estimator (exp_a + exp_b)**0.25
      real(dp), allocatable :: effexps_bra_denom(:)

   contains
      procedure :: init
      procedure :: schwarz_estimate
      procedure :: qvl_estimate
      procedure :: debug_estimate
   end type t_df_screener

   interface t_df_screener
      module procedure :: t_df_screener_constructor
   end interface t_df_screener

contains

   type(t_df_screener) function t_df_screener_constructor(shells, shells_aux, thresh_ws) &
      result(screener)
      type(t_shells), intent(in) :: shells
      type(t_shells), intent(in) :: shells_aux
      real(dp), intent(in) :: thresh_ws

      screener%shells = shells
      screener%shells_aux = shells_aux
      screener%thresh_ws = thresh_ws
      ! sum_1^n = n * (n + 1) / 2
      screener%nbra = (shells%nshells*(shells%nshells + 1))/2
      screener%naux = shells_aux%nshells

      ! Extents
      allocate (screener%extents_bra(screener%nbra))
      allocate (screener%extents_aux(screener%naux))

      ! Centers
      allocate (screener%centers_bra(3, screener%nbra))
      allocate (screener%centers_aux(3, screener%naux))

      ! Schwarz integrals
      !
      ! They will contain the square roots of the norms of the integrals
      ! over shell pairs ((C|C), schwarz_aux) or shell quadruplets ((ab|ab), schwarz_bra).
      allocate (screener%schwarz_bra(screener%nbra))
      allocate (screener%schwarz_aux(screener%naux))

      ! Effective exponents; just the minimum exponent of each shell
      allocate (screener%effexps_bra(screener%nbra))
      allocate (screener%effexps_aux(screener%naux))

      ! Common factor beta, depending on L and an exponent; eq. (28) in [4]
      allocate (screener%betas_aux(screener%naux))

      ! (effexp_a + effexp_b)**0.25
      allocate (screener%effexps_bra_denom(screener%nbra))
   end function t_df_screener_constructor

   subroutine init(self)
      class(t_df_screener), intent(in out) :: self
      type(t_shell) :: shell_a, shell_b
      integer :: i, j, k

      ! Loop over all auxiliary shells and set their centers and calculate their extents
      do i = 1, self%naux
         shell_a = self%shells_aux%shells(i)
         self%centers_aux(:, i) = shell_a%center
         self%extents_aux(i) = auxiliary_extent(self%shells_aux%get_exps(i), 1d-8)
         ! Effective exponent (smallest exponent)
         self%effexps_aux(i) = minval(self%shells_aux%get_exps(i))
         self%betas_aux(i) = self%effexps_aux(i)**(-(2d0*shell_a%L + 3d0)/4d0) &
                             *sqrt(real(double_factorial(2*shell_a%L - 1), dp))
      end do

      ! Loop over all unique principal shell pairs and calculate their
      ! (contracted) centers & extents.
      ! Store them in lower triangular form:
      !   (1, 1)
      !   (2, 1), (2, 2)
      !   (3, 1), (3, 2), (3, 3)
      !   ...
      ! Stored this way, 2d indices (i, j) are mapped to a 1d index via (i * (i-1) / 2 + j)
      do i = 1, self%shells%nshells
         shell_a = self%shells%shells(i)
         ! Effective exponent (smallest exponent)
         self%effexps_bra(i) = minval(self%shells%get_exps(i))
         do j = 1, i
            ! bra index in 1d list
            k = i*(i - 1)/2 + j
            shell_b = self%shells%shells(j)
            self%centers_bra(:, k) = contracted_center( &
                                     self%shells%get_exps(i), self%shells%get_coeffs(i), shell_a%center, &
                                     self%shells%get_exps(j), self%shells%get_coeffs(j), shell_b%center &
                                     )
            self%extents_bra(k) = contracted_extent( &
                                  self%shells%get_exps(i), self%shells%get_coeffs(i), shell_a%center, &
                                  self%shells%get_exps(j), self%shells%get_coeffs(j), shell_b%center, &
                                  self%thresh_ws &
                                  )
            self%effexps_bra_denom(k) = (self%effexps_bra(i) + self%effexps_bra(j))**0.25d0
         end do
      end do

      ! Calculate Schwarz integrals; see eq. (B1) in [3]
      !
      !                           bra            ket
      ! ||(Li Lj | Lk)||²_2 <= ||(Li Lj | Li Lj)||_2 ||(Lk | Lk)||_2
      !
      ! ||(Li Lj | Lk)||_2 <= sqrt(||(Li Lj | Li Lj)||_2 ||(Lk | Lk)||_2)
      !
      ! This requires 4-center-2-electron integrals for the bra and
      ! 2-center-2-electron integrals for the ket/auxiliary shell.
      !
      ! Here, the square roots of the norms will be kept

      ! 4-center-2-electron integrals
      self%schwarz_bra = intor_schwarz_bra(self%shells)

      ! 2-center-2-electron integrals
      self%schwarz_aux = intor_schwarz_aux(self%shells_aux)
   end subroutine init

   ! Extent of a (contracted) auxiliary function. See Eq. (34) in [1].
   ! Typcial values for prec will be between 1e-9 (more precise) and
   ! 1e-6 (less precise). See section C. in [1] for a discussion.
   real(dp) function auxiliary_extent(exponents, prec, min_extent)
      real(dp), intent(in) :: exponents(:)
      real(dp), intent(in) :: prec
      real(dp), optional, intent(in) :: min_extent
      real(dp) :: act_min_extent
      real(dp) :: prec_term
      real(dp) :: sqrt_arg
      real(dp) :: tmp_ext
      integer :: i

      if (present(min_extent)) then
         act_min_extent = min_extent
      else
         act_min_extent = 0.01d0
      end if

      auxiliary_extent = 0
      prec_term = -log(prec)

      do i = 1, size(exponents)
         sqrt_arg = (prec_term + 0.5d0*log(exponents(i)))/exponents(i)
         if (sqrt_arg < 0) then
            tmp_ext = act_min_extent
         else
            tmp_ext = sqrt(sqrt_arg)
         end if
         auxiliary_extent = max(auxiliary_extent, tmp_ext)
      end do
   end function auxiliary_extent

   real(dp) function contracted_extent(axs, das, A, bxs, dbs, B, thresh_ws)
      real(dp), intent(in) :: axs(:), bxs(:)
      real(dp), intent(in) :: das(:), dbs(:)
      real(dp), intent(in) :: A(3), B(3)
      real(dp), intent(in) :: thresh_ws
      real(dp) :: erfi
      real(dp) :: P(3), R(3)
      real(dp) :: ax, bx, px
      real(dp) :: tmp_ext
      real(dp) :: contr_prim_dist
      integer :: i, j

      ! Depending on the paper either erf⁻¹(1.0 - thresh) or erfc⁻¹(thresh) is used.
      ! Both approaches yield the same number.
      erfi = erfinv(1.0 - thresh_ws)
      R = contracted_center(axs, das, A, bxs, dbs, B)

      contracted_extent = 0
      do i = 1, size(axs)
         ax = axs(i)
         do j = 1, size(bxs)
            bx = bxs(j)
            ! Below, some amount of recalculation is involved, as the primitive centers
            ! are already calculated in 'contracted_center()'.
            ! Eq. (B1) in [2]
            px = ax + bx
            P = (ax*A + bx*B)/px
            ! Eq. (B2) in [2]
            tmp_ext = sqrt(2d0/px)*erfi
            contr_prim_dist = norm2(R - P)
            ! Eq. (B4) in [2]
            contracted_extent = max(contracted_extent, tmp_ext + contr_prim_dist)
         end do
      end do
   end function contracted_extent

   ! Center of a contracted charge distribution. See eq. (B2) in [2].
   function contracted_center(axs, das, A, bxs, dbs, B) result(R)
      real(dp), intent(in) :: axs(:), bxs(:)
      real(dp), intent(in) :: das(:), dbs(:)
      real(dp), intent(in) :: A(3), B(3)
      ! Exponents
      real(dp) :: ax, bx
      ! Contraction coefficients and realted quantities
      real(dp) :: da, db, abs_dab, denom
      ! Quantities related to the Gaussian product theorem
      real(dp) :: px, P(3)
      integer :: i, j
      real(dp) :: R(3)

      R = 0
      denom = 0
      do i = 1, size(axs)
         ax = axs(i)
         da = das(i)
         do j = 1, size(bxs)
            bx = bxs(j)
            db = dbs(j)
            px = ax + bx
            ! Center of primitive charge distribution
            P = (ax*A + bx*B)/px
            abs_dab = abs(da*db)
            R = R + abs_dab*P
            denom = denom + abs_dab
         end do
      end do
      R = R/denom
   end function contracted_center

   real(dp) function qvl_estimate(self, i, j, k)
      class(t_df_screener), intent(in) :: self
      integer(i4), intent(in) :: i, j, k

      real(dp) :: R
      real(dp) :: ext_ij, ext_k, ext
      real(dp) :: schwarz_bra, schwarz_aux
      integer(i4) :: bra_index

      ! Map 2d index (i, j) onto 1d list storing the lower triangular matrix
      ! of a symmetric matrix in row-major order.
      bra_index = i*(i - 1)/2 + j
      ! Distance between bra charge distribution and center of auxiliary shell
      R = norm2(self%centers_bra(:, bra_index) - self%centers_aux(:, k))
      ! Extents
      ext_ij = self%extents_bra(bra_index)
      ext_k = self%extents_aux(k)
      ext = ext_ij + ext_k

      schwarz_bra = self%schwarz_bra(bra_index)
      schwarz_aux = self%schwarz_aux(k)

      ! Plain Schwarz estimate
      if (R <= ext) then
         qvl_estimate = schwarz_bra*schwarz_aux
      ! QVL estimate
      else
         qvl_estimate = schwarz_bra*&
                        PI*sqrt(2d0)*self%betas_aux(k) &
                        /(self%effexps_bra_denom(bra_index)*R**(self%shells_aux%shells(k)%L + 1))
      end if

   end function qvl_estimate

   real(dp) function schwarz_estimate(self, i, j, k)
      class(t_df_screener), intent(in) :: self
      integer(i4), intent(in) :: i, j, k
      real(dp) :: schwarz_bra, schwarz_aux
      integer(i4) :: bra_index

      ! Map 2d index (i, j) onto 1d list storing the lower triangular matrix
      ! of a symmetric matrix in row-major order.
      bra_index = i*(i - 1)/2 + j
      schwarz_bra = self%schwarz_bra(bra_index)
      schwarz_aux = self%schwarz_aux(k)
      schwarz_estimate = schwarz_bra*schwarz_aux
   end function schwarz_estimate

   subroutine debug_estimate(self, i, j, k, schwarz_est, qvl_est, R)
      class(t_df_screener), intent(in) :: self
      integer(i4), intent(in) :: i, j, k
      real(dp), intent(out) :: schwarz_est, qvl_est, R

      integer(i4) :: bra_index

      schwarz_est = self%schwarz_estimate(i, j, k)
      qvl_est = self%qvl_estimate(i, j, k)
      bra_index = i*(i - 1)/2 + j
      R = norm2(self%centers_bra(:, bra_index) - self%centers_aux(:, k))
   end subroutine debug_estimate
end module mod_pa_screening
