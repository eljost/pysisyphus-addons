! [1]   https://doi.org/10.1063/1.1567253
!       Fast evaluation of the Coulomb potential for electron densities
!       using multipole accelerated resolution of identity approximation
!       Sierka, Hogekamp, Ahlrichs, 2003
! [2]   https://doi.org/10.1063/1.3693908
!       Distance-dependent Schwarz-based integral estimates for
!       two-electron integrals: Reliable tightness vs. rigorous upper bounds

module mod_pa_screening
  use mod_pa_constants, only: dp, i4
  use mod_pa_erfinv, only: erfinv
  use mod_pa_shells, only: t_shells

  implicit none

  type t_df_screener
    type(t_shells) :: shells
    type(t_shells) :: shells_aux
    real(dp) :: thresh_ws

    contains
      procedure :: init
      procedure :: triple_is_significant
  end type t_df_screener

  contains

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
        sqrt_arg = (prec_term + 0.5d0 * log(exponents(i))) / exponents(i)
        if (sqrt_arg < 0) then
          tmp_ext = act_min_extent
        else
          tmp_ext = sqrt(sqrt_arg)
        endif
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
          P = (ax * A + bx * B) / px
          ! Eq. (B2) in [2]
          tmp_ext = sqrt(2d0 / px) * erfi
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
          P = (ax * A + bx * B) / px
          abs_dab = abs(da * db)
          R = R + abs_dab * P
          denom = denom + abs_dab
        end do
      end do
      R = R / denom
    end function contracted_center

    subroutine init(self)
      class(t_df_screener), intent(in) :: self
        ! Calculate contracted/primitive centers for shells & aux_shells
        ! Calculate extents of principal basis shell pairs
        ! Calculate extents of auxiliary functions
        ! Calculate Schwarz integrals
    end subroutine init

    logical function triple_is_significant(self, i, j, k)
      class(t_df_screener), intent(in) :: self
      integer(i4), intent(in) :: i, j, k
      triple_is_significant = .true.
    end function triple_is_significant
end module mod_pa_screening
