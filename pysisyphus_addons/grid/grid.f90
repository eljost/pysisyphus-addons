module mod_pa_grid
   use ieee_arithmetic
   use mod_pa_constants, only: i4, dp
   use mod_pa_shells, only: t_shell, t_shells

   implicit none

    real(dp), parameter :: CUTOFF = -36.0d0
    real(dp), parameter :: NORM2 = 0.5773502691896258d0
    real(dp), parameter :: NORM3 = 0.2581988897471611d0
    real(dp), parameter :: NORM4 = 0.09759000729485332d0
    real(dp), parameter :: NORM22 = 0.3333333333333333d0

   contains

      subroutine cart_gto3d_rel(La, axs, das, RA, RA2, result)
        integer(i4), intent(in) :: La
        real(dp), intent(in) :: axs(:), das(:)
        real(dp), intent(in) :: RA(3), RA2(3)
        real(dp), intent(in out) :: result(:)
        real(dp) :: dx, dy, dz, dx2, dy2, dz2, dist, exparg, expterm
        real(dp) :: dx3, dy3, dz3, dx4, dy4, dz4
        integer(i4) :: i

        ! Initialize result array
        result = 0

        dx = RA(1)
        dy = RA(2)
        dz = RA(3)
        dx2 = RA2(1)
        dy2 = RA2(2)
        dz2 = RA2(3)
        dist = dx2 + dy2 + dz2

        if (La > 2) then
            dx3 = dx2 * dx
            dy3 = dy2 * dy
            dz3 = dz2 * dz
        end if
        if (La > 3) then
            dx4 = dx2 * dx2
            dy4 = dy2 * dy2
            dz4 = dz2 * dz2
        end if

        ! Loop over primitives
        do i = 1, size(axs)
            exparg = -axs(i) * dist
            if (exparg < CUTOFF) cycle
            ! Multipliy w/ contraction coefficient
            expterm = das(i) * exp(exparg)
            select case (La)
                case (0)
                    result(1) = result(1) + expterm
                case (1)
                    result(1) = result(1) + dx * expterm
                    result(2) = result(2) + dy * expterm
                    result(3) = result(3) + dz * expterm
                case (2)
                    result(1) = result(1) + NORM2 * dx2 * expterm
                    result(2) = result(2) + dx * dy * expterm
                    result(3) = result(3) + dx * dz * expterm
                    result(4) = result(4) + NORM2 * dy2 * expterm
                    result(5) = result(5) + dy * dz * expterm
                    result(6) = result(6) + NORM2 * dz2 * expterm
                case (3)
                    result(1) = result(1) + NORM3 * dx3 * expterm
                    result(2) = result(2) + NORM2 * dx2 * dy * expterm
                    result(3) = result(3) + NORM2 * dx2 * dz * expterm
                    result(4) = result(4) + NORM2 * dx * dy2 * expterm
                    result(5) = result(5) + dx * dy * dz * expterm
                    result(6) = result(6) + NORM2 * dx * dz2 * expterm
                    result(7) = result(7) + NORM3 * dy3 * expterm
                    result(8) = result(8) + NORM2 * dy2 * dz * expterm
                    result(9) = result(9) + NORM2 * dy * dz2 * expterm
                    result(10) = result(10) + NORM3 * dz3 * expterm
                case (4)
                    result(1) = result(1) + NORM4 * dx4 * expterm
                    result(2) = result(2) + NORM3 * dx3 * dy * expterm
                    result(3) = result(3) + NORM3 * dx3 * dz * expterm
                    result(4) = result(4) + NORM22 * dx2 * dy2 * expterm
                    result(5) = result(5) + NORM2 * dx2 * dy * dz * expterm
                    result(6) = result(6) + NORM22 * dx2 * dz2 * expterm
                    result(7) = result(7) + NORM3 * dx * dy3 * expterm
                    result(8) = result(8) + NORM2 * dx * dy2 * dz * expterm
                    result(9) = result(9) + NORM2 * dx * dy * dz2 * expterm
                    result(10) = result(10) + NORM3 * dx * dz3 * expterm
                    result(11) = result(11) + NORM4 * dy4 * expterm
                    result(12) = result(12) + NORM3 * dy3 * dz * expterm
                    result(13) = result(13) + NORM22 * dy2 * dz2 * expterm
                    result(14) = result(14) + NORM3 * dy * dz3 * expterm
                    result(15) = result(15) + NORM4 * dz4 * expterm
                case default
                    ! TODO: implement general formula for higher angular momenta
                    result = ieee_value(result, ieee_signaling_nan)
            end select
        end do
        ! End loop over primitives
    end subroutine cart_gto3d_rel

   subroutine eval_density(shells, grid3d, density, grid_dens, blk_size, thresh, accumulate)
      type(t_shells), intent(in) :: shells
      ! Grid points for density evaluation with shape (3, npoints)
      real(dp), intent(in) :: grid3d(:, :)
      ! Cartesian AO density matrix in pysisyphus-order, shape (ncartbfs, ncartbfs)
      real(dp), intent(in) :: density(:, :)
      ! Arry holding the density on the grid points with shape (npoints, )
      real(dp), intent(in out) :: grid_dens(:)

      ! Optional arguments
      !
      ! Block size; decrease for increased accuracy
      integer(i4), optional, intent(in) :: blk_size
      ! Threshold used for estimation whether a basis function contributes to a given
      ! grid point. Increase for increased accuracy.
      real(dp), optional, intent(in) :: thresh
      ! Boolean flag that controls wheter grid_dens is zeroed or not. When set to
      ! false, its default value, then grid_dens is zerod before density evaluation.
      logical, optional, intent(in) :: accumulate
      ! Actual values
      integer(i4) :: act_blk_size
      real(dp) :: act_thresh
      logical :: act_accumulate

      ! Number of points and number of blocks
      integer(i4) :: npoints, nblks
      ! Indices and sizes related to a given block
      integer(i4) :: cur_blk, cur_blk_start, cur_blk_end, cur_blk_size
      ! Various indices:
      !  i: Indexes grid points
      !  a: Indexes shells
      !  nu, mu: Index basis functions
      integer(i4) :: i, a, nu, mu
      ! Array holding Cartesian basis function values for a block
      real(dp), allocatable :: chis_cart(:, :), chis_mean(:), scratch(:, :)
      type(t_shell) :: shell
      ! Indexing pointing on the current basis function center A, used to calculate RA.
      integer(i4) :: cur_center_ind
      ! R: grid point
      ! RA: distance between R and basis function center A
      ! RA2: square of RA
      real(dp) :: R(3), RA(3), RA2(3)
      ! Mean value of a basis function on grid points in a block. Used to estimate
      ! the contribution of basis function pairs.
      real(dp) :: mean
      ! Factor that takes symmetry of density matrix into account
      real(dp) :: factor
      ! Logical array indicating which basis function contributed to a given grid point
      logical, allocatable :: nu_contributed(:)
      ! Estimated contribution of a basis function pair at a given grid point
      real(dp) :: contrib
      ! Element of the density matrix
      real(dp) :: dens_nm

      ! Handle optional arguments and set their default values
      if (present(blk_size)) then
         act_blk_size = blk_size
      else
         act_blk_size = 100
      end if
      if (present(thresh)) then
         act_thresh = thresh
      else
         act_thresh = 1e-8
      end if
      if (present(accumulate)) then
         act_accumulate = accumulate
      else
         act_accumulate = .false.
      end if

      ! Zero density on grid when we don't accumulate
      if (.not. accumulate) then
         grid_dens = 0
      end if

      allocate(chis_cart(blk_size, shells%ncartbfs))
      allocate(chis_mean(shells%ncartbfs))
      allocate(scratch(blk_size, shells%ncartbfs))
      allocate(nu_contributed(shells%ncartbfs))

      ! Determine number of points and number of blocks
      npoints = size(grid3d, 2)
      nblks = ceiling(real(npoints, dp) / real(blk_size, dp))

      ! Loop over all blocks
      do cur_blk = 1, nblks
         ! Determine start index, end index and size of a block. The last
         ! block may have a different size, because there may not be enough points left.
         cur_blk_start = (cur_blk - 1) * blk_size + 1
         cur_blk_end = min(cur_blk_start + blk_size - 1, npoints)
         cur_blk_size = cur_blk_end - cur_blk_start + 1

         ! Loop over all grid points in block
         do i = 1, cur_blk_size
            R = grid3d(:, cur_blk_start + i - 1)
            ! Loop over all shells
            do a = 1, shells%nshells
               shell = shells%shells(a)
               ! Check if the center of the shell moved compared to the previous shell.
               ! If so, update RA and RA2. 
               if ((shell%center_ind /= cur_center_ind) .or. (a == 1)) then
                  RA = R - shell%center
                  RA2 = RA**2
                  cur_center_ind = shell%center_ind
               end if
               ! Evaluate basis function shell at the grid point R
               call cart_gto3d_rel(shell%L, shells%get_exps(a), shells%get_coeffs(a), RA, RA2, &
                                   chis_cart(i, shell%cart_index:shell%cart_index_end))
            end do
            ! End loop over shells
         end do
         ! End loop over grid points in block

         ! TODO: Convert to spherical basis functions?!

         ! Calculate mean values of basis functions in the current block
         do nu = 1, shells%ncartbfs
            mean = 0
            do i = 1, cur_blk_size
               mean = mean + abs(chis_cart(i, nu))
            end do
            chis_mean(nu) = mean / cur_blk_size
         end do

         scratch = 0
         nu_contributed = .false.
         ! Loop over all basis function pairs
         do nu = 1, shells%ncartbfs
            do mu = nu, shells%ncartbfs
               if (nu == mu) then
                  factor = 1
               else
                  factor = 2
               end if
               ! Density matrix is symmetric
               dens_nm = density(mu, nu)
               ! Estimate contribution of basis function to current block
               contrib = factor * dens_nm * chis_mean(nu) * chis_mean(mu)
               if (abs(contrib) >= thresh) then
                  scratch(:, nu) = scratch(:, nu) + factor * dens_nm * chis_cart(:, mu)
                  nu_contributed(nu) = .true.
               end if
            end do
         end do
         ! End loop over all basis function pairs

         do nu = 1, shells%ncartbfs
            if (.not. nu_contributed(nu)) then
               cycle
            end if
            do i=1, cur_blk_size
            !grid_dens(cur_blk_start:cur_blk_end) = &
               !grid_dens(cur_blk_start:cur_blk_end) + scratch(i, nu) * chis_cart(i, nu)
            grid_dens(cur_blk_start + i - 1) = &
               grid_dens(cur_blk_start + i - 1) + scratch(i, nu) * chis_cart(i, nu)
            end do
         end do
      end do
      ! End loop over all blocks
   end subroutine eval_density

end module mod_pa_grid
