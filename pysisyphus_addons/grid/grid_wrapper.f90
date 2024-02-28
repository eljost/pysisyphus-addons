module mod_pa_grid_wrapper
  use iso_c_binding, only: c_bool, c_int32_t, c_double

  use mod_pa_shells, only: t_shells
  use mod_pa_grid, only: eval_density

  implicit none

  contains
    subroutine f_eval_density(&
      nshells, ndata, bas_centers, bas_spec, bas_data, &
      npoints, grid3d, nbfs, density, grid_dens, blk_size, thresh, accumulate) &
      bind(c, name="f_eval_density")

      ! AO basis
      integer(c_int32_t), value, intent(in) :: nshells, ndata
      integer(c_int32_t), intent(in) ::  bas_centers(nshells, 3)
      integer(c_int32_t), intent(in) ::  bas_spec(nshells, 5)
      real(c_double), intent(in) :: bas_data(ndata)
      ! Number of grid points and their coordinates
      integer(c_int32_t), value, intent(in) ::  npoints
      real(c_double), intent(in) :: grid3d(3, npoints)
      integer(c_int32_t), value, intent(in) ::  nbfs
      ! Density matrix and density matrix on grid
      real(c_double), intent(in) :: density(nbfs, nbfs)
      real(c_double), intent(in out) :: grid_dens(npoints)
      ! Optional arguments
      integer(c_int32_t), value, intent(in) :: blk_size
      real(c_double), value, intent(in) :: thresh
      logical(c_bool), value, intent(in) :: accumulate
      logical :: f_accumulate
      ! Shells
      type(t_shells) :: shells

      ! Convert 1-byte c boolean to 4-byte Fortran boolean
      f_accumulate = logical(accumulate)
      ! Construct shells and evaluate densities
      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      call eval_density(shells, grid3d, density, grid_dens, blk_size, thresh, f_accumulate)

    end subroutine f_eval_density
end module mod_pa_grid_wrapper
