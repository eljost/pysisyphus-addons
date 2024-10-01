module mod_pa_grid_wrapper
  use iso_c_binding, only: c_bool, c_int32_t, c_double

  use mod_pa_shells, only: t_shells
  use mod_pa_grid, only: eval_densities

  implicit none

  contains
    subroutine f_eval_densities(&
      nshells, ndata, bas_centers, bas_spec, bas_data, &
      npoints, grid3d, ndens, nbfs, densities, grid_densities, blk_size, thresh) &
      bind(c, name="f_eval_densities")

      ! AO basis
      integer(c_int32_t), value, intent(in) :: nshells, ndata
      integer(c_int32_t), intent(in) ::  bas_centers(nshells, 3)
      integer(c_int32_t), intent(in) ::  bas_spec(nshells, 5)
      real(c_double), intent(in) :: bas_data(ndata)
      ! Number of grid points and their coordinates
      integer(c_int32_t), value, intent(in) ::  npoints
      real(c_double), intent(in) :: grid3d(3, npoints)
      ! Number of densities and basis functions in the grid
      integer(c_int32_t), value, intent(in) ::  ndens
      integer(c_int32_t), value, intent(in) ::  nbfs
      ! Density matrix and density matrix on grid
      real(c_double), intent(in) :: densities(ndens, nbfs, nbfs)
      real(c_double), intent(in out) :: grid_densities(ndens, npoints)
      ! Optional arguments
      integer(c_int32_t), value, intent(in) :: blk_size
      real(c_double), value, intent(in) :: thresh
      ! Shells
      type(t_shells) :: shells

      ! Construct shells and evaluate densities
      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      call eval_densities(shells, grid3d, densities, grid_densities, blk_size, thresh)

    end subroutine f_eval_densities
end module mod_pa_grid_wrapper
