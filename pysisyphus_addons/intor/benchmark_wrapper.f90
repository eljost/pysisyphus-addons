module mod_pa_intor_benchmark_wrapper
   use iso_c_binding, only: c_int32_t, c_double

   use mod_pa_constants, only: dp
   use mod_pa_shells, only: t_shells
   use mod_pa_intor_benchmark, only: benchmark_int3c2e

   implicit none

contains

   subroutine f_benchmark_int3c2e( &
      nshells, ndata, bas_centers, bas_spec, bas_data, &
      nshells_aux, ndata_aux, bas_centers_aux, bas_spec_aux, bas_data_aux) &
      bind(c, name="f_benchmark_int3c2e")
      ! AO basis
      integer(c_int32_t), value, intent(in) :: nshells, ndata
      integer(c_int32_t), intent(in) ::  bas_centers(nshells, 3)
      integer(c_int32_t), intent(in) ::  bas_spec(nshells, 5)
      real(c_double), intent(in) :: bas_data(ndata)
      ! Auxiliary AO basis
      integer(c_int32_t), value, intent(in) :: nshells_aux, ndata_aux
      integer(c_int32_t), intent(in) ::  bas_centers_aux(nshells_aux, 3)
      integer(c_int32_t), intent(in) ::  bas_spec_aux(nshells_aux, 5)
      real(c_double), intent(in) :: bas_data_aux(ndata_aux)

      type(t_shells) :: shells
      type(t_shells) :: shells_aux

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      shells_aux = t_shells(nshells_aux, bas_centers_aux, bas_spec_aux, bas_data_aux)

      call benchmark_int3c2e(shells, shells_aux)

   end subroutine f_benchmark_int3c2e
end module mod_pa_intor_benchmark_wrapper
