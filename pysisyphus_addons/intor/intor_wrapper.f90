module mod_pa_intor_wrapper

   use iso_c_binding, only: c_int32_t, c_double

   use mod_pa_shells, only: t_shells
   use mod_pa_intor, only: intor_schwarz

   implicit none

contains

   subroutine f_int_schwarz(nshells, ndata, bas_centers, bas_spec, bas_data, &
                                    nbfs, integrals) bind(c, name="f_int_schwarz")
      integer(c_int32_t), value, intent(in) :: nshells, ndata
      integer(c_int32_t), intent(in) ::  bas_centers(nshells, 3)
      integer(c_int32_t), intent(in) ::  bas_spec(nshells, 5)
      real(c_double), intent(in) :: bas_data(ndata)
      integer(c_int32_t), value, intent(in) :: nbfs
      real(c_double), intent(in out) :: integrals(nbfs, nbfs)
      type(t_shells) :: shells

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      call shells%print_shells()

      integrals = intor_schwarz(shells)

   end subroutine f_int_schwarz

end module mod_pa_intor_wrapper
