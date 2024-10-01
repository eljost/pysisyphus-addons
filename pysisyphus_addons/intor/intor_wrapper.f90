module mod_pa_intor_wrapper

   use iso_c_binding, only: c_int32_t, c_double

   use mod_pa_shells, only: t_shells
   use mod_pa_intor, only: intor_schwarz_aux, intor_schwarz_bra

   implicit none

contains

   subroutine f_int_schwarz_aux(nshells, ndata, bas_centers, bas_spec, bas_data, &
                                integrals) bind(c, name="f_int_schwarz_aux")
      integer(c_int32_t), value, intent(in) :: nshells, ndata
      integer(c_int32_t), intent(in) ::  bas_centers(nshells, 3)
      integer(c_int32_t), intent(in) ::  bas_spec(nshells, 5)
      real(c_double), intent(in) :: bas_data(ndata)
      real(c_double), intent(in out) :: integrals(nshells)
      type(t_shells) :: shells

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      integrals = intor_schwarz_aux(shells)

   end subroutine f_int_schwarz_aux

   subroutine f_int_schwarz_bra(nshells, ndata, bas_centers, bas_spec, bas_data, &
                                integrals) bind(c, name="f_int_schwarz_bra")
      integer(c_int32_t), value, intent(in) :: nshells, ndata
      integer(c_int32_t), intent(in) ::  bas_centers(nshells, 3)
      integer(c_int32_t), intent(in) ::  bas_spec(nshells, 5)
      real(c_double), intent(in) :: bas_data(ndata)
      real(c_double), intent(in out) :: integrals(nshells * (nshells + 1) / 2)
      type(t_shells) :: shells

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      integrals = intor_schwarz_bra(shells)

   end subroutine f_int_schwarz_bra

end module mod_pa_intor_wrapper
