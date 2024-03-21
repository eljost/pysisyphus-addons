module mod_pa_diabatization_intor_wrapper

   use iso_c_binding, only: c_int32_t, c_double

   use mod_pa_shells, only: t_shells
   use mod_pa_diabatization, only: contract_coulomb_densities_2d, contract_coulomb_densities_4d

   implicit none

contains

   subroutine f_contract_coulomb_densities_2d( &
      nshells, ndata, bas_centers, bas_spec, bas_data, &
      nshells_aux, ndata_aux, bas_centers_aux, bas_spec_aux, bas_data_aux, &
      ndens, nbfs, densities, naux, df_tensor &
      ) bind(c, name="f_contract_coulomb_densities_2d")
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
      ! Densities
      integer(c_int32_t), value, intent(in) :: ndens
      integer(c_int32_t), value, intent(in) :: nbfs
      real(c_double), intent(in) :: densities(ndens, nbfs, nbfs)
      ! DF tensor
      integer(c_int32_t), value, intent(in) :: naux
      real(c_double), intent(in out) :: df_tensor(naux, ndens)

      type(t_shells) :: shells
      type(t_shells) :: shells_aux

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      shells_aux = t_shells(nshells_aux, bas_centers_aux, bas_spec_aux, bas_data_aux)

      call contract_coulomb_densities_2d(shells, shells_aux, densities, df_tensor)

   end subroutine f_contract_coulomb_densities_2d

   subroutine f_contract_coulomb_densities_4d( &
      nshells, ndata, bas_centers, bas_spec, bas_data, &
      nshells_aux, ndata_aux, bas_centers_aux, bas_spec_aux, bas_data_aux, &
      ndens, nbfs, densities, nstates, coulomb_tensor &
      ) bind(c, name="f_contract_coulomb_densities_4d")
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
      ! Densities
      integer(c_int32_t), value, intent(in) :: ndens
      integer(c_int32_t), value, intent(in) :: nbfs
      real(c_double), intent(in) :: densities(ndens, nbfs, nbfs)
      ! Coulomb tensor
      integer(c_int32_t), value, intent(in) :: nstates
      real(c_double), intent(in out) :: coulomb_tensor(nstates, nstates, nstates, nstates)

      type(t_shells) :: shells
      type(t_shells) :: shells_aux

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      shells_aux = t_shells(nshells_aux, bas_centers_aux, bas_spec_aux, bas_data_aux)

      call contract_coulomb_densities_4d(shells, shells_aux, densities, coulomb_tensor)

   end subroutine f_contract_coulomb_densities_4d
end module mod_pa_diabatization_intor_wrapper
