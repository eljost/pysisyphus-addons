module mod_pa_std_intor_wrapper

   use iso_c_binding, only: c_int32_t, c_double

   use mod_pa_shells, only: t_shells
   use mod_pa_std_intor, only: intor_eri2c, intor_df1c

   implicit none

contains

   subroutine f_intor_eri2c( &
      nshells, ndata, bas_centers, bas_spec, bas_data, &
      nshells_aux, ndata_aux, bas_centers_aux, bas_spec_aux, bas_data_aux, &
      nbfs, eri2c_tensor) bind(c, name="f_intor_eri2c")
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
      ! ERI-2center-tensor
      integer(c_int32_t), value, intent(in) :: nbfs
      real(c_double), intent(in out) :: eri2c_tensor(nbfs, nbfs)
      ! Shells
      type(t_shells) :: shells, shells_aux

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      shells_aux = t_shells(nshells_aux, bas_centers_aux, bas_spec_aux, bas_data_aux)
      call intor_eri2c(shells, shells_aux, eri2c_tensor)
   end subroutine f_intor_eri2c

   subroutine f_intor_df1c( &
      nshells, ndata, bas_centers, bas_spec, bas_data, &
      nshells_aux, ndata_aux, bas_centers_aux, bas_spec_aux, bas_data_aux, &
      nbfs, nbfs_aux, df1c_tensor) bind(c, name="f_intor_df1c")
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
      ! ERI-2center-tensor
      integer(c_int32_t), value, intent(in) :: nbfs
      integer(c_int32_t), value, intent(in) :: nbfs_aux
      real(c_double), intent(in out) :: df1c_tensor(nbfs, nbfs_aux)
      ! Shells
      type(t_shells) :: shells, shells_aux

      shells = t_shells(nshells, bas_centers, bas_spec, bas_data)
      shells_aux = t_shells(nshells_aux, bas_centers_aux, bas_spec_aux, bas_data_aux)
      call intor_df1c(shells, shells_aux, df1c_tensor)
   end subroutine f_intor_df1c

end module mod_pa_std_intor_wrapper
