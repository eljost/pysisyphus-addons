module mod_pa_intor

   use mod_pa_constants, only: dp, i4
   use mod_pa_shells, only: t_shell, t_shells
   use mod_pa_init, only: init
   use mod_pa_int_2c2e, only: int_2c2e
   use mod_pa_int_4c2e, only: int_schwarz

   implicit none

   interface
      subroutine int2c_sub(La, Lb, ax, da, A, bx, db, B, R, res)
         import :: dp

         integer, intent(in) :: La, Lb
         real(dp), intent(in) :: ax(:), bx(:)
         real(dp), intent(in) :: da(:)
         real(dp), intent(in) :: db(:)
         real(dp), intent(in) :: A(3), B(3)
         real(dp), intent(in) :: R(3)
         !real(dp), intent(out), dimension(:) :: res
         real(dp), intent(out) :: res(:)
      end subroutine int2c_sub
   end interface

contains

   subroutine symmetrize_triu(nbfs, integrals)
      integer, intent(in) :: nbfs
      real(dp), intent(in out) :: integrals(:, :)
      integer :: i, j

      do i = 1, nbfs
         do j = i + 1, nbfs
            ! TODO: do it like in Psi4 hermitivitize()?
            integrals(j, i) = integrals(i, j)
         end do
      end do
   end subroutine symmetrize_triu

   function intor2d_sph_symmetric(shells, fncpntr, R) result(integrals)
      type(t_shells), intent(in) :: shells
      procedure(int2c_sub), pointer, intent(in) :: fncpntr
      real(dp), intent(in), optional :: R(3)

      integer(i4) :: a, b
      integer(i4) :: i, j, k
      real(dp), dimension(shells%nsphbfs, shells%nsphbfs) :: integrals
      real(dp), allocatable :: pres(:)
      type(t_shell) :: shell_a, shell_b

      call init()

      ! Initialize integral array with zeros
      integrals = 0
!$omp parallel do default(none) shared(shells, fncpntr, R, integrals) &
!$omp& private(shell_a, b, shell_b, pres, k, i, j) schedule(dynamic)
      ! Loop over shells
      do a = 1, shells%nshells
         shell_a = shells%shells(a)
         do b = a, shells%nshells
            shell_b = shells%shells(b)
            ! TODO: Only allocate once for largest L-pair?!
            ! Or allocate a 2d pres with one row per thread
            allocate (pres(shell_a%sph_size*shell_b%sph_size))

            ! Calculate integral for shell pair
            call fncpntr(shell_a%L, shell_b%L, &
                         shells%get_exps(a), shells%get_coeffs(a), shell_a%center, &
                         shells%get_exps(b), shells%get_coeffs(b), shell_b%center, &
                         R, pres)

            ! Set values/integrals from 1d array pres in 2d integrals array
            k = 1
            do j=shell_b%sph_index, shell_b%sph_index_end
               do i=shell_a%sph_index, shell_a%sph_index_end
                  integrals(i, j) = pres(k)
                  k = k + 1
               end do
            end do
            deallocate (pres)
         end do
      end do
!$omp end parallel do

      call symmetrize_triu(shells%nsphbfs, integrals)
   end function intor2d_sph_symmetric

   function int2c2e_sph(shells) result (integrals)
      type(t_shells), intent(in) :: shells
      real(dp) :: integrals(shells%nsphbfs, shells%nsphbfs)
      procedure(int2c_sub), pointer :: fncpntr

      fncpntr => int_2c2e
      integrals = intor2d_sph_symmetric(shells, fncpntr)
   end function int2c2e_sph

   function intor_schwarz(shells) result(schwarz_integrals)
      type(t_shells), intent(in) :: shells

      integer(i4) :: a, b
      integer(i4) :: i, j, k
      real(dp), dimension(shells%nsphbfs, shells%nsphbfs) :: schwarz_integrals
      real(dp), allocatable :: pres(:)
      type(t_shell) :: shell_a, shell_b

      call init()

      ! Initialize integral array with zeros
      schwarz_integrals = 0
!$omp parallel do default(none) shared(shells, schwarz_integrals) &
!$omp& private(shell_a, b, shell_b, pres, k, i, j) schedule(dynamic)
      ! Loop over shells
      do a = 1, shells%nshells
         shell_a = shells%shells(a)
         do b = a, shells%nshells
            shell_b = shells%shells(b)
            allocate (pres(shell_a%sph_size * shell_b%sph_size))

            ! Calculate integral for shell pair
            call int_schwarz(shell_a%L, shell_b%L, &
                         shells%get_exps(a), shells%get_coeffs(a), shell_a%center, &
                         shells%get_exps(b), shells%get_coeffs(b), shell_b%center, &
                         pres)

            ! Set values/integrals from 1d array pres in 2d integrals array
            k = 1
            do j=shell_b%sph_index, shell_b%sph_index_end
               do i=shell_a%sph_index, shell_a%sph_index_end
                  schwarz_integrals(i, j) = pres(k)
                  k = k + 1
               end do
            end do
            deallocate (pres)
         end do
      end do
!$omp end parallel do

      call symmetrize_triu(shells%nsphbfs, schwarz_integrals)
   end function intor_schwarz

end module mod_pa_intor
