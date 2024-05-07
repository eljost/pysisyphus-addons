module mod_pa_std_intor
   use mod_pa_constants, only: dp, i4
   use mod_pa_helpers, only: print_array
   use mod_pa_shells, only: t_shell, t_shells
   use mod_pa_init, only: init
   use mod_pa_linalg, only: inv_chol
   use mod_pa_int_3c2e, only: int_3c2e
   use mod_pa_intor, only: int2c2e_sph
   use mod_pa_screening, only: t_df_screener

   implicit none

contains

   subroutine intor_df1c(shells, shells_aux, df1c_tensor)

      ! Principal basis and auxiliary basis
      type(t_shells), intent(in) :: shells, shells_aux
      ! 2d-array w/ shape (nints, naux) holding the (aa|C) integrals.
      real(dp) :: df1c_tensor(shells%nsphbfs, shells_aux%nsphbfs)

      ! Intermediate variables
      type(t_shell) :: shell_a, shell_aux_c
      ! Screener for density fitting integrals
      type(t_df_screener) :: screener
      ! 1d-array holding a batch of DF-integrals (aa|C)
      real(dp), allocatable :: integral_batch(:)
      ! nints: Total number of products in the bra of (aa|C)
      ! naux: Number of auxiliary basis functions
      integer(i4) :: nbra, naux
      ! Density fitting metric; will finally contain (P|Q)^(-0.5) and is absorbed
      ! into 'integrals'.
      real(dp), allocatable :: metric(:, :)
      ! Temporary variable holding the estimate of an integral batch; used for screening.
      real(dp) :: int_estimate
      integer :: a, c, i, j, k, batch_ind, a_pointer
      ! Number of skipped integral batches.
      integer(i4) :: skipped
      ! Temporary variables holding (squared) sizes of shell a
      integer(i4) :: size_a, size2_a
      ! Size of integral batch (aa|c)
      integer(i4) :: batch_size

      ! Initialize integral procedures
      call init()

      ! Calculate number of basis function products in bra of (aa|C)
      nbra = shells%nsphbfs
      naux = shells_aux%nsphbfs

      ! Allocate 1d-array holding integral batches once w/ maximum required size.
      allocate (integral_batch((2*shells%L_max + 1)**2*(2*shells_aux%L_max + 1)))

      ! Calculate density fitting metric, do a Cholesky decomposition and invert
      ! the lower triangular matrix.
      metric = int2c2e_sph(shells_aux)
      call inv_chol(metric)

      ! Construct & initialize screener for DF integrals
      screener = t_df_screener(shells, shells_aux, 1d-4)
      call screener%init()

      skipped = 0
      a_pointer = 1

      do a = 1, shells%nshells
         shell_a = shells%shells(a)
         size_a = shell_a%sph_size
         size2_a = size_a**2
         do c = 1, shells_aux%nshells
            shell_aux_c = shells_aux%shells(c)
            ! Size of current integral batch (aa|C)
            batch_size = size2_a*shell_aux_c%sph_size

            ! Skip current shell triple when the integrals are estimated as negligible.
            int_estimate = screener%schwarz_estimate(a, a, c)

            if (int_estimate < 1e-10) then
               skipped = skipped + 1
               cycle
            end if

            ! Calculate integrals over current shell pair/triple (aa|C)
            call int_3c2e(shell_a%L, shell_a%L, shell_aux_c%L, &
                          shells%get_exps(a), shells%get_coeffs(a), shell_a%center, &
                          shells%get_exps(a), shells%get_coeffs(a), shell_a%center, &
                          shells_aux%get_exps(c), shells_aux%get_coeffs(c), shell_aux_c%center, &
                          integral_batch(1:batch_size))

            ! Put integral batch at the appropriate places in 'integrals'
            batch_ind = 0
            do i = 1, shell_a%sph_size
               do j = 1, shell_a%sph_size
                  ! Loop over auxiliary basis functions
                  do k = shell_aux_c%sph_index, shell_aux_c%sph_index_end
                     batch_ind = batch_ind + 1
                     if (i /= j) then
                        cycle
                     end if
                     df1c_tensor(shell_a%sph_index+i-1, k) = integral_batch(batch_ind)
                  end do
                  ! End loop over auxiliary basis functions
               end do
            end do
         end do
         ! End loop over auxiliary shells
         
         ! Increase a_pointer after we iterated over all auxiliary basis functions for a given a
      end do
      ! End loop over principal shells

      ! Absorb DF metric (P|Q)^-0.5 into DF tensor (aa|C).
      ! Multiply integrals from the right w/ transposed lower triangular matrix 'metric'
      call dtrmm("R", "L", "T", "N", nbra, naux, 1d0, metric, naux, df1c_tensor, nbra)

   end subroutine intor_df1c

   subroutine intor_eri2c(shells, shells_aux, eri2c_tensor)
      ! Calculates (aa|bb) ERIs w/ density fitting and stores them in a 2d array.
      ! required for the eXact simplified TD-DFT.

      ! Principal basis and auxiliary basis
      type(t_shells), intent(in) :: shells, shells_aux
      ! 2d-array holding the ERIs
      real(dp), intent(in out) :: eri2c_tensor(:, :)
      real(dp) :: df1c_tensor(shells%nsphbfs, shells_aux%nsphbfs)
      integer :: nbfs, naux

      nbfs = shells%nsphbfs
      naux = shells_aux%nsphbfs

      call intor_df1c(shells, shells_aux, df1c_tensor)

      ! Construct "full"  ERI-tensor (aa|bb) from integrals.T @ integrals
      call dgemm("N", "T", nbfs, nbfs, naux, 1d0, df1c_tensor, nbfs, df1c_tensor, nbfs, 0d0, &
                 eri2c_tensor, nbfs)
   end subroutine intor_eri2c
end module mod_pa_std_intor
