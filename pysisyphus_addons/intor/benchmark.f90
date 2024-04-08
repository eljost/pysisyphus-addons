module mod_pa_intor_benchmark

   use mod_pa_constants, only: i4, dp
   use mod_pa_linalg, only: matrix_powerh
   use mod_pa_shells, only: t_shell, t_shells
   use mod_pa_init, only: init
   use mod_pa_screening, only: t_df_screener
   use mod_pa_intor, only: int2c2e_sph
   use mod_pa_int_3c2e, only: int_3c2e

   implicit none

contains

   subroutine benchmark_int3c2e(shells, shells_aux)
      ! AO basis
      type(t_shells), intent(in) :: shells
      ! Auxiliary AO basis for density fitting
      type(t_shells), intent(in) :: shells_aux

      ! Intermediate variables
      !
      ! Number of auxiliary basis functions
      integer :: naux
      ! Shell triplet for the density fitting 3-center-2-electron integral
      type(t_shell) :: shell_a, shell_b, shell_aux_c
      ! Index variables
      integer(i4) :: a, b, c
      ! 2-center-2-electron integrals (P|Q); later updated to contain its inverted square root (P|Q)^-0.5
      real(dp) :: metric(shells_aux%nsphbfs, shells_aux%nsphbfs)
      ! 1d-array holding the 3-center-2-electron integrals
      real(dp), allocatable :: integrals(:)
      ! Timing related
      real(dp) :: start_time, end_time, dur_time
      type(t_df_screener) :: screener
      real(dp) ::  schwarz_est, qvl_est, R_est, actual
      integer :: file_unit
      integer(i4) :: ntriplets, ntriplets_skipped
      logical :: bra_concentric, bra_small

      ntriplets = shells%nshells*(shells%nshells + 1)/2*shells_aux%nshells
      ntriplets_skipped = 0
      open (newunit=file_unit, file='int_estimates.bin', access="stream")

      ! Determine number of auxiliary basis functions, number of present densities and calculate
      ! number states from the number densities.
      naux = shells_aux%nsphbfs

      print '("Main basis:", i6, " shells, ", i6, " basis functions")', shells%nshells, shells%nsphbfs
      print '("Aux. basis:", i6, " shells, ", i6, " basis functions")', &
         shells_aux%nshells, shells_aux%nsphbfs

      ! Initialize integral pointers, if not already done
      call init()

      screener = t_df_screener(shells, shells_aux, 1d-4)
      call cpu_time(start_time)
      call screener%init()
      call cpu_time(end_time)
      dur_time = end_time - start_time
      print '("Initializing the DF-integral screener took ", f6.2, " s")', dur_time

      call cpu_time(start_time)
      ! Calculate density fitting metric (P|Q) using the auxilary basis and ...
      metric = int2c2e_sph(shells_aux)
      call cpu_time(end_time)
      dur_time = end_time - start_time
      print '("Calculation of density fitting metric (P|Q) took", f8.4, " s.")', end_time - start_time

      ! ... raise it to the power -0.5 to calculate (P|Q)^-0.5
      call cpu_time(start_time)
      call matrix_powerh(metric, -0.5d0)
      call cpu_time(end_time)
      dur_time = end_time - start_time
      print '("Calculation of (P|Q)**-0.5 took", f8.4, " s.")', end_time - start_time

      call cpu_time(start_time)
      ! Contract densities with 3-center-2-electron integrals to form gamma_P
      !
      ! Loop over all shells in the principal AO basis
      do a = 1, shells%nshells
         shell_a = shells%shells(a)
         do b = 1, a
            shell_b = shells%shells(b)
            ! TODO: the loop over the auxiliary shells should probably be the outermost loop,
            ! as this would allow for simple parallelization.
            !
            ! Loop over all shells in the auxiliary basis
            bra_concentric = norm2(shell_a%center - shell_b%center) <= 1d-6
            bra_small = screener%schwarz_bra(a*(a - 1)/2 + b) <= 1d-10
            if (bra_concentric .or. bra_small) then
               cycle
            endif

            do c = 1, shells_aux%nshells
               shell_aux_c = shells_aux%shells(c)

               ! Skip current shell triple when the integrals are estimated as negligible.
               call screener%debug_estimate(a, b, c, schwarz_est, qvl_est, R_est)

               if (qvl_est < 1e-10) then
                  ntriplets_skipped = ntriplets_skipped + 1
               end if

               ! TODO: allocate once for the biggest possible size and w/ additional dimensions for
               ! different threads?!
               allocate (integrals(shell_a%sph_size*shell_b%sph_size*shell_aux_c%sph_size))

               ! Calculate integrals over current shell triple and store them in the array 'integrals'
               call int_3c2e(shell_a%L, shell_b%L, shell_aux_c%L, &
                             shells%get_exps(a), shells%get_coeffs(a), shell_a%center, &
                             shells%get_exps(b), shells%get_coeffs(b), shell_b%center, &
                             shells_aux%get_exps(c), shells_aux%get_coeffs(c), shell_aux_c%center, &
                             integrals)

               actual = norm2(integrals)
               write (file_unit) schwarz_est, qvl_est, actual, R_est

               deallocate (integrals)
            end do
         end do
      end do
      ! End of density contraction
      call cpu_time(end_time)
      dur_time = end_time - start_time
      print '("Density contraction took", f8.4, " s.")', end_time - start_time
      print '("Skipped ", i12, " of ", i12, " integral triplets (", f8.2, "%)")', &
         ntriplets_skipped, ntriplets, 100*real(ntriplets_skipped, dp)/real(ntriplets, dp)

      close (file_unit)

   end subroutine benchmark_int3c2e
end module mod_pa_intor_benchmark
