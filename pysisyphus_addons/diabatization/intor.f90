! [1]    https://doi.org/10.1039/B204199P
!        A fully direct RI-HF algorithm: Implementation, optimised auxiliary
!        basis sets, demonstration of accuracy and eﬃciency
!        Weigend, 2002
! [2]    https://doi.org/10.1021/jp101235a
!        Predicting Accurate Electronic Excitation Transfer Rates via Marcus Theory
!        with Boys or Edmiston−Ruedenberg Localized Diabatization
!        Subotnik, Vura-Weis, Sodt, Ratner, 2010

module mod_pa_diabatization
   use mod_pa_constants, only: i4, dp
   use mod_pa_linalg, only: matrix_powerh
   use mod_pa_shells, only: t_shell, t_shells
   use mod_pa_init, only: init
   use mod_pa_screening, only: t_df_screener
   use mod_pa_intor, only: int2c2e_sph
   use mod_pa_int_3c2e, only: int_3c2e

   implicit none

contains

   ! Contract (transition) densitites in the AO basis w/ the 4d-Coulomb tensor using
   ! density fitting.
   subroutine contract_coulomb_densities_2d(shells, shells_aux, densities, df_tensor)
      ! AO basis
      type(t_shells), intent(in) :: shells
      ! Auxiliary AO basis for density fitting
      type(t_shells), intent(in) :: shells_aux
      ! (Multiple) densities in the AO basis; shape is (ndens, nao, nao)
      real(dp), intent(in) :: densities(:, :, :)
      ! 2d-array holding the contraction gamma_P with the inverse square root of the metric
      ! shape(naux, ndens)
      ! work corresponds to B^JK_Y (eq. 33) in [2].
      real(dp), intent(in out):: df_tensor(:, :)

      ! Intermediate variables
      !
      ! Number of densities, number of auxiliary basis functions
      integer :: ndens, naux
      ! Shell triplet for the density fitting 3-center-2-electron integral
      type(t_shell) :: shell_a, shell_b, shell_aux_c
      ! Index variables
      integer(i4) :: a, b, c
      integer(i4) :: i, j, k
      integer(i4) :: dens_ind, int_ind
      ! 2-center-2-electron integrals (P|Q); later updated to contain its inverted square root (P|Q)^-0.5
      real(dp) :: metric(shells_aux%nsphbfs, shells_aux%nsphbfs)
      ! 1d-array holding the 3-center-2-electron integrals
      real(dp), allocatable :: integrals(:)
      integer(i4) :: size_abc
      ! 2d-array holding the contraction of the densities with the 3-center-2-electron integrals.
      ! shape(naux, ndens)
      ! gamma_P corresponds to R^JK_X (eq. 32) in [2].
      real(dp), allocatable :: gamma_P(:, :)
      ! Timing related
      real(dp) :: start_time, end_time, dur_time
      ! Screener
      type(t_df_screener) :: screener
      ! Norm-estimate for integral batch
      real(dp) ::  int_estimate
      integer(i4) :: ntriplets, ntriplets_skipped

      ntriplets = shells%nshells * (shells%nshells + 1) / 2 * shells_aux%nshells
      ntriplets_skipped = 0

      ! Determine number of auxiliary basis functions, number of present densities and calculate
      ! number states from the number densities.
      naux = shells_aux%nsphbfs
      ndens = size(densities, dim=1)

      print '("Main basis:", i6, " shells, ", i6, " basis functions")', shells%nshells, shells%nsphbfs
      print '("Aux. basis:", i6, " shells, ", i6, " basis functions")', &
         shells_aux%nshells, shells_aux%nsphbfs
      print '(" Densities:", i6)', ndens

      allocate (gamma_P(naux, ndens))
      allocate (integrals((2 * shells%L_max + 1)**2 * (2*shells_aux%L_max + 1)))

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

      ! Initialize matrices
      gamma_P = 0
      df_tensor = 0

      call cpu_time(start_time)
      ! Contract densities with 3-center-2-electron integrals to form gamma_P
      !
      ! Loop over all shells in the auxiliary basis
      do c = 1, shells_aux%nshells
         shell_aux_c = shells_aux%shells(c)
         ! Loop over all shells in the principal AO basis
         do a = 1, shells%nshells
            shell_a = shells%shells(a)
            do b = 1, a
               shell_b = shells%shells(b)
               size_abc = shell_a%sph_size * shell_b%sph_size * shell_aux_c%sph_size

               ! Skip current shell triple when the integrals are estimated as negligible.
               !int_estimate = screener%qvl_estimate(a, b, c)
               int_estimate = screener%schwarz_estimate(a, b, c)

               if (int_estimate < 1e-10) then
                  ntriplets_skipped = ntriplets_skipped + 1
                  cycle
               endif

               ! Calculate integrals over current shell triple and store them in the array 'integrals'
               call int_3c2e(shell_a%L, shell_b%L, shell_aux_c%L, &
                             shells%get_exps(a), shells%get_coeffs(a), shell_a%center, &
                             shells%get_exps(b), shells%get_coeffs(b), shell_b%center, &
                             shells_aux%get_exps(c), shells_aux%get_coeffs(c), shell_aux_c%center, &
                             integrals(1:size_abc))

               ! Contract densities with integrals in a direct fashion
               !
               ! Loop over all present densities
               do dens_ind = 1, ndens
                  int_ind = 1
                  ! Loop over basis functions
                  do i = shell_a%sph_index, shell_a%sph_index_end
                     do j = shell_b%sph_index, shell_b%sph_index_end
                        ! Loop over auxiliary basis functions
                        do k = shell_aux_c%sph_index, shell_aux_c%sph_index_end
                           gamma_P(k, dens_ind) = gamma_P(k, dens_ind) + integrals(int_ind) &
                                                   *densities(dens_ind, i, j)
                           ! Take symmetry of the density matrix into account; off-diagonal
                           ! elements are counted twice.
                           if (a .ne. b) then
                              gamma_P(k, dens_ind) = gamma_P(k, dens_ind) + integrals(int_ind) &
                                                      *densities(dens_ind, j, i)
                           end if
                           int_ind = int_ind + 1
                        end do
                        ! End loop over auxiliary basis functions
                     end do
                  end do
                  ! End loops over basis functions
               end do
               ! End loops over densities
            end do
         end do
      end do
      ! End of density contraction

      deallocate (integrals)

      call cpu_time(end_time)
      dur_time = end_time - start_time
      print '("Density contraction took", f8.4, " s.")', end_time - start_time
      print '("Screened out ", i12, " of ", i12, " integral triplets (", f8.2, "%)")', &
         ntriplets_skipped, ntriplets, 100 * real(ntriplets_skipped, dp) / real(ntriplets, dp)

      ! Contract gamma_P with inverse square root (P|Q)^-1/2 of the metric (P|Q)
      ! The result is stored in 'df_tensor'.
      ! This contraction corresponds to eq. (33) in [2]
      call dgemm("N", "N", naux, ndens, naux, 1d0, metric, naux, gamma_P, naux, &
                 0d0, df_tensor, naux)

   end subroutine contract_coulomb_densities_2d

   ! Contract (transition) densitites in the AO basis w/ the 4d-Coulomb tensor using
   ! density fitting.
   subroutine contract_coulomb_densities_4d(shells, shells_aux, densities, coulomb_tensor)
      ! AO basis
      type(t_shells), intent(in) :: shells
      ! Auxiliary AO basis for density fitting
      type(t_shells), intent(in) :: shells_aux
      ! (Multiple) densities in the AO basis; shape is (ndens, nao, nao)
      real(dp), intent(in) :: densities(:, :, :)
      ! 4d Coulomb-tensor; shape (nstates, nstates, nstates, nstates)
      ! coulomb_tensor corresponds to R_JKLM (eq. 34) in [2]
      real(dp), intent(in out) :: coulomb_tensor(:, :, :, :)

      ! Intermediate variables
      !
      ! Number of densities, number of states, number of auxiliary basis functions
      integer :: ndens, nstates, naux
      integer(i4) :: i, j, k, l, l_max, m
      integer(i4) :: dens_ind1, dens_ind2
      ! 2-center-2-electron integrals (P|Q); later updated to contain its inverted square root (P|Q)^-0.5
      ! work corresponds to B^JK_Y (eq. 33) in [2].
      real(dp), allocatable :: df_tensor(:, :)
      ! Element of the Coulomb tensor
      real(dp) :: val

      ! Determine number of auxiliary basis functions, number of present densities and calculate
      ! number states from the number densities.
      naux = shells_aux%nsphbfs
      ndens = size(densities, dim=1)
      nstates = int(sqrt(8*real(ndens, dp) + 1)/2d0 - 0.5d0)

      allocate(df_tensor(naux, ndens))
      call contract_coulomb_densities_2d(shells, shells_aux, densities, df_tensor)

      ! Loop over work array and construct 4d Coulomb-tensor from it.
      ! Take 8-fold symmetry into account.
      ! Corresponds to eq. (34) in [2].
      dens_ind1 = 1
      do i = 1, nstates
         do j = 1, i
            dens_ind2 = 1
            do k = 1, i
               if (k == i) then
                  l_max = j
               else
                  l_max = k
               end if
               do l = 1, l_max
                  val = 0
                  do m = 1, naux
                     val = val + df_tensor(m, dens_ind1)*df_tensor(m, dens_ind2)
                  end do
                  coulomb_tensor(i, j, k, l) = val
                  coulomb_tensor(j, i, k, l) = val
                  coulomb_tensor(i, j, l, k) = val
                  coulomb_tensor(j, i, l, k) = val
                  coulomb_tensor(k, l, i, j) = val
                  coulomb_tensor(k, l, j, i) = val
                  coulomb_tensor(l, k, i, j) = val
                  coulomb_tensor(l, k, j, i) = val
                  dens_ind2 = dens_ind2 + 1
               end do
            end do
            dens_ind1 = dens_ind1 + 1
         end do
      end do
   end subroutine contract_coulomb_densities_4d
end module mod_pa_diabatization
