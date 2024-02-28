module mod_pa_shells

   use iso_fortran_env, only: error_unit, int32

   use mod_pa_constants, only: dp
   use mod_pa_normalization, only: norm_cgto_coeffs

   implicit none

   type :: t_shell
      ! Index of shell center.
      integer :: center_ind
      real(dp), dimension(3) :: center
      ! Atomic number.
      integer :: atomic_num
      ! Total angular momentum.
      integer :: L
      ! Number of primitive Gaussians contributing to this contracted shell.
      integer :: nprims
      ! Pointer, pointing to coefficients and exponents in the respective arrays in t_shells.
      integer :: pntr
      integer :: pntr_end
      ! Number of Cartesian basis functions in the shell
      integer :: cart_size
      ! Starting and ending indices of the shell in Cartesian integral matrices.
      integer :: cart_index
      integer :: cart_index_end
      ! Number of spherical basis functions in the shell
      integer :: sph_size
      ! Starting and ending indices of the shell in spherical integral matrices.
      integer :: sph_index
      integer :: sph_index_end

   end type t_shell

   interface t_shell
      module procedure :: t_shell_constructor
   end interface t_shell

   ! Container, representing multiple contracted Gaussian shells.
   type t_shells
      type(t_shell), allocatable :: shells(:)
      real(dp), allocatable :: exponents(:)
      real(dp), allocatable :: coefficients(:)
      ! Number of t_shell objects.
      integer :: nshells
      ! Number of Cartesian basis functions.
      integer :: ncartbfs
      ! Number of spherical basis functions.
      integer :: nsphbfs

      contains
      procedure :: print_shell
      procedure :: print_shells
      procedure :: get_coeffs
      procedure :: get_exps
   end type t_shells

   interface t_shells
      module procedure :: t_shells_constructor
   end interface t_shells

contains
   ! Custom t_shell constructor, that calculates several values, so they don't have to be
   ! supplied explicitly.
   type(t_shell) function t_shell_constructor( &
      center_ind, center, atomic_num, L, nprims, pntr, cart_index, sph_index) result(res)

      integer, intent(in) :: center_ind
      real(dp), dimension(3), intent(in) :: center
      integer, intent(in) :: atomic_num
      integer, intent(in) :: L
      integer, intent(in) :: nprims
      integer, intent(in) :: pntr
      integer, intent(in) :: cart_index
      integer, intent(in) :: sph_index

      res%center_ind = center_ind
      res%center = center
      res%atomic_num = atomic_num
      res%L = L
      res%nprims = nprims
      res%pntr = pntr
      res%pntr_end = pntr + nprims - 1
      res%cart_size = (L + 2)*(L + 1)/2
      res%sph_size = 2*L + 1
      res%cart_index = cart_index
      res%cart_index_end = cart_index + res%cart_size - 1
      res%sph_index = sph_index
      res%sph_index_end = sph_index + res%sph_size - 1
   end function t_shell_constructor

   function t_shells_constructor(nshells, bas_centers, bas_spec, bas_data) result(shells)
      ! Number of shells
      integer(int32), intent(in) :: nshells
      ! bas_centers
      ! One entry per shell, integer array.
      !   center_ind, atomic number, pointer to center coordinates in bas_data (3 integers)
      integer(int32), intent(in) ::  bas_centers(nshells, 3)
      ! bas_spec
      ! One entry per shell, integer array.
      !   shell_ind, total angmom, N_pgto, N_cgto, \
      !   pointer to contraction coefficients and exponents in bas_data \
      !   (2*N_pgto floats)
      integer(int32), intent(in) ::  bas_spec(nshells, 5)
      real(dp), intent(in) :: bas_data(:)

      type(t_shell) :: shell
      type(t_shells) :: shells
      ! Number of Cartesian and spherical basis functions
      integer :: ncartbfs, nsphbfs
      ! Cartesian & spherical starting index of the respective shell
      integer :: cart_index, sph_index
      ! Cartesian center of a shell
      real(kind=dp) :: center(3)
      ! Center index, atomic number and total angular momentum of a shell
      integer :: center_ind, atomic_num, L
      ! cur_pntr is read from bas_centers & bas_spec and used to access bas_data
      integer :: cur_pntr
      ! prim_pntr is incremented by the number of primitives after each shell is read
      ! and points to contraction coefficients and exponents in the shells object.
      integer :: prim_pntr, prim_pntr_end
      ! Number of primitives
      integer :: nprims
      integer :: i

      ! Count total number of primitives
      nprims = 0
      do i = 1, nshells
         nprims = nprims + bas_spec(i, 3)
      end do

      allocate (shells%shells(nshells))
      allocate (shells%coefficients(nprims))
      allocate (shells%exponents(nprims))

      prim_pntr = 1
      cart_index = 1
      sph_index = 1
      ncartbfs = 0
      nsphbfs = 0
      do i = 1, nshells
         center_ind = bas_centers(i, 1)
         atomic_num = bas_centers(i, 2)
         cur_pntr = bas_centers(i, 3)
         center = bas_data(cur_pntr:cur_pntr+2)

         L = bas_spec(i, 2)
         nprims = bas_spec(i, 3)
         cur_pntr = bas_spec(i, 5)
         ! Create new shell ...
         shell = t_shell(center_ind, center, atomic_num, L, nprims, &
                         prim_pntr, cart_index, sph_index)
         ! ... and store it
         shells%shells(i) = shell
         prim_pntr_end = prim_pntr + nprims - 1
         ! Read contraction coefficients and orbital exponents
         shells%coefficients(prim_pntr:prim_pntr_end) = bas_data(cur_pntr:cur_pntr+nprims-1)
         cur_pntr = cur_pntr + nprims
         shells%exponents(prim_pntr:prim_pntr_end) = bas_data(cur_pntr:cur_pntr+nprims-1)
         ! Apply coefficient & exponent dependent normalization to coefficients.
         ! Angular momentum dependent (lmn) dependent normalization factors are taken
         ! into account in the integral library.
         call norm_cgto_coeffs(shells%coefficients(prim_pntr:prim_pntr_end), &
                               shells%exponents(prim_pntr:prim_pntr_end), L)
         ! Increment pointers & indices
         prim_pntr = prim_pntr + nprims
         ncartbfs = ncartbfs + shell%cart_size
         nsphbfs = nsphbfs + shell%sph_size
         cart_index = cart_index + shell%cart_size
         sph_index = sph_index + shell%sph_size
      end do

      ! Set remainig attributes on the shells derived type
      shells%nshells = nshells
      shells%ncartbfs = ncartbfs
      shells%nsphbfs = nsphbfs
   end function t_shells_constructor

   ! Helper procedure to print all shells in t_shells
   subroutine print_shells(self)
      class(t_shells), intent(in) :: self
      integer :: i

      do i = 1, self%nshells
        call self%print_shell(i)
      end do
   end subroutine print_shells

   ! Helper procedure to print single shell t_shell
   subroutine print_shell(self, i)
      class(t_shells), intent(in) :: self
      integer, intent(in) :: i
      type(t_shell) :: shell

      shell = self%shells(i)
      print *, "Shell", i
      print '(A,A,I2,A,I5,A,F12.6,F12.6,F12.6)', char(9), "L=", shell%L, &
          ", nprims=", shell%nprims, ", at", shell%center
      print "(A,*(F12.6,:',',1x))", "coeffs", self%coefficients(shell%pntr:shell%pntr_end)
      print "(A,*(F14.6,:',',1x))", "exps", self%exponents(shell%pntr:shell%pntr_end)
   end subroutine print_shell

   function get_coeffs(self, i) result(coeffs)
      class(t_shells), intent(in) :: self
      integer, intent(in) :: i
      real(dp) :: coeffs(self%shells(i)%nprims)

      coeffs = self%coefficients(self%shells(i)%pntr:self%shells(i)%pntr_end)
   end function get_coeffs

   function get_exps(self, i) result(exps)
      class(t_shells), intent(in) :: self
      integer, intent(in) :: i
      real(dp) :: exps(self%shells(i)%nprims)

      exps = self%exponents(self%shells(i)%pntr:self%shells(i)%pntr_end)
   end function get_exps

end module mod_pa_shells
