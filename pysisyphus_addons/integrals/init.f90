module mod_pa_init

   use mod_pa_boys, only: boys_init
   use mod_pa_int_2c2e, only: int_2c2e_init
   use mod_pa_int_3c2e, only: int_3c2e_init

   implicit none

   logical :: initialized = .false.

contains
   subroutine init() bind(c, name="f_init")
      if (initialized) then
         print *, "Already initialized!"
         return
      end if


      call boys_init()
      call int_2c2e_init()
      call int_3c2e_init()
      initialized = .true.
      print *, "Initialized integral modules."

   end subroutine init

   end module mod_pa_init
