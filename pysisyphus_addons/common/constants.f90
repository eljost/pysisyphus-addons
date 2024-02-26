module mod_pa_constants
   use iso_fortran_env, only: compiler_options

   ! Double precision kind
   integer, parameter :: dp = selected_real_kind(15, 307)
   ! 4 byte integer
   integer, parameter :: i4 = selected_int_kind(9)
   !integer, parameter :: i8 = selected_int_kind(15)

   ! TODO: check if d0 at the end of the literals matters ... I guess it shouldn't?!
   real(dp), parameter :: PI = 4*atan(1.0d0)
   real(dp), parameter :: PI_4 = PI/4d0
   logical :: with_openmp = index(compiler_options(),'openmp') /= 0
end module mod_pa_constants
