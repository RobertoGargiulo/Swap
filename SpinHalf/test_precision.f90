program prova

  use iso_fortran_env
  use iso_c_binding
  implicit none

  integer :: i
  integer (c_int) :: i_int   
  integer (c_short) :: i_short
  integer (c_long) :: i_long
  integer (c_long_long) :: i_long_long
  integer (c_size_t) :: i_size_t
  integer (c_int8_t) :: i_int8_t
  integer (c_int16_t) :: i_int16_t
  integer (c_int32_t) :: i_int32_t
  integer (c_int64_t) :: i_int64_t
  integer (c_int128_t) :: i_int128_t
  integer (kind=int8) :: i8


  real (c_double) :: r_double
  real (c_float) :: r_float
  real (c_long_double) :: r_long_double
  real (c_float128) :: r_float128


  do i = 0, 4
    i_size_t = 2**(8*2**i-1)-1
    print "(A,I0,A,I0)", "2**(", 8*2**i-1, ") - 1 = ", i_size_t
  enddo

  i_size_t = 1
  do i = 1, 100
    i_size_t = i_size_t * 2 !_c_size_t
    print *, i, i_size_t-1
  enddo

  print *, "Maximum values and errors (epsilon): "
  print "(A18,I0)", "max_c_int = ", huge(i_int)
  print "(A18,I0)", "max_c_short = ", huge(i_short)
  print "(A18,I0)", "max_c_long = ", huge(i_long)
  print "(A18,I0)", "max_c_long_long = ", huge(i_long_long)
  print "(A18,I0)", "max_c_size_t = ", huge(i_size_t)
  print "(A18,I0)", "max_c_int8_t = ", huge(i_int8_t)
  print "(A18,I0)", "max_c_int64_t = ", huge(i_int64_t)
  print "(A18,I0)", "max_c_int128_t = ", huge(i_int128_t)


  print *, "Maximum, minimum values and errors (epsilon): "
  print *, "c_double: max = ", huge(r_double), "min = ", tiny(r_double), &
    & "eps = ", epsilon(r_double)
  print *,  "log(max) = ", log(huge(r_double)), "log(min) = ", log(tiny(r_double)), "log(eps) = ", log(epsilon(r_double))

  r_double = 0
  do i = 1, 1000000
    r_double = 1*r_double + epsilon(r_double)
  enddo
  r_double = r_double / 1000000

  print *, r_double, log(r_double), log(r_double) - log(r_double)
  print *, "c_float: max = ", huge(r_float), "min = ", tiny(r_float), &
    & "eps = ", epsilon(r_float)
  print *,  "log(max) = ", log(huge(r_float)), "log(min) = ", log(tiny(r_float)), "log(eps) = ", log(epsilon(r_float))

  print *, "max_c_long_double = ", huge(r_long_double), "min_c_long_double = ", tiny(r_long_double), &
    & "eps_c_long_double = ", epsilon(r_long_double)
  print *, "max_c_float128 = ", huge(r_float128), "min_c_float128 = ", tiny(r_float128), &
    & "eps_c_float128 = ", epsilon(r_float128)

  r_double = 0
  print *, "log(0) = ", log(r_double)
  
  print *, "pi: "
  print *, 4*atan(1.0), 4*datan(1.d0), 4.d0*atan(1.d0), 4.d0*datan(1.d0), &
    & 4._c_double*atan(1._c_size_t), 4._c_size_t*datan(1._c_size_t)

  print *, "Complex precisions: "
  print *, "c_double = ", cmplx(1,0,c_double), cmplx(1,c_double), cmplx(1,kind=c_double)

  print *, "c_double_complex = ", cmplx(1,0,c_double_complex), &
    & cmplx(1,c_double_complex), cmplx(1,kind=c_double_complex)

  print *, "c_float_complex = ", cmplx(1,0,c_float_complex), & 
    & cmplx(1,c_float_complex), cmplx(1,kind=c_float_complex)

  print *, "c_long_double_complex = ", cmplx(1,0,c_long_double_complex), &
    & cmplx(1,c_long_double_complex), cmplx(1,kind=c_long_double_complex)

  99 format(A,*(F26.20,SP,F26.20,'i'))


end program
