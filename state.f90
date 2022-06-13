program pstate

  use genmat
  use printing
  use exponentiate
  use iso_c_binding

  implicit none

  complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real(c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, i
  complex (c_double_complex) :: alpha, beta

  complex (c_double_complex), dimension(:), allocatable :: state
  integer (c_int), dimension(:), allocatable :: config


  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin

  allocate(state(dim), config(nspin))

  write (*,*) "alpha"
  read (*,*) alpha
  print*,""

  write (*,*) "beta"
  read (*,*) beta
  print*,""

  call buildStaggState(nspin, dim, alpha, beta, state)

  call printvec(dim, state, 'C')

 print *, mag_stag_z(nspin, dim, state) 

end program pstate
