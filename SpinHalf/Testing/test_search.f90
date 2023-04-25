program test

  use sorts, only: sort => dpquicksort
  use functions, only: search => binsearch_rightmost
  use iso_c_binding, only: dp => c_double, ip => c_int
  implicit none

  real (dp), parameter :: pi = 4*datan(1.0_dp)
  real (dp), allocatable :: array(:)
  real (dp) :: val
  integer (ip) :: indx, i, n

  write (*,*) "Size of array:"
  read (*,*) n
  print *, ""

  allocate(array(n))

  call random_number(array)
  call random_number(val)
  array = 2*pi*(array - 0.5)
  val = 2*pi*(val-0.5)
  call sort(array)
  
  print "(*(A20))", "i", "array(i)", "val"
  do i = 1, size(array)
    print *, i, array(i), val
  enddo

  indx = search(val, array)

  print "(*(A20))", "indx", "array(indx)"
  print *, indx, array(indx)

  print *, "Check: "
  print "(*(A20))", "i", "array(i) - val"
  do i = 1, size(array)
    print *, i, array(i) - val
  enddo

end program
