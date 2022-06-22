program prova
  
  use genmat
  use exponentiate
  use iso_c_binding
  implicit none

  integer (c_int) :: nspin, dim, i, j, k, n
  integer (c_int), dimension(:), allocatable :: states

  print *, "nspin ="
  read (*,*) nspin
  dim = 2**nspin
  allocate(states(dim))

  call zero_mag_states(nspin, dim, states)

  n= product((/(i,i=1,nspin)/)) / (product((/(i,i=1,nspin/2)/)))**2

  do i = 1, n
    print *, states(i), i
  enddo





end program prova
