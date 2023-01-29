program test_Sz

  use functions, only: basis_Sz, dimSpin1_Sz
  use iso_c_binding, only: ip => c_int, dp => c_double, cp => c_double_complex
  implicit none

  integer (ip), parameter :: dimSpin1 = 3
  integer (ip) :: nspin, dim_Sz, Sz, dim
  integer (ip), allocatable :: states(:)

  write (*,*) "nspin = "
  read (*,*) nspin
  print *, ""

  dim = 0
  do Sz = -nspin, nspin
    dim_Sz = dimSpin1_Sz(nspin, Sz)
    dim = dim + dim_Sz
    print *, "Sz = ", Sz, ", dim_Sz = ", dim_Sz
    allocate(states(dim_Sz))
    call basis_Sz(nspin, dim_Sz, Sz, states)
    deallocate(states)
    print *, ""
  enddo
  print *, "dim = ", dim, ", dimSpin1**nspin = ", dimSpin1**nspin


end program
