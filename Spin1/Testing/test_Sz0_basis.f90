program prova

  use genmat
  use printing
  !use exponentiate
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)
  integer (c_int), parameter :: dimSpin1 = 3

  integer (c_int)     ::  nspin, dim, dim_Sz0
  integer (c_int), allocatable :: config(:), states(:), k_vec(:)
  integer (c_int) :: i, j, k, n, k1, k2, k3, Sz
  complex (c_double_complex), allocatable :: state(:)
  complex (c_double_complex) :: alpha(dimSpin1)

  logical :: SELECT
  EXTERNAL SELECT

  
  allocate(k_vec(3))
  do nspin = 2, 16, 2
    do k1 = 0, nspin/2
      k2 = k1
      k3 = nspin - (k1 + k2)
      k_vec(1) = k1
      k_vec(2) = k2
      k_vec(3) = k3
      !i = multinom(n,k_vec)
      !print*, "multinom(n,k) = ", multinom(nspin,k_vec), "nspin/k_vec = ", nspin, k_vec
    enddo
  enddo

  do nspin = 1, 12
    i = dimSpin1_Sz0(nspin)
    print *, "dim_Sz0 = ", i, nspin, i/3*3
  enddo


  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  dim = dimSpin1**nspin
  dim_Sz0 = dimSpin1_Sz0(nspin)

  allocate(states(dim_Sz0), config(nspin))
  call zero_mag_states(nspin, dim_Sz0, states)




end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
