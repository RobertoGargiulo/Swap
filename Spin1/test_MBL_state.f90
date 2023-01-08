program prova

  use genmat
  use printing
  use exponentiate
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)
  integer (c_int), parameter :: dimSpin1 = 3

  integer (c_int)     ::  nspin, dim
  integer (c_int), allocatable :: config(:), states(:)
  integer (c_int) :: i, j, k
  real (c_double) :: imp(3,3), T1
  real (c_double), allocatable :: Jxy(:), HS(:,:), H(:,:), Vz(:), hz(:)
  real (c_double), allocatable :: E(:), W_r(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:), PH(:), W(:,:), U(:,:)
  complex (c_double_complex), allocatable :: state(:)
  complex (c_double_complex) :: alpha(dimSpin1)

  logical :: SELECT
  EXTERNAL SELECT

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  dim = 3**nspin

  allocate(Vz(nspin-1), Jxy(nspin-1), hz(nspin))
  allocate(H(dim,dim))

  Vz(:) = 1
  Jxy(:) = 1
  hz(:) = 1
  call buildHMBL(nspin, dim, Jxy, Vz, hz, H)
  call printmat(dim, H, 'R')
  

  allocate(state(dim))
  alpha(1) = 1
  alpha(2) = 1
  alpha(3) = 1
  call buildProdState(nspin, dim, alpha, state)

  call printvec(dim, state, 'R')

  allocate(config(nspin))
  do i = 1, dim
    call decode(i-1,nspin,config)
    print *, state(i), config(:)
  enddo


end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
