program prova

  use genmat
  use printing
  use exponentiate
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim
  integer (c_int), allocatable :: config(:), states
  integer (c_int) :: i, j, k
  real (c_double) :: imp(3,3), T1
  real (c_double), allocatable :: Jxy(:), HS(:,:), H(:,:), Vz(:), hz(:)
  real (c_double), allocatable :: E(:), W_r(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:), PH(:), W(:,:), U(:,:)

  logical :: SELECT
  EXTERNAL SELECT

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  dim = 3**nspin

!  allocate(HS(dim,dim))
!
!  call buildHSwap(nspin, dim, HS)
!
!  allocate(E(dim), W_r(dim,dim), USwap(dim,dim))
!  !H2 = H2 + identity(dim)
!  call diagSYM( 'V', dim, HS, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
!  deallocate(HS)
!
!  T1 = pi/2
!  call expSYM( dim, -C_UNIT*T1, E, W_r, USwap )
!  deallocate(E, W_r)


  allocate(Vz(nspin-1), Jxy(nspin-1), hz(nspin))
  allocate(H(dim,dim))

  Vz(:) = 1
  Jxy(:) = 1
  hz(:) = 1
  call buildHMBL(nspin, dim, Jxy, Vz, hz, H)
  call printmat(dim, H, 'R')
  



end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
