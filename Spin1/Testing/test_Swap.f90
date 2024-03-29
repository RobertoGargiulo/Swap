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
  real (c_double), allocatable :: h(:), Jxy(:), H(:,:), v(:)
  real (c_double), allocatable :: E(:), W_r(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:), PH(:), W(:,:)

  logical :: SELECT
  EXTERNAL SELECT

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  dim = 3**nspin

  call buildHSwap(nspin, dim, H)

  allocate(E(dim), W_r(dim,dim), USwap(dim,dim))
  !H2 = H2 + identity(dim)
  call diagSYM( 'V', dim, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  
  !print*, "E = "
  !print "(*(F5.2,2X))", E(:)

  T1 = pi/2
  call expSYM( dim, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)
  !print*, "USwap = "
  !call printmat(dim, USwap, 'C')

  allocate(PH(dim),W(dim,dim))
  call diagUN( SELECT, dim, USwap, PH, W)

  print*, "arg(PH) = "
  print "(*(F5.2,2X))", atan2(dreal(PH(:)),dimag(PH(:)) )

  USwap = matmul(USwap,USwap)
  !print*, "USwap^2 = "
  !call printmat(dim, USwap, 'C')
  USwap = abs(USwap)
  !print*, "abs(USwap^2) = "
  !call printmat(dim, USwap, 'C')
  !if (all( abs(USwap) - identity(dim) <= 1.0e-6)) then
  !  print*, "USwap^2 = c * Id"
  !else
  !  print*, "USwap^2 != c * Id"
  !endif

  deallocate(PH,W,USwap)




end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
