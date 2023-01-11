program prova

  use genmat
  use printing
  use exponentiate
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)
  integer (c_int), parameter :: dimSpin1 = 3

  integer (c_int)     ::  nspin, dim, dim_Sz0
  integer (c_int), allocatable :: config(:), states(:), k_vec(:)
  integer (c_int) :: i, j, k, n, k1, k2, k3, Sz
  complex (c_double_complex), allocatable :: state(:)
  complex (c_double_complex) :: alpha(dimSpin1)
  real (c_double), allocatable :: H(:,:), hz(:), Jxy(:), Vz(:)
  real (c_double), allocatable :: E(:), W_r(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:), PH(:), W(:,:)

  integer(c_int) :: count_beginning, count_end, count_rate

  logical :: SELECT
  EXTERNAL SELECT


  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  call system_clock(count_beginning, count_rate)

  dim = dimSpin1**nspin
  dim_Sz0 = dimSpin1_Sz0(nspin)


  !----------- Full Hilbert Space --------------------------!

  allocate(H(dim,dim))
  !allocate(H(dim,dim), hz(nspin), Jxy(nspin-1), Vz(nspin-1))
  !call random_number(hz)
  !hz = 2*hz
  !call random_number(Jxy)
  !Jxy = 2*(Jxy-0.5)
  !Vz = 0
  !call buildHMBL(nspin, dim, Jxy, Vz, hz, H)
  !print *, "H = "
  !call printmat(dim, H, 'R')

  call buildHSwap(nspin, dim, H)
  print *, "H = "
  call printnzmat(dim, H, 'R')
  print*, ""

  allocate(E(dim), W_r(dim,dim), USwap(dim,dim))
  !H2 = H2 + identity(dim)
  call diagSYM( 'V', dim, H, E, W_r)
  deallocate(H)
  
  print*, "E = "
  print "(*(F5.2,2X))", E(:)
  print *, ""

  call expSYM( dim, -C_UNIT*(pi/2), E, W_r, USwap )
  deallocate(E, W_r)
  !print*, "USwap = "
  !call printmat(dim, USwap, 'C')

  allocate(PH(dim),W(dim,dim))
  call diagUN( SELECT, dim, USwap, PH, W)

  print*, "arg(PH) = "
  print "(*(F5.2,2X))", atan2(dreal(PH(:)),dimag(PH(:)) )
  print *, ""

  USwap = matmul(USwap,USwap)
  !print*, "USwap^2 = "
  !call printmat(dim, USwap, 'C')
  USwap = abs(USwap)
  !print*, "abs(USwap^2) = "
  !call printmat(dim, USwap, 'C')
  if (all( abs(USwap) - identity(dim) <= 1.0e-6)) then
    print*, "USwap^2 = c * Id"
    print *, ""
  else
    print*, "USwap^2 != c * Id"
    print *, ""
  endif

  deallocate(PH,W,USwap)

  !deallocate(H)


  !----------------- Sz0  Subspace ------------------!

  allocate(H(dim_Sz0,dim_Sz0))
  !call buildSz0_HMBL(nspin, dim_Sz0, Jxy, Vz, hz, H)
  !print *, "H0 = "
  !call printmat(dim_Sz0, H, 'R')

  call buildSz0_HSwap(nspin, dim_Sz0, H)
  print *, "H0 = "
  !call printmat(dim_Sz0, H, 'R')
  call printnzmat(dim_Sz0, H, 'R')
  if (all(H == transpose(H))) then
    write(*,*) "H = H^T"
    print *, ""
  else
    write(*,*) "H != H^T"
    print *, ""
  endif


  allocate(E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  !H2 = H2 + identity(dim)
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
  deallocate(H)
  print *, "E = "
  write(*,"(*(F5.2,2X))") E(:)
  print *, ""

  call expSYM( dim_Sz0, -C_UNIT*pi/2.d0, E, W_r, USwap )
  deallocate(E, W_r)
  print*, "USwap = "
  call printnzmat(dim_Sz0, USwap, 'A')


  allocate(PH(dim_Sz0),W(dim_Sz0,dim_Sz0))
  call diagUN( SELECT, dim_Sz0, USwap, PH, W)

  print*, "arg(PH) = "
  write(*,"(*(F5.2,2X))") atan2(dreal(PH(:)),dimag(PH(:)) )
  print *, ""

  USwap = matmul(USwap,USwap)
  !print*, "USwap^2 = "
  !call printmat(dim_Sz0, USwap, 'C')
  USwap = abs(USwap)
  !print*, "abs(USwap^2) = "
  !call printmat(dim_Sz0, USwap, 'C')
  if (all( abs(USwap) - identity(dim_Sz0) <= 1.0e-6)) then
    write(*,*) "USwap^2 = c * Id"
    print *, ""
  else
    write(*,*) "USwap^2 != c * Id"
    print *, ""
  endif

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")




end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
