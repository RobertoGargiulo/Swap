program prova

  use genmat
  use printing
  use exponentiate

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim
  integer (c_int), allocatable :: config(:)
  integer (c_int) :: i, j, k
  real (c_double) :: imp(3,3), T1
  real (c_double), allocatable :: h(:), Jxy(:), H1(:,:), H2(:,:), v(:)
  real (c_double), allocatable :: E(:), W_r(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:)

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

  dim = 3**nspin

  allocate(config(nspin))
  allocate(h(nspin), H1(dim,dim))
  allocate(Jxy(nspin-1), H2(dim,dim), v(dim))
  
!  do i = 0, dim-1
!    call decode(i, nspin, config)
!    print*, config(:), i
!  enddo
!  
!
!  imp = reshape((/ (( (2-i)*KDelta(i,j), j=1,3),i=1,3) /), shape(imp))
!
!  do i = 1, 3
!    print*, imp(i,:)
!  enddo

  h(:) = 1
  Jxy(:) = 1

  call buildOneBody(nspin, dim, h, H1)

  call buildTwoBody(nspin, dim, Jxy, H2)

  !do i = 1, dim
  !  print *, sum(H2(i,:))
  !enddo
  
  call buildTwoBodyZZ(nspin, dim, Jxy, H2)


  call buildHSwap(nspin, dim, H2)

  allocate(E(dim), W_r(dim,dim), USwap(dim,dim))
  H2 = H2 - identity(dim)
  call diagSYM( 'V', dim, H2, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  H2 = matmul(H2,H2)
  !print *, "HSwap^2 = "
  !call printmat(dim,H2,'R')
  

  if(all(H2 == identity(dim))) then
    write(*,*) "HSwap^2 = Id" 
  else
    write(*,*) "HSwap^2 != Id"
  endif

  deallocate(H2)
  
  print*, "E = "
  print "(*(F5.2,X))", E(:)

  T1 = pi/2
  call expSYM( dim, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)

  !print*, "USwap = "
  !call printmat(dim, USwap, 'C')

end program
