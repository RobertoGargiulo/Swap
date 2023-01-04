program prova

  use genmat
  use printing

  integer (c_int)     ::  nspin, dim
  integer (c_int), allocatable :: config(:)
  integer (c_int) :: i, j, k
  real (c_double) :: imp(3,3)
  real (c_double), allocatable :: h(:), Jxy(:), H1(:,:), H2(:,:), v(:)

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

  do i = 1, dim
    print *, sum(H2(i,:))
  enddo
  

end program
