program entropy

  use MBL
  use genmat
  use printing
  use iso_c_binding
  implicit none

  complex (c_double_complex), allocatable :: psi(:), rho(:,:)
  real (c_double), allocatable :: vec1(:), vec2(:)
  real (c_double) :: EE, MI
  integer (c_int) :: nspin, dim, i, j, nspin_A, dim_A
  integer (c_int) :: count_beginning, count_rate, count_end
  real (c_double) :: num


do nspin = 2, 10, 2

  
  print *, "" 
  print *, "nspin = ", nspin
  dim = 2**nspin
  !read (*,*) dim

  call system_clock(count_beginning, count_rate)

  allocate(rho(dim,dim), psi(dim), vec1(dim), vec2(dim))

  psi = 0
  !call random_number(vec1)
  !call random_number(vec2)
  !psi = vec1 + (0.d0,1.d0)*vec2
  
  !j = 1
  call random_number(num)
  i = CEILING(dim*num)
  psi(i) = 1
  call random_number(num)
  j = CEILING(dim*num)
  do while(j == i) 
    call random_number(num)
    j = CEILING(dim*num)
  enddo
  psi(j) = 1
  psi = psi/sqrt(dot_product(psi,psi))
  do i = 1, dim
    if(abs(psi(i)) > 1.0e-10) print *, psi(i), i
  enddo
  !print *, psi(i), psi(j)

  !rho = 0
  !do i = 1, dim
  !  do j = 1, dim
  !    rho(i,j) = psi(i) * dconjg(psi(j))
  !  enddo
  !enddo
  !call printmat(dim, rho, 'C')
  !EE = 0
  !print *, EE, log(real(dim))

  !print *, "Entanglement of random pure state: " 
  !EE = entanglement(dim, rho)
  !print *, EE, log(real(dim))

  print *, "     nspinA                  S_A                       S_B                      S_AB                        MI"
  do nspin_A = 1, nspin/2
    dim_A = 2**nspin_A
    MI = mutual_information(nspin, nspin_A, dim, dim_A, psi)
  enddo
  call take_time(count_rate, count_beginning, count_end, 'T', "Program")
  print *, ""

  deallocate(rho, psi, vec1, vec2)
enddo
  
 
end program entropy
