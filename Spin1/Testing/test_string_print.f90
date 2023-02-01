program prova2

  use iso_c_binding, dp => c_double
  implicit none

  character(len=300) :: string, s2

  integer :: i, j, nspin, steps

  real(dp), allocatable :: sigmaz_avg(:,:)
  real(dp) :: T0 = 1.0

  nspin = 4
  steps = 10

  allocate(sigmaz_avg(steps,nspin))

  call random_number(sigmaz_avg)

  string = ""
  do i = 1, nspin
    write(s2, "(I0)") i
    string(26*(i-1)+2:26*i+2) = trim("sigma_z^")//trim(s2)
    !print *, string
  enddo
  print "(3X,A4,20X,A)", "j*T0", trim(string)
  !print "(*(A26))", "sigmaz_1", "sigmaz_2"
  do j = 1, steps
    print *, j*T0, sigmaz_avg(j,:)
  enddo

end program

    
