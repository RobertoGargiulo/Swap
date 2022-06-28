program swap

  use genmat
  use exponentiate
  use printing
  use omp_lib
  use iso_c_binding
  implicit none


  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  iteration, n_iterations, steps
  integer (c_int)     ::  i, j, k, p, nthreads

  
  real (c_double) :: rand, two_norm
  real (c_double), dimension(:), allocatable :: mag, avg, sigma


  integer(c_int) :: count_beginning, count_end, count_rate, day, month, year, date(8), time_min
  real (c_double) :: time_s
  character(len=200) :: filestring
  character(len=8) :: time_string

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  call system_clock(count_beginning, count_rate)

  !Allocate Floquet and MBL Operators
  allocate(mag(steps), avg(steps), sigma(steps))


  avg = 0
  sigma = 0
  !iteration = 1
  !$OMP PARALLEL
  call init_random_seed()
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "In parallel? ", omp_in_parallel()
  !$OMP DO REDUCTION(+: avg, sigma) PRIVATE(mag)
  do iteration = 1, n_iterations
    
    do j = 1, steps
      !call random_number(rand)
      mag(j) = iteration
      !p_avg(j) = mag(j)
      avg(j) = avg(j) + mag(j)
      sigma(j) = sigma(j) + mag(j)**2
      print "(3X,I0,2X,I0,2X,F5.2,2X,F5.2,2X,F5.2,2X,I0)", iteration, j, mag(j), avg(j), sigma(j), omp_get_thread_num()
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !avg = avg/n_iterations
  !sigma = sqrt(sigma/n_iterations - avg**2)/sqrt(real(n_iterations))
  print *, ""
  
  do j = 1, steps
    print *, avg(j), sigma(j)
  enddo

  
  call system_clock(count_end)

  time_s = real(count_end - count_beginning) / real(count_rate)
  time_min = int(time_s/60)

  print "(A,1X,I4,A,2X,F15.10,A)", "Elapsed Time: ", time_min, "min", time_s - 60*time_min, "sec"
  print *, ""

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
