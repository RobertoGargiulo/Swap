program swap

  use genmat
  use exponentiate
  use omp_lib
  use iso_c_binding
  !use general
  implicit none

  complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, n_iterations
  integer (c_int)     ::  i, j, k, p, tid, nthreads

  
  real (c_double) :: rand, mag, two_norm

  real (c_double), dimension(:,:), allocatable :: H, W_r
  real (c_double), dimension(:), allocatable :: E
  complex (c_double_complex), dimension(:,:), allocatable :: U


  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, day, month, year, date(8), time_min
  real (c_double) :: time_s
  character(len=200) :: filestring
  character(len=8) :: time_string


  !Parametri Modello: J, h_x, h_z, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine, Steps di Evoluzione, Stato Iniziale

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin


  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""


  call system_clock(count_beginning, count_rate)

  !Allocate Floquet and MBL Operators
  allocate(H(dim,dim), E(dim), W_r(dim, dim), U(dim,dim))


  !iteration = 1
  !$OMP PARALLEL DO PRIVATE(H, E, W_r, two_norm, nspin, dim)
  do iteration = 1, n_iterations
    

    call random_number(H)
    H = (H + transpose(H))/2
    call diagSYM( 'V', dim, H, E, W_r )
    call expSYM(dim, -C_UNIT, E, W_r, U )
    call random_number(rand)
    two_norm = rand * norm2(H)

 
    !!$OMP ordered
    print *,  "Norm of H = ", two_norm
    !print *, "Max size of thread team: ", omp_get_max_threads()
    !print *, "Size of Thread team: ", omp_get_num_threads()
    !print *, "Thread ID: ", omp_get_thread_num()
    print *, "In parallel? ", omp_in_parallel()
    print *, ""
    !!$omp end ordered

  enddo
  !$OMP END PARALLEL DO

  deallocate(H)
  
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
