program swap

  use exponentiate
  use genmat
  use printing
  use omp_lib
  use iso_c_binding
  !use general
  implicit none

  complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, steps, n_iterations
  integer (c_int)     ::  i, j, k, p, tid, nthreads
  integer (c_int)     ::  unit_mag, unit_ph, unit_w, unit_avg
  integer (c_short), dimension(:), allocatable  :: base_state, config

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  
  real (c_double) :: norm
  real (c_double), dimension(:), allocatable :: avg, sigma
  real (c_double), dimension(:,:), allocatable :: imb
  complex (c_double_complex) :: alpha, beta

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r, H_MBL

  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:), allocatable :: state, init_state
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap, U_MBL

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, time_min
  real (c_double) :: time_s
  character(len=200) :: filestring
  character(len=8) :: time_string


  !Parametri Modello: J, V, h_z, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine, Steps di Evoluzione, Stato Iniziale

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""



  !Standard Values
  T0 = 1
  J_coupling = 1
  V_coupling = J_coupling
  hz_coupling = J_coupling
  kick = 0.1
  T1 = pi/4 + kick

  write (*,*) "Time Step/Period T0"
  read (*,*) T0
  print*,""
  
  write (*,*) "Longitudinal Interaction Constant J * ZZ"
  read (*,*) J_coupling
  print*,""

  write (*,*) "Transverse Interaction Constant V * (XX + YY)"
  read (*,*) V_coupling
  print*,""

  write (*,*) "Longitudinal Field h_z * Z"
  read (*,*) hz_coupling
  print*,""
  !---Read below for distributions of J, V, hz
  

  T1 = pi/4 + kick

    !Coefficienti dello stato iniziale |psi> = (alpha|up>+beta|down>)^L
  alpha = 1
  beta = 0

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------
  write(filestring,92) "data/magnetizations/Clean_MBL_Imbalance_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, ".txt"
  !open(newunit=unit_mag,file=filestring)

  write(filestring,92) "data/magnetizations/Clean_MBL_OMP_AVG_FLUCT_Imbalance_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, ".txt"
  open(newunit=unit_avg,file=filestring)
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A,F4.2, A)

  !BUILD INITIAL STATE (of type staggered)
  allocate(init_state(dim))
  call buildStaggState(nspin, dim, alpha, beta, init_state)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_z(nspin))

  !Allocate Floquet and MBL Operators
  allocate(U(dim,dim), H(dim,dim), E(dim), W_r(dim,dim))

  !Allocate initial and generic state
  allocate(state(dim))

  !Allocate observables and averages
  allocate( avg(steps), sigma(steps), imb(n_iterations,steps) )

  !Allocate for Eigenvalues/Eigenvectors
  !allocate(PH(dim), W(dim,dim))

  imb = 0
  avg = 0
  sigma = 0
  !$OMP PARALLEL
  call init_random_seed() 
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+:avg, sigma) private(iteration, h_z, H, E, W_r, U, state, norm, j)
  do iteration = 1, n_iterations
    
    !if (mod(iteration,10)==0) then 
    !  print *, "iteration = ", iteration
    !endif

    !-------------------------------------------------
    !PARAMETERS
  
    Jint = -J_coupling

    Vint = -V_coupling
  
    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]
    !---------------------------------------------------
  
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildHMBL( nspin, dim, Jint, Vint, h_z, H )

    call diagSYM( 'V', dim, H, E, W_r )
    call expSYM( dim, -C_UNIT*T0, E, W_r, U )
  
    !EVOLUTION OF INITIAL STATE and COMPUTATION OF MAGNETIZATION 
 
    state = init_state
    norm = dot_product(state,state)
    j = 1 
    imb(iteration,j) = imbalance(nspin, dim, state) 
    avg(j) = avg(j) + imbalance(nspin, dim, state)
    sigma(j) = sigma(j) + imbalance(nspin, dim, state)**2
    print *, imbalance(nspin, dim, state), j*T0, norm

    do j = 2, steps
      state = matmul(U,state)
      norm = dot_product(state,state)
      state = state / sqrt(norm)
      imb(iteration,j) = imbalance(nspin, dim, state) 
      avg(j) = avg(j) + imbalance(nspin, dim, state) 
      sigma(j) = sigma(j) + imbalance(nspin, dim, state)**2
      print *, imbalance(nspin, dim, state), j*T0, norm
    enddo
    !print *, ""
 

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  avg = avg/n_iterations
  sigma = sqrt(sigma/n_iterations - avg**2)/sqrt(real(n_iterations))
  do j = 1, steps
    write(unit_avg,*) avg(j), sigma(j), j*T0
  enddo
  !do iteration = 1, n_iterations
  !  j = 1
  !  write(unit_mag,*) "iteration = ", iteration
  !  write(unit_mag,*) imb(iteration,j), j*T0
  !  do j = 2, steps
  !    write(unit_mag,*) imb(iteration,j), j*T0
  !  enddo
  !enddo



  deallocate(Jint, Vint, h_z)
  deallocate(E, W_r, H)
  deallocate(U)
  !deallocate(USwap)
  !deallocate(PH, W)

  close(unit_mag)

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
