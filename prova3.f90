program swap

  use exp_sparse
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
  integer (c_int)     ::  i, j, k, p, tid, nthreads, nz_dim, krylov_dim
  integer (c_int)     ::  unit_mag, unit_ph, unit_w
  integer (c_short), dimension(:), allocatable  :: base_state, config

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  
  real (c_double) :: mag, norm, time
  complex (c_double_complex) :: alpha, beta

  integer (c_int), dimension(:), allocatable :: ROWS, COLS
  real (c_double), dimension(:), allocatable :: E, H_sparse
  real (c_double), dimension(:,:), allocatable :: H, W_r, H_MBL

  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:), allocatable :: state, init_state, state_i, state_f
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap, U_MBL

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, day, month, year, date(8), time_min
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

  write (*,*) "Maximum Krylov Dimension"
  read (*,*) krylov_dim
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
  
  write (*,*) "Longitudinal Interaction Constant -J * (XX + YY)"
  read (*,*) J_coupling
  print*,""

  write (*,*) "Transverse Interaction Constant -V * ZZ"
  read (*,*) V_coupling
  print*,""

  write (*,*) "Longitudinal Field h_z * Z"
  read (*,*) hz_coupling
  print*,""
  !---Read below for distributions of J, V, hx, hz
  

  T1 = pi/4 + kick

  !Coefficienti dello stato iniziale |psi> = (alpha|up>+beta|down>)^L
  alpha = 1
  beta = 0

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !BUILD INITIAL STATE (of type staggered)
  allocate(init_state(dim))
  call buildStaggState(nspin, dim, alpha, beta, init_state)

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  !------------- NO DRIVING ---------- Uncomment following lines to use the driving
!  allocate(H(dim,dim), E(dim), W_r(dim,dim), USwap(dim,dim))
!  call buildHSwap(nspin, dim, H)
!  call diagSYM( 'V', dim, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
!  deallocate(H)
!  call expSYM( dim, -C_UNIT*T1, E, W_r, USwap) 
!  deallocate(E, W_r)
  
  !print *, "USwap = "
  !call printmat(dim, USwap, 'R')

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_z(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(U(dim,dim), H(dim,dim), E(dim), W_r(dim,dim), H_MBL(dim,dim))
  nz_dim = (nspin+1)*dim/2
  allocate(H_sparse(nz_dim), ROWS(nz_dim), COLS(nz_dim))

  !Allocate initial and generic state
  !allocate(state(dim))
  allocate(state_i(dim), state_f(dim))

  state_i = init_state
  state_f = 0

  !Allocate for Eigenvalues/Eigenvectors
  !allocate(PH(dim), W(dim,dim))

  !$OMP PARALLEL 
  call init_random_seed()
  print *, "Max size of thread team: ", omp_get_max_threads()
  print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(h_z, H_sparse, ROWS, COLS, state_i, state_f, norm, j, time )
  do iteration = 1, n_iterations
    
    !if (mod(iteration,10)==0) then 
    !  print *, "iteration = ", iteration
    !endif

    !print *, "Size of Thread team: ", omp_get_num_threads()
    !print *, "Thread ID: ", omp_get_thread_num()
    !print *, "Number of processors: ", omp_get_num_procs()

    !-------------------------------------------------
    !PARAMETERS
  
    !call random_number(Jint)
    Jint = -J_coupling

    !call random_number(Vint)
    Vint = -V_coupling
  
    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]
    !---------------------------------------------------
  
    !BUILD FLOQUET (EVOLUTION) OPERATOR
      
    call buildSPARSE_HMBL(nspin, dim, Jint, Vint, h_z, H_sparse, ROWS, COLS)



    !call buildHMBL( nspin, dim, Jint, Vint, h_z, H_MBL )

    !call printmat(dim, H_MBL, 'R')

    !call diagSYM('V', dim, H_MBL, E, W_r)
    !call expSYM(dim, -C_UNIT*T0, E, W_r, U)

  
    !EVOLUTION OF INITIAL STATE and COMPUTATION OF MAGNETIZATION 
 
    !state = init_state
    !norm = dot_product(state,state)
    !j = 1 
    !print *, imbalance(nspin, dim, state), mag_z(nspin, dim, state),  j*T0, norm
  
    !do j = 2, steps
    !  state = matmul(U,state)
    !  norm = dot_product(state,state)
    !  state = state / sqrt(norm)
    !  print *, imbalance(nspin, dim, state), mag_z(nspin, dim, state), j*T0, norm
    !enddo
    !print *, ""
 
    state_i = init_state
    norm = dot_product(state_i, state_i)
    j = 1
    time = j*T0
    !print *, "Imbalance", "Magnetization", "Time", "Norm"
    !print *, imbalance(nspin, dim, state_i), mag_z(nspin, dim, state_i), time, norm
    do j = 2, steps
      call evolve(dim, nz_dim, krylov_dim, ROWS, COLS, -C_UNIT*dcmplx(H_sparse), state_i, T0, state_f)
      state_i = state_f
      norm = dot_product(state_i, state_i)
      time = j*T0
      state_i = state_i / sqrt( dot_product(state_i, state_i) )
      !print *, imbalance(nspin, dim, state_i), mag_z(nspin, dim, state_i), time, norm
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !deallocate(Jint, Vint, h_z)
  !deallocate(E, W_r, H)
  !deallocate(U)
  !deallocate(USwap)
  !deallocate(PH, W)

  !close(unit_mag)
  !close(unit_ph)
  !close(unit_w)

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
