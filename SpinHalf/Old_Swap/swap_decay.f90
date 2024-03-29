program swap_decay

  use exponentiate
  use genmat
  use printing
  use MBL
  use omp_lib
  use sorts
  use iso_c_binding
  !use general
  implicit none

  !complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  !complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, n_iterations
  integer (c_long)    ::  steps, j
  integer (c_int)     ::  dim_Sz0, i, l
  integer (c_int)     ::  unit_decay
  integer (c_int)     ::  start

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  
  real (c_double) :: norm

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap, UF

  complex(c_double_complex), dimension(:), allocatable :: init_state, state, state2

  integer(c_int) :: count_beginning, count_end, count_rate
  character(len=200) :: filestring, filestring_Neel, filestring_Nayak

  real (c_double), allocatable, dimension(:) :: t_decay
  real (c_double)   :: imb_p, t_decay_avg, sigma_t_decay
  integer (c_int)   :: idecay, idecay2, n_decays
  integer (c_long)  :: n_periods

  real (c_double)   :: IMB, LI
  integer (c_int)   :: idx_state

  integer (c_int), allocatable :: idxSz0(:), config(:)
  complex (c_double_complex) :: alpha, beta

  !Parametri Modello: J, V, h_z, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine, Steps di Evoluzione, Stato Iniziale

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin
  dim_Sz0 = binom(nspin,nspin/2)

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  write (*,*) "Number of Periods"
  read (*,*) n_periods
  print*,""

  write (*,*) "Period T0"
  read (*,*) T0
  print*,""
  
  write (*,*) "Transverse Interaction Constant -J * (XX + YY)"
  read (*,*) J_coupling
  print*,""

  write (*,*) "Longitudinal Interaction Constant -V * ZZ"
  read (*,*) V_coupling
  print*,""

  write (*,*) "Longitudinal Field h_z * Z"
  read (*,*) hz_coupling
  print*,""

  write (*,*) "Perturbation on Kick T1 = pi/4 + kick"
  read (*,*) kick
  print*,""
  !---Read below for distributions of J, V, hz
  
  T1 = pi/4 + kick


  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring_Neel,93) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_Neel_Decay_Times_Imbalance_nspin", &
    & nspin, "_steps", steps, "_nperiods", n_periods, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_Nayak,93) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_Nayak_Decay_Times_Imbalance_nspin", &
    & nspin, "_steps", steps, "_nperiods", n_periods, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  open(newunit=unit_decay,file=filestring_Neel)
 
  93  format(A,I0, A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)

  !------------------------------------------------

  !BUILD INITIAL STATE (of type staggered)
  allocate(init_state(dim_Sz0))

  alpha = cos(pi/8)
  beta = sin(pi/8)
  !call buildNayakState_Sz0(nspin, dim_Sz0, alpha, beta, init_state)
  call buildNeelState_Sz0(nspin, dim_Sz0, init_state)
  
  !call printvec(dim_Sz0, init_state, 'A')

  allocate(idxSz0(dim_Sz0), config(nspin))
  call zero_mag_states(nspin, dim_Sz0, idxSz0)
  do i = 1, dim_Sz0
    if (abs(init_state(i))**2 > 1.0e-6) then
      l = idxSz0(i)
      call decode(l, nspin, config)
      print "( 4X,F8.4, 4X,I6, 4X,*(I0) )", abs(init_state(i))**2, l, config(:)
    endif
  enddo

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  call buildSz0_HSwap(nspin, dim_Sz0, H)
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
  deallocate(H)
  call expSYM( dim_Sz0, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)
  
  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_z(nspin))

  !Allocate Floquet and MBL Operators
  allocate(U(dim_Sz0,dim_Sz0), H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))

  !Allocate initial and generic state
  allocate( state(dim_Sz0) )

  !Allocate observables and averages
  allocate( t_decay(n_iterations) )

  !Allocate for Eigenvalues/Eigenvectors

  n_decays = 0
  t_decay = steps
  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+: n_decays) private(iteration, h_z, Vint, norm, j, &
  !$OMP & state, H, E, W_r, U, PH, W, &
  !$OMP & imb_p, idecay)
  !, UF)
  do iteration = 1, n_iterations
    
    if (n_iterations < 10) then
     print *, "iteration = ", iteration
    else
      if (mod(iteration,n_iterations/10)==0) then 
        print *, "iteration = ", iteration
      endif
    endif

    !-------------------------------------------------
    !PARAMETERS
  
    !call random_number(Jint)
    !Jint = 2*J_coupling*(Jint - 0.5) !Jint in [-J,J]
    Jint = -J_coupling

    Vint = -V_coupling
    call random_number(Vint)
    Vint = -V_coupling + V_coupling*(Vint - 0.5) !Jint in [-V,V]
  
    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]
  
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    !print *, "Building Dense Operator"
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    call expSYM( dim_Sz0, -C_UNIT*n_periods*T0, E, W_r, U )
    U = matmul(USwap,U) !U -> UF

    state = init_state
    norm = real(dot_product(state,state))
    j = 1

    idecay = 0
    do j = 2, steps
      imb_p = imbalance_Sz0(nspin, dim_Sz0, state)
      !print *, imb_p
      state = matmul(U,state)
      norm = real(dot_product(state,state))
      state = state / sqrt(norm)
      if (idecay == 0) then
        if(imb_p*imbalance_Sz0(nspin, dim_Sz0, state) > 0) then
          t_decay(iteration) = j*n_periods*T0
          idecay = 1
          n_decays = n_decays + 1
          exit
        endif
      endif

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 


  t_decay_avg = sum(t_decay) / n_iterations
  sigma_t_decay = sqrt( (sum(t_decay**2)/n_iterations - t_decay_avg**2) / n_iterations )
  do iteration = 1, n_iterations
    write(unit_decay,*) t_decay(iteration), iteration
  enddo
  write(unit_decay,*) "Decay Time ( t*, sigma(t*), n_decays )"
  write(unit_decay,*) t_decay_avg, sigma_t_decay, n_decays
  close(unit_decay)


  print *, "Decay Time ( t*, sigma(t*), n_decays)"
  print*, t_decay_avg, sigma_t_decay, n_decays
    

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

  deallocate(Jint, Vint, h_z)
  deallocate(H,E,W_r,U)
  deallocate(USwap)


  !call take_time(count_rate, count_beginning, count_end)

end program swap_decay
