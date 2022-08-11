program swap

  use exp_sparse
  use exponentiate
  use genmat
  use printing
  use MBL
  use omp_lib
  use sorts
  use iso_c_binding
  !use general
  implicit none

  complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, steps, n_iterations, dim_eff
  integer (c_int)     ::  i, j, k, p, nz_dim, dim_Sz0, nz_Sz0_dim
  integer (c_int)     ::  unit_mag, unit_ph, unit_w, unit_avg
  integer (c_int)     ::  start

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  real(c_double) :: t_avg, t_sigma, r_dis_avg, r_dis_sigma, t_avg2, t_sigma2
  
  real (c_double) :: norm
  real (c_double), dimension(:), allocatable :: avg, sigma, avg2, sigma2, r_avg, r_sigma
  complex (c_double_complex) :: alpha, beta

  integer (c_int), dimension(:), allocatable :: ROWS, COLS
  real (c_double), dimension(:), allocatable :: H_sparse, E
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap, UF

  complex(c_double_complex), dimension(:), allocatable :: init_state, state_i, state_f, state, state2

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, time_min, count1, count2
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
  dim_Sz0 = binom(nspin,nspin/2)
  nz_Sz0_dim = non_zero_HMBL_Sz0(nspin)

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
  kick = 0
  T1 = pi/4 + kick

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

    !Coefficienti dello stato iniziale |psi> = (alpha|up>+beta|down>)^L
  alpha = 1
  beta = 0

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring,93) "data/magnetizations/Sz0_DENSE_SWAP_hz_Disorder_AVG_FLUCT_Imbalance_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"
  open(newunit=unit_avg,file=filestring)

  91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  93  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)
!
 
  !------------------------------------------------

  !BUILD INITIAL STATE (of type staggered)
  allocate(init_state(dim_Sz0))
  call buildStaggState_Sz0(nspin, dim_Sz0, init_state)

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  call buildSz0_HSwap(nspin, dim_Sz0, H)
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim_Sz0, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)
  
  !print *, "USwap = "
  !call printmat(dim_Sz0, USwap, 'C')
  !print *, "USwap*USwap^dagger = "
  !call printmat(dim_Sz0, matmul(USwap,USwap), 'C')

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_z(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(U(dim_Sz0,dim_Sz0), H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))

  !Allocate initial and generic state
  allocate( state(dim_Sz0) )

  !Allocate observables and averages
  allocate( avg(steps), sigma(steps))
  allocate( r_avg(n_iterations), r_sigma(n_iterations))
  allocate( avg2(steps), sigma2(steps))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))

  !state_i = init_state
  !state_f = 0
  avg = 0
  sigma = 0
  !avg2 = 0
  !sigma2 = 0
  r_avg = 0
  r_sigma = 0
  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+:avg, sigma) private(iteration, h_z, norm, j, &
  !$OMP & state, H, E, W_r, U, PH, W, state2, UF)
  do iteration = 1, n_iterations
    
    if (mod(iteration,10)==0) then 
      print *, "iteration = ", iteration
    endif

    !-------------------------------------------------
    !PARAMETERS
  
    !call random_number(Jint)
    !Jint = 2*J_coupling*(Jint - 0.5) !Jint in [-J,J]
    Jint = -J_coupling

    !call random_number(Vint)
    !Vint = 2*V_coupling*(Vint - 0.5) !Jint in [-V,V]
    Vint = -V_coupling
  
    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]
  
!    write (*,*) "Jint = ", Jint(:)
!    write (*,*) "Vint = ", Vint(:)
!    write (*,*) "h_z = ", h_z(:)
!    print *, ""
  
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    !print *, "Building Dense Operator"
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    UF = matmul(USwap,U)
    !U = USwap
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    E = real(C_UNIT*log(PH))
    call dpquicksort(E)
    call gap_ratio(dim_Sz0, E, r_avg(iteration), r_sigma(iteration))

    state = init_state
    norm = dot_product(state,state)
    j = 1
    avg(j) = avg(j) + imbalance_Sz0(nspin, dim_Sz0, state)
    sigma(j) = sigma(j) + imbalance_Sz0(nspin, dim_Sz0, state)**2

    state2 = init_state
    norm = dot_product(state2,state2)
    j = 1
    avg2(j) = avg2(j) + imbalance_Sz0(nspin, dim_Sz0, state2)
    sigma2(j) = sigma2(j) + imbalance_Sz0(nspin, dim_Sz0, state2)**2
    !if(iteration==1) then
      !print *, imbalance_Sz0(nspin, dim_Sz0, state), imbalance_Sz0(nspin, dim_Sz0, state2), j, norm
    !endif

    do j = 2, steps
      state = matmul(UF,state)
      norm = dot_product(state,state)
      state = state / sqrt(norm)
      avg(j) = avg(j) + imbalance_Sz0(nspin, dim_Sz0, state) 
      sigma(j) = sigma(j) + imbalance_Sz0(nspin, dim_Sz0, state)**2

      state2 = matmul(U,state2)
      norm = dot_product(state2,state2)
      state2 = state2 / sqrt(norm)
      avg2(j) = avg2(j) + imbalance_Sz0(nspin, dim_Sz0, state2) 
      sigma2(j) = sigma2(j) + imbalance_Sz0(nspin, dim_Sz0, state2)**2
      !if (iteration==1) then
      !  print *, imbalance_Sz0(nspin, dim_Sz0, state), j, norm
      !endif
      !print *, imbalance_Sz0(nspin, dim_Sz0, state), imbalance_Sz0(nspin, dim_Sz0, state2), j, norm
    enddo
    !print *, ""
    !call take_time(count_rate, count2, count1, 'T', "Dense Evolution")
    !print *, "End of Dense Evolution"
    !print *, ""
    !print *, imbalance_Sz0(nspin, dim_Sz0, state), j, norm

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 


  avg = avg/n_iterations
  sigma = sqrt(sigma/n_iterations - avg**2)/sqrt(real(n_iterations))
  avg2 = avg2/n_iterations
  sigma2 = sqrt(sigma2/n_iterations - avg2**2)/sqrt(real(n_iterations))
  do j = 1, steps
    write(unit_avg,*) j*T0, avg(j), sigma(j), avg2(j), sigma2(j)
  enddo

  call time_avg('F', n_iterations, 1, r_avg, r_sigma, r_dis_avg, r_dis_sigma)
  start = 1 !int(100/T0) !The average starts from the step for which 100 = start*T0
  call time_avg('T', steps, start, avg, sigma, t_avg, t_sigma)
  call time_avg('F', steps, start, avg2, sigma2, t_avg2, t_sigma2)

  write(unit_avg,*) "Time Averages and Errors of Imbalance"
  write(unit_avg,*) t_avg, t_sigma, t_avg2, t_sigma2
  write(unit_avg,*) "Average and Variance of Gap Ratio (over the spectrum and then disorder)"
  write(unit_avg,*) r_dis_avg, r_dis_sigma

  print *,"Time Averages and Errors"
  print *, t_avg, t_sigma, t_avg2, t_sigma2
  print *, "Average and Variance of Gap Ratio (over the spectrum and then disorder)"
  print *, r_dis_avg, r_dis_sigma

  call take_time(count_rate, count_beginning, count1, 'T', "Program: ")

  deallocate(Jint, Vint, h_z)
  deallocate(avg, sigma)
  deallocate(H,E,W_r,U)
  !deallocate(USwap)
  !deallocate(PH, W)

  !close(unit_avg)

  !call take_time(count_rate, count_beginning, count_end)

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
