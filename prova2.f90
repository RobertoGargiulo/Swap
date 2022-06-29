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
  integer (c_int)     ::  i, j, k, p, nz_dim, krylov_dim
  integer (c_int)     ::  unit_mag, unit_ph, unit_w, unit_avg
  integer (c_int), dimension(:), allocatable  :: base_state, config

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  
  real (c_double) :: norm, time
  real (c_double), dimension(:), allocatable :: avg, sigma, avg2, sigma2
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
  
  write (*,*) "Transverse Interaction Constant -J * (XX + YY)"
  read (*,*) J_coupling
  print*,""

  write (*,*) "Longitudinal Interaction Constant -V * ZZ"
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

  !DATA FILES
  
  !write(filestring,92) "data/magnetizations/Clean_MBL_Imbalance_nspin", nspin, "_steps", steps, &
  !  &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, &
  !  & "_no_kick", kick, ".txt"
  !open(newunit=unit_mag,file=filestring)

  write(filestring,92) "data/magnetizations/Clean_MBL_OMP_vs_SERIAL_AVG_FLUCT_Imbalance_nspin", nspin, "_steps", steps, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, ".txt"
  open(newunit=unit_avg,file=filestring)

  !EIGENVALUES/EIGENVECTORS
!  write(filestring,92) "data/eigenvalues/Swap_PH_nspin", nspin, "_steps", steps, &
!  &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, "_kick", kick, ".txt"
!  !open(newunit=unit_ph, file=filestring)
!  
!  write(filestring,92) "data/eigenvalues/Swap_W_nspin", nspin, "_steps", steps, &
!  &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, "_kick", kick, ".txt"
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A,F4.2, A)
!
!  !open(newunit=unit_w, file=filestring)
 
  !------------------------------------------------

  !BUILD INITIAL STATE (of type staggered)
  allocate(init_state(dim))
  call buildStaggState(nspin, dim, alpha, beta, init_state)

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  !------------- NO DRIVING ---------- Uncomment following lines to use the driving
!  allocate(H(dim,dim), E(dim), W_r(dim,dim), USwap(dim,dim))
!  call buildHSwap(nspin, dim, H)
!  call diagSYM( 'V', dim, H, E, W_r)
!!  print *, "HSwap = "
!!  call printmat(dim, H, 'R')
!  deallocate(H)
!  call expSYM( dim, -C_UNIT*T1, E, W_r, USwap) 
!  deallocate(E, W_r)
  
  !print *, "USwap = "
  !call printmat(dim, USwap, 'R')

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_z(nspin))

  !Allocate Floquet and MBL Operators
  allocate(U(dim,dim), H(dim,dim), E(dim), W_r(dim,dim))
  nz_dim = (nspin+1)*dim/2
  allocate(H_sparse(nz_dim), ROWS(nz_dim), COLS(nz_dim))

  !Allocate initial and generic state
  allocate( state(dim), state_i(dim), state_f(dim))

  !Allocate observables and averages
  allocate( avg(steps), sigma(steps), avg2(steps), sigma2(steps) )

  !Allocate for Eigenvalues/Eigenvectors
  !allocate(PH(dim), W(dim,dim))

  state_i = init_state
  state_f = 0
  avg = 0
  sigma = 0
  avg2 = 0
  sigma2 = 0
  !$OMP PARALLEL
  call init_random_seed() 
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+:avg, sigma, avg2, sigma2) private(iteration, h_z, H, E, W_r, U, state, norm, j, time, ROWS, COLS, H_sparse)
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
  
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildHMBL( nspin, dim, Jint, Vint, h_z, H )
    !print *, "H_MBL = "
    !call printmat(dim, H,'R')

    call diagSYM( 'V', dim, H, E, W_r )
    call expSYM( dim, -C_UNIT*T0, E, W_r, U )
    !print *, "U_MBL = "
    !call printmat(dim, U_MBL,'C')
    call buildSPARSE_HMBL(nspin, dim, Jint, Vint, h_z, H_sparse, ROWS, COLS)

    !U = U_MBL !NO DRIVING
    !U = matmul(U_MBL, USwap)
  
    !print *, "UF = "
    !call printmat(dim, U,'C')
    !-------------------------------------------------

    !DIAGONALIZE FLOQUET OPERATOR
!    call diagUN( SELECT, dim, U, PH, W)
  
!    print *, "Eigenvalues of U_F "
!    print "(*(/f15.10spf15.10' i'))", PH(:)
!    print *,""
!  
!    print *, "Eigenvectors of U_F "
!    do i = 1,dim
!      write (*,"(*('|',f5.2spf5.2' i'))") W(i,:)
!    enddo
!    print *,""
    !------------------------------------------------
  
    !PRINT Eigenvalues/Eigenvectors to file
    !call writevec(unit_ph,dim,PH,'C')
    !call writemat(unit_w,dim,W,'C')
  
    !-----------------------------------------------
  
    !EVOLUTION OF INITIAL STATE and COMPUTATION OF MAGNETIZATION 
 
    state = init_state
    norm = dot_product(state,state)
    j = 1
    time = j*T0
    avg(j) = avg(j) + imbalance(nspin, dim, state)
    sigma(j) = sigma(j) + imbalance(nspin, dim, state)**2
    !print *, imbalance(nspin, dim, state), j, norm
  
    do j = 2, steps
      state = matmul(U,state)
      norm = dot_product(state,state)
      state = state / sqrt(norm)
      time = j*T0
      avg(j) = avg(j) + imbalance(nspin, dim, state) 
      sigma(j) = sigma(j) + imbalance(nspin, dim, state)**2
      !print *, imbalance(nspin, dim, state), j, norm
    enddo
    
    state_i = init_state
    norm = dot_product(state_i, state_i)
    j = 1
    avg2(j) = avg2(j) + imbalance(nspin, dim, state_i)
    sigma2(j) = sigma2(j) + imbalance(nspin, dim, state_i)**2
    !print *, "Imbalance", "Magnetization", "Time", "Norm"
    !print *, imbalance(nspin, dim, state_i), mag_z(nspin, dim, state_i), time, norm
    do j = 2, steps
      call evolve(dim, nz_dim, krylov_dim, ROWS, COLS, -C_UNIT*dcmplx(H_sparse), state_i, T0, state_f)
      state_i = state_f
      norm = dot_product(state_i, state_i)
      state_i = state_i / sqrt( dot_product(state_i, state_i) )
      avg2(j) = avg2(j) + imbalance(nspin, dim, state_i)
      sigma2(j) = sigma2(j) + imbalance(nspin, dim, state_i)**2
      !print *, imbalance(nspin, dim, state_i), mag_z(nspin, dim, state_i), time, norm
    enddo
    !print *, ""
 

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  avg = avg/n_iterations
  sigma = sqrt(sigma/n_iterations - avg**2)/sqrt(real(n_iterations))
  avg2 = avg2/n_iterations
  sigma2 = sqrt(sigma2/n_iterations - avg2**2)/sqrt(real(n_iterations))
  do j = 1, steps
    write(unit_avg,*) avg(j), sigma(j), avg2(j), sigma2(j), j*T0
  enddo

  deallocate(Jint, Vint, h_z)
  deallocate(E, W_r, H)
  deallocate(U)
  deallocate(avg, sigma, avg2, sigma2)
  !deallocate(USwap)
  !deallocate(PH, W)

  close(unit_avg)

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
