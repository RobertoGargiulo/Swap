program swap

  use exp_sparse
  use exponentiate
  use genmat
  use printing
  use MBL
  use omp_lib
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
  real(c_double) :: t_avg, t_sigma, r_dis_avg, r_dis_sigma
  
  real (c_double) :: norm
  real (c_double), dimension(:), allocatable :: avg, sigma, avg2, sigma2, r_avg, r_sigma
  complex (c_double_complex) :: alpha, beta

  integer (c_int), dimension(:), allocatable :: ROWS, COLS
  real (c_double), dimension(:), allocatable :: H_sparse, E
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:,:), allocatable :: U, W

  complex(c_double_complex), dimension(:), allocatable :: init_state, state_i, state_f, state

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
  

  write(filestring,91) "data/magnetizations/Sz0_DENSE_MBL_hz_Disorder_AVG_FLUCT_Imbalance_nspin", &
    & nspin, "_steps", steps, "_time_step", T0, &
    &  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, ".txt"
  open(newunit=unit_avg,file=filestring)

  !EIGENVALUES/EIGENVECTORS
  !write(filestring,92) "data/eigen/Sz0_DENSE_MBL_PH_nspin", nspin, "_steps", steps, &
  !&  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, ".txt"
  !open(newunit=unit_ph, file=filestring)
  
  !write(filestring,92) "data/eigen/Sz0_DENSE_MBL_W_nspin", nspin, "_steps", steps, &
  !&  "_iterations", n_iterations, "_J", J_coupling, "_V", V_coupling, "_hz", hz_coupling, ".txt"
  !open(newunit=unit_w, file=filestring)

  91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  !If V >= 10 (or hz) we can use
  ! hz_coupling --> int(hz_coupling), hz_coupling-int(hz_coupling)
  !   92  format(A,I0, A,I0, A,I0, A,F4.2, A,I2.2,F0.2, A,F4.2, A, I0, A)
  !   
!
 
  !------------------------------------------------

  !BUILD INITIAL STATE (of type staggered)
  allocate(init_state(dim_Sz0))
  call buildStaggState_Sz0(nspin, dim_Sz0, init_state)

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
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(U(dim_Sz0,dim_Sz0), H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))

  !Allocate initial and generic state
  allocate( state(dim_Sz0) )

  !Allocate observables and averages
  allocate( avg(steps), sigma(steps))
  allocate( r_avg(n_iterations), r_sigma(n_iterations))
  !allocate( avg2(steps), sigma2(steps))

  !Allocate for Eigenvalues/Eigenvectors
  !allocate(PH(dim_Sz0))
  !allocate(W(dim_Sz0,dim_Sz0))

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
  !$OMP & state, H, E, W_r, U, PH, W)
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
    !PRINT Eigenvalues/Eigenvectors to file
    !write(unit_ph,*) "iteration = ", iteration
    !call writevec(unit_ph,dim_Sz0,E,'R')
    !write(unit_w,*) "iteration = ", iteration
    !call writemat(unit_w,dim_Sz0,W_r,'R')
    call gap_ratio(dim_Sz0, E, r_avg(iteration), r_sigma(iteration))

    state = init_state
    norm = dot_product(state,state)
    j = 1
    avg(j) = avg(j) + imbalance_Sz0(nspin, dim_Sz0, state)
    sigma(j) = sigma(j) + imbalance_Sz0(nspin, dim_Sz0, state)**2
    !print *, imbalance(nspin, dim, state), j, norm
  
    do j = 2, steps
      state = matmul(U,state)
      norm = dot_product(state,state)
      state = state / sqrt(norm)
      avg(j) = avg(j) + imbalance_Sz0(nspin, dim_Sz0, state) 
      sigma(j) = sigma(j) + imbalance_Sz0(nspin, dim_Sz0, state)**2
      !print *, imbalance(nspin, dim, state), j, norm
    enddo
    !print *, ""
    !call take_time(count_rate, count2, count1, 'T', "Dense Evolution")
    !print *, "End of Dense Evolution"
    !print *, ""
 

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 


  avg = avg/n_iterations
  sigma = sqrt(sigma/n_iterations - avg**2)/sqrt(real(n_iterations))
  !avg2 = avg2/n_iterations
  !sigma2 = sqrt(sigma2/n_iterations - avg2**2)/sqrt(real(n_iterations))
  do j = 1, steps
    write(unit_avg,*) avg(j), sigma(j), j*T0
    !write(unit_avg,*) avg(j), sigma(j), avg2(j), sigma2(j), j*T0
    !print *, avg(j), sigma(j), avg2(j), sigma2(j), j*T0
  enddo
  
  call time_avg(n_iterations, 1, r_avg, r_sigma, r_dis_avg, r_dis_sigma)
  start = int(100/T0) !The average starts from the step for which 100 = start*T0
  call time_avg(steps, start, avg, sigma, t_avg, t_sigma)
  !call time_avg(steps, start, avg2, sigma2, t_avg2, t_sigma2)

  write(unit_avg,*) "Time Averages and Errors of Imbalance"
  write(unit_avg,*) t_avg, t_sigma!, t_avg2, t_sigma2
  write(unit_avg,*) "Average and Variance of Gap Ratio (over the spectrum and then disorder)"
  write(unit_avg,*) r_dis_avg, r_dis_sigma

  print *,"Time Averages and Errors"
  print *, t_avg, t_sigma!, t_avg2, t_sigma2
  print *, "Average and Variance of Gap Ratio (over the spectrum and then disorder)"
  print *, r_dis_avg, r_dis_sigma

  call take_time(count_rate, count_beginning, count1, 'T', "Program: ")

  deallocate(Jint, Vint, h_z)
  deallocate(avg, sigma)
  deallocate(H,E,W_r,U)
  !deallocate(USwap)
  !deallocate(PH, W)

  close(unit_avg)

  !call take_time(count_rate, count_beginning, count_end)

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end

subroutine take_time(count_rate, count_start, count_end, opt, filestring)
  implicit none
  integer, intent(in) :: count_rate, count_start
  integer, intent(out) :: count_end
  character :: opt*1
  character (*) :: filestring

  real :: time_s
  integer :: time_min

  call system_clock(count_end)

  time_s = real(count_end - count_start) / real(count_rate)
  time_min = int(time_s/60)
  if(opt == 'T') then
    print "(A,A,A,1X,I4,A,2X,F15.10,A)", "Elapsed Time for ", filestring, ": ", time_min, "min", time_s - 60*time_min, "sec"
  else if(opt == 'F') then
    print *, ""
  endif
  !print *, ""
end subroutine
