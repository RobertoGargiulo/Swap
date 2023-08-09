subroutine find_degeneracies( n, energies, unq, deg )

  use iso_c_binding
  implicit none

  integer (c_int), intent(in) :: n
  real (c_double), intent(in) :: energies(n) !energies has to be sorted already
  real (c_double), intent(out)  :: unq(n)
  integer (c_int), intent(out)  :: deg(n)

  integer :: i, j, tot_deg
  real (c_double) :: tol = 1.0e-6


  unq = 0
  deg = 1

  j = 1
  unq(j) = energies(j)
  do i = 2, n
    if ( abs(energies(i) - energies(i-1)) > tol ) then
      j = j + 1
      unq(j) = energies(i)
      !deg(j-1) = 1
    else
      deg(j) = deg(j) + 1
    endif
  enddo
  deg(j+1:n) = 0

  !do i = 1, n
  !  print *, i, unq(i), deg(i)
  !enddo

  tot_deg = 0
  do i = 1, n
    if (deg(i)>1) then
      print *, i, unq(i), deg(i)
      tot_deg = tot_deg + deg(i)
    endif
  enddo
  print *, "Total fraction of states with degeneracies: ", real(tot_deg)/real(n)
  
end subroutine

program swap

  use exp_sparse
  use exponentiate
  use genmat
  use printing
  use MBL
  use omp_lib
  use sorts
  !use sorting
  !use stdlib_sorting_sort_index
  use iso_c_binding
  !use general
  implicit none

  complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, n_iterations
  integer (c_int)     ::  i, j, k, p, dim_Sz0
  integer (c_int)     ::  unit_ph

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  real(c_double) :: r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2
  
  real (c_double), dimension(:), allocatable :: r_avg, r_sq, r_avg2, r_sq2

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r, QE, E_MBL
  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap
  real (c_double), allocatable :: uE(:), QE_exact(:)
  integer (c_int), allocatable :: deg(:)

  real (c_double), allocatable :: log_avg(:), log_sq(:), log_pair_avg(:), log_near_avg(:)
  real (c_double) :: log_dis_avg, log_dis_sigma, log_pair_dis_avg, log_near_dis_avg

  integer (c_int), allocatable :: idx(:), config(:), states(:)

  logical :: SELECT
  EXTERNAL SELECT
  EXTERNAL find_degeneracies

  integer(c_int) :: count_beginning, count_end, count_rate, count1, count2
  character(len=200) :: filestring



  !Parametri Modello: J, V, h_z, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin
  dim_Sz0 = binom(nspin,nspin/2)

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
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

    !Coefficienti dello stato iniziale |psi> = (alpha|up>+beta|down>)^L

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring,93) "data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Quasi_Energies_nspin", &
    & nspin, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"
  !open(newunit=unit_ph,file=filestring)

  !91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  !92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  93  format(A,I0, A,F4.2, A,I0, A,F4.2, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)
 
  !------------------------------------------------

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  call buildSz0_HSwap(nspin, dim_Sz0, H)
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim_Sz0, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_z(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))
  allocate(U(dim_Sz0,dim_Sz0))
  !allocate(idx(dim_Sz0))

  !Allocate observables and averages
  allocate( r_avg(n_iterations), r_sq(n_iterations))
  allocate( r_avg2(n_iterations), r_sq2(n_iterations))
  
  !Allocate gap vectors 
  allocate( log_avg(n_iterations), log_sq(n_iterations), log_pair_avg(n_iterations), log_near_avg(n_iterations))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))
  allocate(E_MBL(n_iterations,dim_Sz0),QE(n_iterations,dim_Sz0))

  allocate(config(nspin), states(dim_Sz0))

  allocate(uE(dim_Sz0), deg(dim_Sz0), QE_exact(dim_Sz0))


  call zero_mag_states(nspin, dim_Sz0, states)

  r_avg = 0
  r_sq = 0
  r_avg2 = 0
  r_sq2 = 0

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(iteration, h_z, Vint, H, E, W_r, &
  !$OMP & U, W, PH)
  do iteration = 1, n_iterations
    
    if (n_iterations < 10) then
      print *, "iteration = ", iteration
    else
      if (mod(iteration, n_iterations/10)==0) then 
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
  
    !write (*,*) "Jint = ", Jint(:)
    !write (*,*) "Vint = ", Vint(:)
    !write (*,*) "h_z = ", h_z(:)
    !print *, ""
  
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    E_MBL(iteration,:) = E
    call gap_ratio(dim_Sz0, E, r_avg2(iteration), r_sq2(iteration))
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    E = real(C_UNIT*log(PH))
    call dpquicksort(E)
    call gap_ratio(dim_Sz0, E, r_avg(iteration), r_sq(iteration))
    QE(iteration,:) = E

    call exact_quasi_energies_Sz0(nspin, dim_Sz0, Vint, h_z, QE_exact)

    call find_degeneracies( size(E), E, uE, deg)
    print *, sum(deg), dim_Sz0

    call log_gap_difference(dim_Sz0, E, log_pair_avg(iteration), log_near_avg(iteration), log_avg(iteration), log_sq(iteration))
    !print *, iteration, log_pair_avg(iteration), log_near_avg(iteration), log_avg(iteration), log_sq(iteration)

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  !print *, E(:)
  do iteration = 1, n_iterations
    do i = 1, dim_Sz0
      !write (unit_ph,*) iteration, QE(iteration,i), E_MBL(iteration,i)
      !write (*,*) iteration, QE(iteration,i), E_MBL(iteration,i)
      print *, i, QE(iteration,i)
    enddo
  enddo
  !print *, sum(deg), dim_Sz0

  r_dis_avg = sum(r_avg) / n_iterations
  r_dis_sigma = sqrt( ( sum(r_sq)/n_iterations - r_dis_avg**2 ) / n_iterations )
  r_dis_avg2 = sum(r_avg2) / n_iterations
  r_dis_sigma2 = sqrt( ( sum(r_sq2)/n_iterations - r_dis_avg2**2 ) / n_iterations )

  log_dis_avg = sum(log_avg) / n_iterations
  log_dis_sigma = sqrt( ( sum(log_sq)/n_iterations - log_dis_avg**2 ) / n_iterations )
  log_pair_dis_avg = sum(log_pair_avg)/n_iterations
  log_near_dis_avg = sum(log_near_avg)/n_iterations


  print *, "Average and Variance of Gap Ratio (over the spectrum and then disorder)"
  print *, r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2

  print *, "Average of Difference of Logarithm of Gap log(Delta^alpha / Delta_0^alpha)"
  print *, log_dis_avg, log_dis_sigma, log_pair_dis_avg, log_near_dis_avg


  deallocate(Jint, Vint, h_z)
  deallocate(H,E,W_r)
  deallocate(USwap)
  deallocate(U, PH, W)
  
  !close(unit_ph)

  call take_time(count_rate, count_beginning, count1, 'T', "Program")

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end








