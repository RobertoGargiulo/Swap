program exact_LR

  use functions, only: binom, init_random_seed, zero_mag_states, norm, disorder_average, &
    & find_degeneracies
  use exponentiate, only: diagSYM, expSYM, diagUN
  use observables, only: gap_ratio, spectral_pairing => log_gap_difference, &
    & shift_spectral_pairing => log_gap_difference_half_spectrum_shift, &
    & exact_QE => exact_quasi_energies_Sz0_LR, exact_E => exact_energies_Sz0_LR, &
    & gap_difference
  use matrices, only: buildHSwap => buildSz0_HSwap, buildHMBL => buildSz0_HMBL_LR
  use printing, only: take_time, printmat
  use sorts, only: sort => dpquicksort
  use omp_lib
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  implicit none

  complex (dcp), parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), parameter :: C_ONE = dcmplx(1._dp, 0._dp)
  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  real (dp), parameter :: pi = 4._dp * atan(1._dp)

  integer (ip)     ::  nspin, dim, n_disorder
  integer (ip)     ::  i, j, l, k, p, q, dim_Sz0
  integer (ip)     ::  unit_ph

  real(dp), dimension(:), allocatable :: hz
  real(dp), dimension(:,:), allocatable :: Vzz
  real(dp) :: T0, Vzz_coupling, hz_coupling, alpha
  real(dp) :: r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2
  
  real (dp), dimension(:), allocatable :: r_avg, r_sq, r_avg2, r_sq2

  integer (ip), allocatable :: deg(:), idxuE(:)
  real (dp), dimension(:), allocatable :: E
  real (dp), dimension(:,:), allocatable :: QE, E_MBL

  real (dp), allocatable :: log_avg(:), log_sq(:), log_pair_avg(:), log_near_avg(:), &
    & log_pair_sq(:), log_near_sq(:)
  real (dp) :: log_dis_avg, log_dis_sigma, log_pair_dis_avg, log_near_dis_avg, &
    & log_pair_dis_sigma, log_near_dis_sigma
  real (dp), allocatable :: shift_log_avg(:), shift_log_sq(:), shift_log_pair_avg(:), shift_log_near_avg(:), &
    & shift_log_pair_sq(:), shift_log_near_sq(:)
  real (dp) :: shift_log_dis_avg, shift_log_dis_sigma, shift_log_pair_dis_avg, shift_log_near_dis_avg, &
    & shift_log_pair_dis_sigma, shift_log_near_dis_sigma

  real (dp), allocatable :: pair_avg(:), near_avg(:)
  real (dp) :: pair_dis_avg, near_dis_avg

  integer (ip), allocatable :: idx(:), config(:), states(:), tot_deg(:)
  integer (ip) :: tot_deg_avg, max_deg, min_deg

  logical :: SELECT
  EXTERNAL SELECT

  integer(ip) :: count_beginning, count_end, count_rate
  character(len=200) :: filestring


  !Parametri Modello: V, hz, T0, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin
  dim_Sz0 = binom(nspin,nspin/2)

  write (*,*) "Number of Disorder Realizations"
  read (*,*) n_disorder
  print*,""

  write (*,*) "Period T0"
  read (*,*) T0
  print*,""

  write (*,*) "Longitudinal Interaction Constant -V * ZZ"
  read (*,*) Vzz_coupling
  print*,""

  write (*,*) "Longitudinal Field hz * Z"
  read (*,*) hz_coupling
  print*,""

  !---Read below for distributions of V, hz

  write (*,*) "Power-Law Coefficient Vzz = V_{ij} / |i-j|^alpha"
  read (*,*) alpha
  print*,""

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring,93) "data/eigen/quasienergies_exact_Swap_LR_Sz0_nspin", &
    & nspin, "_period", T0, "_n_disorder", n_disorder, &
    & "_Vzz", int(Vzz_coupling), Vzz_coupling-int(Vzz_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), &
    & "_alpha", int(alpha), alpha-int(alpha), ".txt"
  !open(newunit=unit_ph,file=filestring)

  !91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  !92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  93  format(A,I0, A,F4.2, A,I0, A,I0,F0.2, A,I0,F0.2, A,I0,F0.2, A)
 
  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Vzz(nspin-1,nspin), hz(nspin))

  !Allocate observables and averages
  allocate( r_avg(n_disorder), r_sq(n_disorder))
  allocate( r_avg2(n_disorder), r_sq2(n_disorder))
  
  !Allocate gap vectors 
  allocate( log_avg(n_disorder), log_sq(n_disorder), log_pair_avg(n_disorder), log_near_avg(n_disorder))
  allocate( shift_log_avg(n_disorder), shift_log_sq(n_disorder), shift_log_pair_avg(n_disorder), shift_log_near_avg(n_disorder))
  allocate( log_pair_sq(n_disorder), log_near_sq(n_disorder), shift_log_pair_sq(n_disorder), shift_log_near_sq(n_disorder) )
  allocate( pair_avg(n_disorder), near_avg(n_disorder))
  allocate( tot_deg(n_disorder) )

  !Allocate for Eigenvalues/Eigenvectors
  allocate(E(dim_Sz0))
  allocate(E_MBL(n_disorder,dim_Sz0),QE(n_disorder,dim_Sz0))

  allocate(config(nspin), states(dim_Sz0))

  allocate(idxuE(dim_Sz0), deg(dim_Sz0))

  call zero_mag_states(nspin, dim_Sz0, states)

  r_avg = 0
  r_sq = 0
  r_avg2 = 0
  r_sq2 = 0

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(i, hz, Vzz, E, idxuE, deg) 
  !, H, E, W_r, &
  !!$OMP & U, W, PH)
  do i = 1, n_disorder
 
    if (n_disorder < 10) then
      print *, "Disorder Realization = ", i
    else
      if (mod(i, n_disorder/10)==0) then 
        print *, "Disorder Realization = ", i
        !call take_time(count_rate, count_beginning, count1, 'F', filestring)
      endif
    endif

    !-------------------------------------------------
    !PARAMETERS
 
    !Vzz = -Vzz_coupling
    call random_number(Vzz)
    Vzz = -Vzz_coupling + Vzz_coupling*(Vzz - 0.5) !Vzz in [-V,V]
    !Vzz = 1
    do k = 1, nspin-1
      Vzz(k,1:k) = 0
      do q = k+1, nspin
        Vzz(k,q) = Vzz(k,q) / ( abs(k-q)**alpha * norm(alpha,nspin) )
      enddo
      !write (*,*) Vzz(k,:)
    enddo
    Vzz = Vzz / norm(alpha, nspin)
 
    call random_number(hz)
    hz = 2*hz_coupling*(hz-0.5) !hz in [-hz_coupling, hz_coupling]
 
    !write (*,*) "Vzz = ", Vzz(:)
    !write (*,*) "hz = ", hz(:)
    !print *, ""
 
    !---------------------------------------------------

    E = exact_E( nspin, dim_Sz0, Vzz, hz)
    call sort(E)
    E_MBL(i,:) = E
    call gap_ratio(E, r_avg2(i), r_sq2(i))

    E = exact_QE( nspin, dim_Sz0, Vzz, hz )
    call sort(E)
    QE(i,:) = E
    call gap_ratio(E, r_avg(i), r_sq(i))
    call spectral_pairing(E, log_pair_avg(i), log_pair_sq(i), log_near_avg(i), log_near_sq(i), log_avg(i), log_sq(i))
    call shift_spectral_pairing(E, shift_log_pair_avg(i), shift_log_pair_sq(i), &
      & shift_log_near_avg(i), shift_log_near_sq(i), shift_log_avg(i), shift_log_sq(i))

    !print *, "Degeneracies of QE_exact:"
    call find_degeneracies( size(E), E, idxuE, deg) 
    !print *, sum(deg, deg>1), dim_Sz0
    tot_deg(i) = sum(deg, deg>1)

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  !call write_info(unit_ph)
  !print "(*(A26))", "Disorder Realization", "l", "QE", "E_MBL"
  !write (unit_ph, "(2(A12),2(A26))") "Disorder Realization", "l", "QE", "E_MBL"
  !do i = 1, n_disorder
  !  do l = 1, dim_Sz0
  !    !write (*,*) i, l, QE(i,l), E_MBL(i,l)
  !    !write (unit_ph,*) i, l, QE(i,l), E_MBL(i,l)
  !  enddo
  !enddo

  call disorder_average(r_avg, r_sq, r_dis_avg, r_dis_sigma)
  call disorder_average(r_avg2, r_sq2, r_dis_avg2, r_dis_sigma2)

  call disorder_average(log_avg, log_sq, log_dis_avg, log_dis_sigma)
  call disorder_average(log_pair_avg, log_pair_sq, log_pair_dis_avg, log_pair_dis_sigma)
  call disorder_average(log_near_avg, log_near_sq, log_near_dis_avg, log_near_dis_sigma)

  call disorder_average( shift_log_avg, shift_log_sq, shift_log_dis_avg, shift_log_dis_sigma)
  call disorder_average( shift_log_pair_avg, shift_log_pair_sq, shift_log_pair_dis_avg, &
    & shift_log_pair_dis_sigma)
  call disorder_average( shift_log_near_avg, shift_log_near_sq, shift_log_near_dis_avg, &
    & shift_log_near_dis_sigma)

  tot_deg_avg = real(sum(tot_deg), kind=dp)/(n_disorder * dim_Sz0)
  min_deg = minval(tot_deg)
  max_deg = maxval(tot_deg)

  print *, "Average of Gap Ratio (over the spectrum and then disorder)"
  print "(*(A26))", "<r>_Flip", "sigma(r)_Flip", "<r>_MBL", "sigma(r)_MBL"
  print *, r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2

  print *, "Average of pi-Logarithmic (pi-)Gap"
  print "(*(A26))", "<log(Delta_pi/Delta_0)>", "sigma(log(Delta_pi/Delta_0))", & 
    & "<log(Delta_pi)>", "sigma(log(Delta_pi))", "<log(Delta_0)>", "sigma(log(Delta_0))"
  print *, log_dis_avg, log_dis_sigma, log_pair_dis_avg, log_pair_dis_sigma, &
   & log_near_dis_avg, log_near_dis_sigma

  print *, "Average of Half Shifted Logarithmic (pi-)Gap"
  print "(*(A26))", "<log(Delta_pi/Delta_0)>", "sigma(log(Delta_pi/Delta_0))", & 
    & "<log(Delta_pi)>", "sigma(log(Delta_pi))", "<log(Delta_0)>", "sigma(log(Delta_0))"
  print *, shift_log_dis_avg, shift_log_dis_sigma, shift_log_pair_dis_avg, &
    & shift_log_pair_dis_sigma, shift_log_near_dis_avg, shift_log_near_dis_sigma

  print *, "Average of Ordinary (pi-)Gap"
  print "(*(A26))", "<Delta_pi>_Swap", "<Delta_0>_Swap"
  print *, pair_dis_avg, near_dis_avg

  print *, "Average fraction of degenerate states (Sum of dimensions of degenerate eigenspaces)"
  print *, tot_deg_avg, min_deg, max_deg

  deallocate(Vzz, hz)
 
  !close(unit_ph)

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end

subroutine write_info(unit_file)

  integer, intent(in) :: unit_file

  write (unit_file,*) "Some info: "
  write (unit_file,*) "Exact Quasi-Energies of Floquet Operator U_F = U_swap e^(-i H)."
  write (unit_file,*) "Spin-1/2 chain with hamiltonian H = sum hz * Z + V * ZZ."
  write (unit_file,*) "Periodic swap, U_swap = exp(-i pi/4 * sum (sigma*sigma) )."
  write (unit_file,*) "V = V_{ij}/|i-j|^alpha, with V_{ij} is taken in [-3V/2, -V/2];  hz is taken in [-hz, hz]."

end subroutine
