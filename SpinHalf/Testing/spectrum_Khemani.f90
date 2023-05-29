program flip

  use functions, only: init_random_seed, find_degeneracies, disorder_average
  use exponentiate, only: diagSYM, expSYM, diagUN, diagHE, expHE
  use observables, only: gap_ratio, spectral_pairing => log_gap_difference, &
    & shift_spectral_pairing => log_gap_difference_half_spectrum_shift
    !& exact_QE => exact_quasi_energies_Sz0_LR, exact_E => exact_energies_Sz0_LR
  use matrices, only: buildHFlip, buildHMBL => buildHKhemani
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
  integer (ip)     ::  i, j, l, k, p, q
  integer (ip)     ::  unit_ph

  real(dp), dimension(:), allocatable :: Vzz, hz, hx, hy
  real(dp) :: T0, T1, Vzz_coupling, hz_coupling, hx_coupling, hy_coupling, kick
  real(dp) :: r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2
  
  real (dp), dimension(:), allocatable :: r_avg, r_sq, r_avg2, r_sq2

  integer (ip), allocatable :: deg(:), idxuE(:)
  real (dp), dimension(:), allocatable :: E, QE_exact, E_exact
  real (dp), dimension(:,:), allocatable :: H_r, W_r, QE, E_MBL
  complex(dcp), dimension(:), allocatable :: PH
  complex(dcp), dimension(:,:), allocatable :: U, W, UFlip, H_c

  real (dp), allocatable :: log_avg(:), log_sq(:), log_pair_avg(:), log_near_avg(:), &
    & log_pair_sq(:), log_near_sq(:)
  real (dp) :: log_dis_avg, log_dis_sigma, log_pair_dis_avg, log_near_dis_avg, &
    & log_pair_dis_sigma, log_near_dis_sigma
  real (dp), allocatable :: shift_log_avg(:), shift_log_sq(:), shift_log_pair_avg(:), shift_log_near_avg(:), &
    & shift_log_pair_sq(:), shift_log_near_sq(:)
  real (dp) :: shift_log_dis_avg, shift_log_dis_sigma, shift_log_pair_dis_avg, shift_log_near_dis_avg, &
    & shift_log_pair_dis_sigma, shift_log_near_dis_sigma

  integer (ip), allocatable :: idx(:), config(:), states(:)

  logical :: SELECT
  EXTERNAL SELECT

  integer(ip) :: count_beginning, count_end, count_rate, count1, count2
  character(len=200) :: filestring


  !Parametri Modello: hx, V, hz, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin

  write (*,*) "Number of Disorder Realizations"
  read (*,*) n_disorder
  print*,""

  write (*,*) "Period T0"
  read (*,*) T0
  print*,""

  write (*,*) "Transverse Field hx * X"
  read (*,*) hx_coupling
  print*,""

  write (*,*) "Transverse Field hy * Y"
  read (*,*) hy_coupling
  print*,""

  write (*,*) "Longitudinal Field hz * Z"
  read (*,*) hz_coupling
  print*,""

  write (*,*) "Longitudinal Interaction Constant -V * ZZ"
  read (*,*) Vzz_coupling
  print*,""

  !---Read below for distributions of hx, V, hz

  write (*,*) "Perturbation on Kick T1 = pi/2 + kick"
  read (*,*) kick
  print*,""

  T1 = pi/2 + kick

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring,93) "data/eigen/quasienergies_Khemani_nspin", &
    & nspin, "_period", T0, "_n_disorder", n_disorder, &
    & "_hx", int(hx_coupling), hx_coupling-int(hx_coupling), "_hy", int(hy_coupling), hy_coupling-int(hy_coupling),&
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), &
    & "_Vzz", int(Vzz_coupling), Vzz_coupling-int(Vzz_coupling), "_kick", kick, ".txt"
  open(newunit=unit_ph,file=filestring)

  !91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  !92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  93  format(A,I0, A,F4.2, A,I0, A,I0,F0.4, A,I0,F0.4, A,I0,F0.4, A,I0,F0.3, A,F5.3, A)
 
  !------------------------------------------------

  !BUILD DRIVING PROTOCOL (NO DISORDER) UFlip = exp(-i*(pi/4 + eps)*HFlip)
  allocate(H_r(dim,dim), E(dim), W_r(dim,dim), UFlip(dim,dim))
  call buildHFlip(nspin, dim, H_r)
  call diagSYM( 'V', dim, H_r, E, W_r)
!  print *, "HFlip = "
!  call printmat(dim, H, 'R')
  deallocate(H_r)
  call expSYM( dim, -C_UNIT*T1, E, W_r, UFlip )
  deallocate(E, W_r)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( hx(nspin), hy(nspin), Vzz(nspin-1), hz(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(H_c(dim,dim), E(dim), W(dim,dim))
  allocate(U(dim,dim))
  !allocate(idx(dim))

  !Allocate observables and averages
  allocate( r_avg(n_disorder), r_sq(n_disorder))
  allocate( r_avg2(n_disorder), r_sq2(n_disorder))
  
  !Allocate gap vectors 
  allocate( log_avg(n_disorder), log_sq(n_disorder), log_pair_avg(n_disorder), log_near_avg(n_disorder))
  allocate( shift_log_avg(n_disorder), shift_log_sq(n_disorder), shift_log_pair_avg(n_disorder), shift_log_near_avg(n_disorder))
  allocate( log_pair_sq(n_disorder), log_near_sq(n_disorder), shift_log_pair_sq(n_disorder), shift_log_near_sq(n_disorder) )

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim))
  allocate(E_MBL(n_disorder,dim),QE(n_disorder,dim))

  allocate(idxuE(dim), deg(dim), QE_exact(dim), E_exact(dim))

  r_avg = 0
  r_sq = 0
  r_avg2 = 0
  r_sq2 = 0

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(i, hx, hy, hz, Vzz, H_c, E, &
  !$OMP & U, W, PH)
  do i = 1, n_disorder
 
    if (n_disorder < 10) then
      print *, "Disorder Realization = ", i
    else
      if (mod(i, n_disorder/10)==0) then 
        print *, "Disorder Realization = ", i
      endif
    endif

    !-------------------------------------------------
    !PARAMETERS
 
    call random_number(hx)
    hx = hx_coupling + 2*hx_coupling*(hx-0.5)
    call random_number(hy)
    hy = hy_coupling + 2*hy_coupling*(hy-0.5)
    call random_number(hz)
    hz = hz_coupling + 2*hz_coupling*(hz-0.5) !hz in [0, 2*hz_coupling]

    call random_number(Vzz)
    Vzz = Vzz_coupling + Vzz_coupling*(Vzz - 0.5) !Vzz in [V-V/2,V+V/2]
 
 
    !write (*,*) "hx = ", hx(:)
    !write (*,*) "hy = ", hy(:)
    !write (*,*) "hz = ", hz(:)
    !write (*,*) "Vzz = ", Vzz(:)
    !print *, ""
 
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildHMBL( nspin, dim, Vzz, hx, hy, hz, H_c )
    call diagHE( 'V', dim, H_c, E, W )
    E_MBL(i,:) = E

    call gap_ratio(E, r_avg2(i), r_sq2(i))

    call expHE( dim, -C_UNIT*T0, E, W, U )
    U = matmul(UFlip,U)
    call diagUN( SELECT, dim, U, PH, W)
    E = real(C_UNIT*log(PH), kind=dp)
    call sort(E)
    QE(i,:) = E

    !print *, "Degeneracies of QE:"
    !call find_degeneracies( size(E), E, idxuE, deg) 
    !print *, sum(deg, deg>1), dim

    call gap_ratio(E, r_avg(i), r_sq(i))
    call spectral_pairing(E, log_pair_avg(i), log_pair_sq(i), log_near_avg(i), log_near_sq(i), log_avg(i), log_sq(i))
    call shift_spectral_pairing(E, shift_log_pair_avg(i), shift_log_pair_sq(i), &
      & shift_log_near_avg(i), shift_log_near_sq(i), shift_log_avg(i), shift_log_sq(i))
    !print *, i, log_pair_avg(i), log_near_avg(i), log_avg(i), log_sq(i)

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  call write_info(unit_ph)
  !print "(*(A26))", "Disorder Realization", "l", "QE", "E_MBL"
  write (unit_ph, "(2(A12),2(A26))") "Disorder Realization", "l", "QE", "E_MBL"
  do i = 1, n_disorder
    do l = 1, dim
      !write (*,*) i, l, QE(i,l), E_MBL(i,l)
      write (unit_ph,*) i, l, QE(i,l), E_MBL(i,l)
    enddo
  enddo

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

  deallocate(hx, Vzz, hz)
  deallocate(H_c,E)
  deallocate(UFlip)
  deallocate(U, PH, W)
 
  close(unit_ph)

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
  write (unit_file,*) "Quasi-Energies of Floquet Operator U_F = U_flip e^(-i H)."
  write (unit_file,*) "Spin-1/2 chain with hamiltonian H = sum V * ZZ + hx * X + hy * Y + hz * Z."
  write (unit_file,*) "Periodic perturbed flip, U_flip = exp(-i(pi/2 + kick) * sum X )."
  write (unit_file,*) "V is taken in [V/2, 3V/2];  the h_i are taken in [0, 2*h_i]."
  write (unit_file,*) "Exact diagonalization of the dense matrix has been used to compute U_F and diagonalize it."

end subroutine
