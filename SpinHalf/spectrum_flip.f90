program flip

  use functions, only: init_random_seed
  use exponentiate, only: diagSYM, expSYM, diagUN
  use observables, only: gap_ratio, spectral_pairing => log_gap_difference !, &
    !& exact_QE => exact_quasi_energies_Sz0_LR, exact_E => exact_energies_Sz0_LR
  use matrices, only: buildHFlip, buildHMBL => buildHNayak
  use printing, only: take_time, printmat
  use sorts, only: sort => dpquicksort
  use omp_lib
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  implicit none

  complex (dcp), parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), parameter :: C_ONE = dcmplx(1._dp, 0._dp)
  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  real (dp), parameter :: pi = 4.d0 * datan(1.d0)

  integer (ip)     ::  nspin, dim, n_disorder
  integer (ip)     ::  i, j, l, k, p, q
  integer (ip)     ::  unit_ph

  real(dp), dimension(:), allocatable :: Vzz, hz, hx
  real(dp) :: T0, T1, Vzz_coupling, hz_coupling, hx_coupling, kick
  real(dp) :: r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2
  
  real (dp), dimension(:), allocatable :: r_avg, r_sq, r_avg2, r_sq2

  integer (ip), allocatable :: deg(:), idxuE(:)
  real (dp), dimension(:), allocatable :: E, QE_exact, E_exact
  real (dp), dimension(:,:), allocatable :: H, W_r, QE, E_MBL
  complex(dcp), dimension(:), allocatable :: PH
  complex(dcp), dimension(:,:), allocatable :: U, W, UFlip

  real (dp), allocatable :: log_avg(:), log_sq(:), log_pair_avg(:), log_near_avg(:)
  real (dp) :: log_dis_avg, log_dis_sigma, log_pair_dis_avg, log_near_dis_avg

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

  write (*,*) "Longitudinal Interaction Constant -V * ZZ"
  read (*,*) Vzz_coupling
  print*,""

  write (*,*) "Longitudinal Field hz * Z"
  read (*,*) hz_coupling
  print*,""

  !---Read below for distributions of hx, V, hz

  write (*,*) "Perturbation on Kick T1 = pi/2 + kick"
  read (*,*) kick
  print*,""

  T1 = pi/2 + kick

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring,93) "data/eigen/quasienergies_Flip_nspin", &
    & nspin, "_period", T0, "_n_disorder", n_disorder, &
    & "_hx", int(hx_coupling), hx_coupling-int(hx_coupling), "_Vzz", int(Vzz_coupling), Vzz_coupling-int(Vzz_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"
  open(newunit=unit_ph,file=filestring)

  !91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  !92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  93  format(A,I0, A,F4.2, A,I0, A,I0,F0.3, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)
 
  !------------------------------------------------

  !BUILD DRIVING PROTOCOL (NO DISORDER) UFlip = exp(-i*(pi/4 + eps)*HFlip)
  allocate(H(dim,dim), E(dim), W_r(dim,dim), UFlip(dim,dim))
  call buildHFlip(nspin, dim, H)
  call diagSYM( 'V', dim, H, E, W_r)
!  print *, "HFlip = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim, -C_UNIT*T1, E, W_r, UFlip )
  deallocate(E, W_r)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( hx(nspin), Vzz(nspin-1), hz(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(H(dim,dim), E(dim), W_r(dim,dim))
  allocate(U(dim,dim))
  !allocate(idx(dim))

  !Allocate observables and averages
  allocate( r_avg(n_disorder), r_sq(n_disorder))
  allocate( r_avg2(n_disorder), r_sq2(n_disorder))
  
  !Allocate gap vectors 
  allocate( log_avg(n_disorder), log_sq(n_disorder), log_pair_avg(n_disorder), log_near_avg(n_disorder))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim))
  allocate(W(dim,dim))
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
  !$OMP DO private(i, hx, hz, Vzz, H, E, W_r, &
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
    hx = hx_coupling * hx !hx in [0,hx]

    !Vzz = -Vzz_coupling
    call random_number(Vzz)
    Vzz = -Vzz_coupling + Vzz_coupling*(Vzz - 0.5) !Vzz in [-V,V]
 
    call random_number(hz)
    hz = 2*hz_coupling*(hz-0.5) !hz in [-hz_coupling, hz_coupling]
 
    !write (*,*) "hx = ", hx(:)
    !write (*,*) "Vzz = ", Vzz(:)
    !write (*,*) "hz = ", hz(:)
    !print *, ""
 
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildHMBL( nspin, dim, hx, Vzz, hz, H )
    call diagSYM( 'V', dim, H, E, W_r )
    E_MBL(i,:) = E

    call gap_ratio(dim, E, r_avg2(i), r_sq2(i))

    call expSYM( dim, -C_UNIT*T0, E, W_r, U )
    U = matmul(UFlip,U)
    call diagUN( SELECT, dim, U, PH, W)
    E = real(C_UNIT*log(PH), kind(dp))
    call sort(E)
    QE(i,:) = E

    !print *, "Degeneracies of QE:"
    !call find_degeneracies( size(E), E, idxuE, deg) 
    !print *, sum(deg), dim

    !QE_exact = exact_QE(nspin, dim, Vzz, hz)
    !call sort(QE_exact)
    !print *, "Degeneracies of QE_exact:"
    !call find_degeneracies( size(QE_exact), QE_exact, idxuE, deg) 
    !print *, sum(deg), dim

    !E_exact = exact_E(nspin, dim, Vzz, hz)
    !call sort(E_exact)



    !print *, "Quasienergies:"
    !print "(*(A26))", "Disorder Realization", "l", "QE", "QE_exact", "E_MBL"
    !do l = 1, dim
    !  write (*,*) i, l, QE(i,l), QE_exact(l), E_MBL(i,l), E_exact(l)
    !enddo


    call gap_ratio(dim, E, r_avg(i), r_sq(i))
    call spectral_pairing(dim, E, log_pair_avg(i), log_near_avg(i), log_avg(i), log_sq(i))
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

  r_dis_avg = sum(r_avg) / n_disorder
  r_dis_sigma = sqrt( ( sum(r_sq)/n_disorder - r_dis_avg**2 ) / n_disorder )
  r_dis_avg2 = sum(r_avg2) / n_disorder
  r_dis_sigma2 = sqrt( ( sum(r_sq2)/n_disorder - r_dis_avg2**2 ) / n_disorder )

  log_dis_avg = sum(log_avg) / n_disorder
  log_dis_sigma = sqrt( ( sum(log_sq)/n_disorder - log_dis_avg**2 ) / n_disorder )
  log_pair_dis_avg = sum(log_pair_avg)/n_disorder
  log_near_dis_avg = sum(log_near_avg)/n_disorder


  print *, "Average of Gap Ratio (over the spectrum and then disorder)"
  print "(*(A26))", "<r>_Flip", "sigma(r)_Flip", "<r>_MBL", "sigma(r)_MBL"
  print *, r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2

  print *, "Average of Logarithmic (pi-)Gap"
  print "(*(A26))", "<log(Delta_pi/Delta_0)>_Flip", "sigma(log(Delta_pi/Delta_0))_Flip", & 
    & "<log(Delta_pi)>_Flip", "<log(Delta_0)>_Flip"
  print *, log_dis_avg, log_dis_sigma, log_pair_dis_avg, log_near_dis_avg


  deallocate(hx, Vzz, hz)
  deallocate(H,E,W_r)
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
  write (unit_file,*) "Spin-1/2 chain with hamiltonian H = sum hz * Z + V * ZZ + hx * X."
  write (unit_file,*) "Periodic perturbed flip, U_flip = exp(-i(pi/2 + kick) * sum X )."
  write (unit_file,*) "V is taken in [-3V/2, -V/2];  hz is taken in [-hz, hz];  hz is taken in [0, hx]."
  write (unit_file,*) "Exact diagonalization of the dense matrix has been used to compute U_F and diagonalize it."

end subroutine

subroutine find_degeneracies( n, energies, unq, deg )
  !Finds all degeneracies of a real ordered 1D array of length n (up to a fixed tolerance)
  !idxuE contains the indices of the unique values of E
  !deg contains the corresponding degeneracy

  use iso_c_binding
  implicit none

  integer (c_int), intent(in) :: n
  real (c_double), intent(in) :: energies(n) !'energies' has to be sorted already
  !real (c_double), intent(out)  :: unq(n)
  integer (c_int), intent(out)  :: deg(n), unq(n)

  integer :: i, j, tot_deg
  real (c_double) :: tol = 1.0e-10


  unq = 0
  deg = 1

  j = 1
  unq(j) = 1 !energies(j)
  do i = 2, n
    if ( abs(energies(i) - energies(i-1)) > tol ) then
      j = j + 1
      !unq(j) = energies(i)
      unq(j) = i
    else
      deg(j) = deg(j) + 1
    endif
  enddo
  deg(j+1:n) = 0

  tot_deg = 0
  !print *, "Degenerate quasienergies: "
  !print "(*(A15))", "i", "idxunique(i)", "unique(idx)", "deg(i)"
  do i = 1, n
    if (deg(i)>0) then
      !print *, i, unq(i), energies(unq(i)), deg(i)
      tot_deg = tot_deg + merge(deg(i), 0, deg(i) > 1)
    endif
  enddo
  print *, "Total fraction of states with degeneracies: ", real(tot_deg)/real(n)
  print *, "Maximum degeneracy of a given energy: ", maxval(deg)
  print *, ""

end subroutine
