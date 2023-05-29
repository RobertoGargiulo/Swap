program test_LR

  use functions, only: zero_mag_states, decode, binom, &
    & init_random_seed, norm, disorder_average
  use exponentiate, only: diagSYM, expSYM, diagUN
  use observables, only: local_imbalance => local_imbalance_Sz0, &
    & sigmaz_corr => sigmaz_tot_corr_Sz0, IPR
  use matrices, only: buildHSwap => buildSz0_HSwap, buildHMBL => buildSz0_HMBL_LR
  use printing, only: take_time, printmat
  use omp_lib
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  implicit none

  complex (dcp), parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), parameter :: C_ONE = dcmplx(1._dp, 0._dp)
  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  real (dp), parameter :: pi = 4._dp * atan(1._dp)
  real (dp), parameter :: tol = 1.0e-8

  integer (ip)     ::  nspin, dim, n_disorder
  integer (ip)     ::  i, j, l, r, k, p, q, dim_Sz0
  integer (ip)     ::  unit_corr

  real (dp), dimension(:), allocatable :: Jxy, hz
  real (dp), dimension(:,:), allocatable :: Vzz
  real (dp) :: T0, T1, Jxy_coupling, Vzz_coupling, hz_coupling, kick, alpha
  real (dp) :: r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2
  
  real (dp), dimension(:), allocatable :: r_avg, r_sq, r_avg2, r_sq2

  integer (ip), allocatable :: deg(:), idxuE(:)
  real (dp), dimension(:), allocatable :: E, QE_exact, E_exact
  real (dp), dimension(:,:), allocatable :: H, W_r, QE
  complex (dcp), dimension(:), allocatable :: PH, psi_Sz0
  complex (dcp), dimension(:,:), allocatable :: U, W, USwap

  real (dp), allocatable :: CORR(:,:)
  real (dp), allocatable, dimension(:,:) :: LI, IPR_arr
  real (dp) :: CORR_avg, CORR_sigma

  integer (ip), allocatable :: idx(:), config(:), states(:)

  logical :: SELECT
  EXTERNAL SELECT

  integer(ip) :: count_beginning, count_end, count_rate, count1, count2
  character(len=200) :: filestring


  !Parametri Modello: J, V, hz, T0, T1/epsilon, nspin/L
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

  write (*,*) "Transverse Interaction Constant -J * (XX + YY)"
  read (*,*) Jxy_coupling
  print*,""

  write (*,*) "Longitudinal Interaction Constant -V * ZZ"
  read (*,*) Vzz_coupling
  print*,""

  write (*,*) "Longitudinal Field hz * Z"
  read (*,*) hz_coupling
  print*,""

  !---Read below for distributions of J, V, hz

  write (*,*) "Perturbation on Kick T1 = pi/4 + kick"
  read (*,*) kick
  print*,""

  write (*,*) "Power-Law Coefficient Vzz = V_{ij} / |i-j|^alpha"
  read (*,*) alpha
  print*,""
 
  T1 = pi/4 + kick

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring,93) "data/eigen/correlations_Swap_LR_Sz0_nspin", &
    & nspin, "_period", T0, "_n_disorder", n_disorder, &
    & "_Jxy", Jxy_coupling, "_Vzz", int(Vzz_coupling), Vzz_coupling-int(Vzz_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, &
    & "_alpha", int(alpha), alpha-int(alpha), ".txt"
  open(newunit=unit_corr, file=filestring)

  !91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  !92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  93  format(A,I0, A,F4.2, A,I0, A,F7.5, A,I0,F0.2, A,I0,F0.2, A,F5.3, A,I0,F0.2, A)
 
  !------------------------------------------------

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  call buildHSwap(nspin, dim_Sz0, H)
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim_Sz0, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jxy(nspin-1), Vzz(nspin-1,nspin), hz(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))
  allocate(U(dim_Sz0,dim_Sz0))
  !allocate(idx(dim_Sz0))

  allocate(config(nspin), states(dim_Sz0))
  call zero_mag_states(nspin, dim_Sz0, states)

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))
  allocate(QE(n_disorder,dim_Sz0))
  allocate(psi_Sz0(dim_Sz0))

  !Allocate for Entanglement
  allocate(LI(n_disorder,dim_Sz0), IPR_arr(n_disorder,dim_Sz0))
  allocate(CORR(n_disorder,dim_Sz0))

  !$OMP PARALLEL
  call init_random_seed() 
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "Current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(i, hz, Vzz, H, E, W_r, &
  !$OMP & U, W, PH, psi_Sz0, r, l)
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
 
    !call random_number(Jxy)
    !Jxy = 2*Jxy_coupling*(Jxy - 0.5) !Jxy in [-J,J]
    Jxy = -Jxy_coupling

    !Vzz = -Vzz_coupling
    call random_number(Vzz)
    Vzz = -Vzz_coupling + Vzz_coupling*(Vzz - 0.5) !Vzz in [-V-V/2,-V+V/2]
    !Vzz = 1
    do k = 1, nspin-1
      Vzz(k,1:k) = 0
      do q = k+1, nspin
        Vzz(k,q) = Vzz(k,q) / ( abs(k-q)**alpha )
      enddo
      !write (*,*) Vzz(k,:)
    enddo
    Vzz = Vzz / norm(alpha, nspin)
 
    call random_number(hz)
    hz = 2*hz_coupling*(hz-0.5) !hz in [-hz_coupling, hz_coupling]
 
    !write (*,*) "Jxy = ", Jxy(:)
    !write (*,*) "Vzz = ", Vzz(:)
    !write (*,*) "hz = ", hz(:)
    !print *, ""
 
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildHMBL( nspin, dim_Sz0, Jxy, Vzz, hz, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    E = real(C_UNIT*log(PH), kind=dp)
    QE(i,:) = E

    do l = 1, dim_Sz0

      psi_Sz0 = W(1:dim_Sz0,l)

      LI(i,l) = local_imbalance(nspin, dim_Sz0, psi_Sz0)
      IPR_arr(i,l) = IPR(psi_Sz0)
      CORR(i,l) = sigmaz_corr(nspin, dim_Sz0, psi_Sz0)
      

      !---- Print to output -------!!!!
      !print "(A8,4X,A3,4X,A)", "|c_i|^2", "l", "config"
      !do r = 1, dim_Sz0
      !  if( abs(psi_Sz0(r))**2 > 1.0e-4) then
      !    j = states(r)
      !    call decode(j, nspin, config)
      !    print 1, abs(psi_Sz0(r))**2, l, config(:)
      !  endif
      !enddo
      !1 format (F8.5, 4X,I3, 4X,*(I0))

      !write (*, "(A12,*(A26))") "Dis. Realization", "CORR_Z", "LI", "IPR"
      !print *, i, CORR(i,l), LI(i,l), IPR_arr(i,l)
      !print*, ""

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  call write_info(unit_corr)
  !write (*, "(A12,*(A26))") "Dis. Realization", "CORR_Z", "LI", "IPR"
  write (unit_corr, "(A12,*(A26))") "Dis. Realization", "CORR_Z", "LI", "IPR"
  do i = 1, n_disorder
    do l = 1, dim_Sz0
      !print *, i, CORR(i,l), LI(i,l), IPR_arr(i,l), LI(i,l) * nspin * (LI(i,l) * nspin - 1) / 2
      write (unit_corr, *) i, CORR(i,l), LI(i,l), IPR_arr(i,l)
    enddo
  enddo

  !CORR_avg = sum(CORR(:,

  !call disorder_average( CORR(:), CORR(:)**2, CORR_avg, CORR_sigma)
  print *, "Total average correlations:"
  print *,  sum(CORR) / size(CORR) !, sum( LI * nspin * (LI*nspin - 1) / 2 ) / size(LI)
  print *, "Total average normalized correlations:"
  print *,  sum( CORR / ( (LI*nspin)*(LI*nspin-1)/2) ) / size(CORR)

  deallocate(Jxy, Vzz, hz)
  deallocate(H,E,W_r)
  deallocate(USwap)
  deallocate(U, PH, W)
 
  close(unit_corr)

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
  write (unit_file,*) "Eigenvectors of Floquet Operator U_F = U_swap e^(-i H)."
  write (unit_file,*) "Spin-1/2 chain with hamiltonian H = sum hz * Z + V * ZZ + J * (XX + YY)."
  write (unit_file,*) "Periodic perturbed swap, U_swap = exp(-i(pi/4 + kick) * sum (sigma*sigma) )."
  write (unit_file,*) "V = V_{ij}/|i-j|^alpha, with V_{ij} is taken in [-3V/2, -V/2];  hz is taken in [-hz, hz];  J is uniform."
  write (unit_file,*) "Exact diagonalization of the dense matrix has been used to compute U_F and diagonalize it."

end subroutine
