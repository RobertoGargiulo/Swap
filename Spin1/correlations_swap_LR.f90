program test_LR

  use functions, only: dimSpin1_Sz, decode, zero_mag_states, &
    & init_random_seed, norm => normalization_power_law
  use exponentiate, only: diagSYM, expSYM, diagUN
  use observables, only: local_imbalance => local_imbalance_Sz0, &
    & sigmaz_corr => sigmaz_tot_corr_Sz0, IPR, sigmaz2_corr => sigmaz2_tot_corr_Sz0
  use matrices, only: buildHSwap => buildSz_HSwap, buildHMBL => buildSz_HMBL_LR
  use printing, only: take_time, printmat
  use omp_lib
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex, ilp => c_long
  implicit none

  complex (dcp), parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), parameter :: C_ONE = dcmplx(1._dp, 0._dp)
  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  real (dp), parameter :: pi = 4._dp * atan(1._dp)
  real (dp), parameter :: tol = 1.0e-8
  integer (ip), parameter :: dimSpin1 = 3

  integer (ip)     ::  nspin, dim, n_disorder
  integer (ip)     ::  i, l, k, q, dim_Sz, Sz
  integer (ip)     ::  unit_corr

  real (dp), dimension(:), allocatable :: Jxy, hz
  real (dp), dimension(:,:), allocatable :: Vzz
  real (dp) :: T0, T1, Jxy_coupling, Vzz_coupling, hz_coupling, kick, alpha
  
  real (dp), dimension(:), allocatable :: E
  real (dp), dimension(:,:), allocatable :: H, W_r, QE
  complex (dcp), dimension(:), allocatable :: PH, psi_Sz0
  complex (dcp), dimension(:,:), allocatable :: U, W, USwap

  real (dp), allocatable :: CORR(:,:), CORR2(:,:)
  real (dp), allocatable, dimension(:,:) :: LI, IPR_arr

  integer (ip), allocatable :: config(:), states(:)

  logical :: SELECT
  EXTERNAL SELECT
  EXTERNAL write_info

  integer(ilp) :: count_beginning, count_end, count_rate
  character(len=200) :: filestring


  !Parametri Modello: J, V, hz, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""

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

  write (*,*) "Perturbation on Kick T1 = pi/2 + kick"
  read (*,*) kick
  print*,""

  write (*,*) "Power-Law Coefficient Vzz = V_{ij} / |i-j|^alpha"
  read (*,*) alpha
  print*,""
 
  dim = dimSpin1**nspin
  T1 = pi/2 + kick
  Sz = 0 !!!!!!! <------- Magnetization of Subspace
  dim_Sz = dimSpin1_Sz(nspin, Sz)

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
  93 format(A,I0, A,F4.2, A,I0, A,F7.5, A,I0,F0.2, A,I0,F0.2, A,F5.3, A,I0,F0.2, A)
 
  !------------------------------------------------

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz,dim_Sz), E(dim_Sz), W_r(dim_Sz,dim_Sz), USwap(dim_Sz,dim_Sz))
  call buildHSwap(nspin, dim_Sz, Sz, H)
  call diagSYM( 'V', dim_Sz, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim_Sz, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jxy(nspin-1), Vzz(nspin-1,nspin), hz(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(H(dim_Sz,dim_Sz), E(dim_Sz), W_r(dim_Sz,dim_Sz))
  allocate(U(dim_Sz,dim_Sz))
  !allocate(idx(dim_Sz))

  allocate(config(nspin), states(dim_Sz))
  call zero_mag_states(nspin, dim_Sz, states)

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz))
  allocate(W(dim_Sz,dim_Sz))
  allocate(QE(n_disorder,dim_Sz))
  allocate(psi_Sz0(dim_Sz))

  !Allocate for Entanglement
  allocate(LI(n_disorder,dim_Sz), IPR_arr(n_disorder,dim_Sz))
  allocate(CORR(n_disorder,dim_Sz),CORR2(n_disorder,dim_Sz))

  !$OMP PARALLEL
  call init_random_seed() 
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "Current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(i, hz, Vzz, H, E, W_r, &
  !$OMP & U, W, PH, psi_Sz0, l, k, q)
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
    call buildHMBL( nspin, dim_Sz, Sz, Jxy, Vzz, hz, H )
    call diagSYM( 'V', dim_Sz, H, E, W_r )
    call expSYM( dim_Sz, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz, U, PH, W)
    E = real(C_UNIT*log(PH), kind=dp)
    QE(i,:) = E

    do l = 1, dim_Sz

      psi_Sz0 = W(1:dim_Sz,l)

      LI(i,l) = local_imbalance(nspin, dim_Sz, psi_Sz0)
      IPR_arr(i,l) = IPR(psi_Sz0)
      CORR(i,l) = sigmaz_corr(nspin, dim_Sz, psi_Sz0)
      CORR2(i,l) = sigmaz2_corr(nspin, dim_Sz, psi_Sz0)
      

      !---- Print to output -------!!!!
      !print "(A8,4X,A3,4X,A)", "|c_i|^2", "l", "config"
      !do r = 1, dim_Sz
      !  if( abs(psi_Sz0(r))**2 > 1.0e-4) then
      !    j = states(r)
      !    call decode(j, nspin, config)
      !    print 1, abs(psi_Sz0(r))**2, l, config(:)
      !  endif
      !enddo
      !1 format (F8.5, 4X,I3, 4X,*(I0))

      !write (*, "(A12,*(A26))") "Dis. Realization", "CORR_Z", "LI", "IPR"
      !write(*, "(1I12,*(1X,G25.16))"), i, CORR(i,l), CORR2(i,l), LI(i,l), IPR_arr(i,l)
      !print*, ""

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  call write_info(unit_corr)
  !write (*, "(A12,*(A26))") "Dis. Realization", "CORR_Z", "LI", "IPR"
  write (unit_corr, "(A12,*(A26))") "Dis. Realization", "CORR_Z", "CORR_Z^2", "LI", "IPR", "QE"
  do i = 1, n_disorder
    do l = 1, dim_Sz
      !write (*, *) i, CORR(i,l), LI(i,l), IPR_arr(i,l), QE(i,l)
      write (unit_corr, "(1I12,*(1X,G25.16))") i, CORR(i,l), CORR2(i,l), LI(i,l), IPR_arr(i,l), QE(i,l)
    enddo
  enddo


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

  write (unit_file,'(A)') "Some info: "
  write (unit_file,'(A)') "Correlations of the Eigenvectors of Floquet Operator U_F = U_swap e^(-i H)."
  write (unit_file,'(A)') "Specifically: Connected Correlations functions of S_k^z, S_k^z**2, Local Imbalance, IPR, Quasi-energies."
  write (unit_file,'(A)') "Spin-1/2 chain with hamiltonian H = sum hz * Z + V * ZZ + J * (XX + YY)."
  write (unit_file,'(A)') "Periodic perturbed swap, U_swap = exp(-i(pi/2 + kick) * sum (sigma*sigma)^2 + (sigma*sigma) )."
  write (unit_file,'(A)') "V = V_{ij}/|i-j|^alpha, with V_{ij} is taken in [-3V/2, -V/2];  hz is taken in [-hz, hz];  J is uniform."
  write (unit_file,'(A)') "Exact diagonalization of the dense matrix has been used to compute U_F and diagonalize it."

end subroutine
