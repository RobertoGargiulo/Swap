program test_LR

  use functions, only : init_random_seed, dimSpin1_Sz, &
    & project => projectState_FullHS_to_Sz, norm_V => normalization_power_law
  use exponentiate, only: diagSYM, expSYM
  use observables, only: sigmaz_Sz
  use matrices, only: buildHSwap => buildSz_HSwap, buildUMBL => buildSz_UMBL_LR, &
    & print_hamiltonian_Sz, print_unitary_Sz
  use printing, only: take_time, printmat
  use states, only: buildstate => buildUpZeroState, &
    & printstate_Sz, printstate
  use omp_lib
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex, ilp => c_long
  implicit none

  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  real (dp), parameter :: pi = 4._dp * atan(1._dp)
  !real (dp), parameter :: tol = 1.0e-8
  character(len=*), parameter :: name_initial_state = "UpZero"
  integer (ip), parameter :: dimSpin1 = 3

  integer (ip)     ::  nspin, dim, n_disorder
  integer (ip)     ::  i, k, q, t, dim_Sz, Sz
  real (dp) :: norm
  integer (ilp)     :: j, n_pow_periods, n_periods, steps

  real (dp), dimension(:), allocatable :: Jxy, hz
  real (dp), dimension(:,:), allocatable :: Vzz
  real (dp) :: T0, T1, Jxy_coupling, Vzz_coupling, hz_coupling, kick, alpha
  !complex (dcp) :: beta, delta

  real (dp), dimension(:), allocatable :: E
  real (dp), dimension(:,:), allocatable :: H, W_r
  complex (dcp), allocatable :: psi(:), psi_swap(:), psi_Sz(:)
  complex (dcp), dimension(:,:), allocatable :: U, USwap!, UMBL, UMBL2, Id

  integer(ip) :: count_beginning, count_end, count_rate

  character(len=200) :: filestring
  integer (ip) :: unit_decay

  !integer (ip)  :: n_periods
  real (dp), allocatable, dimension(:) :: sigmaz_previous, sigmaz_current, &
    & sigmaz_initial
  integer (ilp), allocatable, dimension(:) :: t_decay
  real (dp)   :: Z_previous, Z_current, t_decay_avg, sigma_t_decay
  integer (ip)   :: idecay, n_decays
  !real (dp) :: time1, time2

  EXTERNAL write_info
  interface
    pure function ones(n) result(arr)
      use iso_c_binding, dp => c_double, ip => c_int
      integer (ip), intent(in) :: n
      real (dp) :: arr(n)
    end function
  end interface

  !Parametri Modello: J, V, hz, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = dimSpin1**nspin

  write (*,*) "Number of Disorder Realizations"
  read (*,*) n_disorder
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  write (*,*) "Power of 2 for the Number of Periods at each step (n_per = 2^(pow)+1)"
  read (*,*) n_pow_periods 
  print*,""
  n_periods = 2**(n_pow_periods)
  j = n_periods * steps
  print *, "Total maximum time: ", j

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
 
  T1 = pi/2 + kick

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring,93) "data/dynamics/decay_sigmaz_Swap_LR_" // trim(name_initial_state) // "_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_n_disorder", n_disorder, "_n_periods", n_periods, &
    & "_Jxy", Jxy_coupling, "_Vzz", int(Vzz_coupling), Vzz_coupling-int(Vzz_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, &
    & "_alpha", int(alpha), alpha-int(alpha), ".txt"
  open(newunit=unit_decay, file=filestring)
  call write_info(unit_decay, trim(name_initial_state))

  93  format(A,I0, A,I0, A,F4.2, A,I0, A,I0, A,F7.5, A,I0,F0.2, A,I0,F0.2, A,F5.3, A,I0,F0.2, A)


  !------------- Initial State ------------------!
  allocate(psi(dim))
  call buildstate(nspin, dim, psi)
  call printstate(nspin, dim, psi, trim(name_initial_state))
  call project(nspin, dim, psi, dim_Sz, Sz, psi_Sz)
  deallocate(psi)
  call printstate_Sz(nspin, dim_Sz, Sz, psi_Sz, trim(name_initial_state))
  allocate(sigmaz_initial(nspin))
  sigmaz_initial = sigmaz_Sz(nspin, dim_Sz, Sz, psi_Sz)

 
  !------------------------------------------------

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz,dim_Sz), E(dim_Sz), W_r(dim_Sz,dim_Sz), USwap(dim_Sz,dim_Sz))
  call buildHSwap(nspin, dim_Sz, Sz, H)
  call diagSYM( 'V', dim_Sz, H, E, W_r)
  !print *, "HSwap = "
  !call print_hamiltonian_Sz(nspin, dim_Sz, Sz, H)
  !call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim_Sz, -C_UNIT*T1, E, W_r, USwap )
  !print *, "USwap = "
  !call print_unitary_Sz(nspin, dim_Sz, Sz, USwap)
  deallocate(E, W_r)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jxy(nspin-1), Vzz(nspin-1,nspin), hz(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H(dim_Sz,dim_Sz), E(dim_Sz), W_r(dim_Sz,dim_Sz))
  allocate(U(dim_Sz,dim_Sz))
  !allocate(UMBL(dim_Sz,dim_Sz), UMBL2(dim_Sz,dim_Sz))
  !allocate(Id(dim_Sz,dim_Sz))
  !Id = 0
  !forall (i=1:dim_Sz) Id(i,i) = 1

  !Allocate for Eigenvalues/Eigenvectors
  allocate(psi_swap(dim_Sz))

  !Allocate for Dynamics
  allocate( t_decay(n_disorder), sigmaz_previous(nspin), sigmaz_current(nspin) )
  n_decays = 0
  t_decay = steps

  print *, "Starting disorder loop."

  !$OMP PARALLEL
  call init_random_seed() 
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "Current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+: n_decays) private(i, j, t, hz, Vzz, norm, &
  !$OMP & psi_swap, U, idecay, &
  !$OMP & sigmaz_previous, sigmaz_current, Z_previous, Z_current )
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
 
    Jxy = -Jxy_coupling

    call random_number(Vzz)
    Vzz = -Vzz_coupling + Vzz_coupling*(Vzz - 0.5) !Vzz in [-V-V/2,-V+V/2]
    do k = 1, nspin-1
      Vzz(k,1:k) = 0
      do q = k+1, nspin
        Vzz(k,q) = Vzz(k,q) / ( abs(k-q)**alpha )
      enddo
    enddo
    Vzz = Vzz / norm_V(alpha, nspin)
 
    call random_number(hz)
    hz = 2*hz_coupling*(hz-0.5) !hz in [-hz_coupling, hz_coupling]
 
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR

    call buildUMBL( nspin, dim_Sz, Sz, Jxy, Vzz, hz, -C_UNIT*T0, U )
    U = matmul(USwap,U)
    do t = 1, n_pow_periods
      U = matmul(U,U)
    enddo
    !print *, "U = "
    !call print_unitary_Sz(nspin, dim_Sz, Sz, U)
    !print *, "Error UU^dagger - Id = "
    !print *, sum(abs( matmul(U,conjg(transpose(U))) - Id)) / size(U)


    psi_swap = psi_Sz
    j = 1
    idecay = 0
    do j = 2, steps

      sigmaz_previous = sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap)
      Z_previous = sum( sign(ones(nspin/2), sigmaz_initial(1::2) - sigmaz_initial(2::2)  ) * &
        & (sigmaz_previous(1::2) - sigmaz_previous(2::2)) )
      !write (*,"(1A, I0, 2X1A,*(G16.10))") "j = ", j-1, "sigma = ", sigmaz_previous(:)
      !print *, Z_previous
      !print *, ""
      !print *, sign(ones(nspin/2), sigmaz_initial(1::2) - sigmaz_initial(2::2)  ) * &
      !  & (sigmaz_previous(1::2) - sigmaz_previous(2::2))
      !print *, ""

      if(mod(j,10)==0) then
        norm = real(dot_product(psi_swap,psi_swap), kind=dp)
        psi_swap = psi_swap/sqrt(norm)
      endif
      psi_swap = matmul(U, psi_swap)
      sigmaz_current = sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap)
      Z_current = sum( sign(ones(nspin/2), sigmaz_initial(1::2) - sigmaz_initial(2::2)  ) * &
        & (sigmaz_current(1::2) - sigmaz_current(2::2)) )
      if (idecay==0) then
        if( (-1)**(n_periods) * Z_previous * Z_current < 0) then
          t_decay(i) = j * n_periods
          idecay = 1
          n_decays = n_decays + 1
          !print *, j, sigmaz_current
          !print *, Z_current
          !print *, sign(ones(nspin/2), sigmaz_initial(1::2) - sigmaz_initial(2::2)  ) * &
          !  & (sigmaz_current(1::2) - sigmaz_current(2::2))
          !print *, ""
          exit
        endif
      endif

    enddo
    !print *, ""

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

 
  t_decay_avg = sum(t_decay) / n_disorder
  sigma_t_decay = sqrt( (sum(t_decay**2)/n_disorder - t_decay_avg**2) / n_disorder)

  write(unit_decay,"(2(A26),A12)") "Dis. Realiz.", "t*"
  do i = 1, n_disorder
    write(unit_decay,*) i, t_decay(i)
  enddo
  write (unit_decay,"(2(A26),A12)") "<t*>", "sigma(t*)", "n_decays"
  write (unit_decay,*) t_decay_avg, sigma_t_decay, n_decays
  close (unit_decay)

  print *, "Average Decay Time and Errors"
  write (*,"(2(A26),A12)") "<t*>", "sigma(t*)", "n_decays"
  print*, t_decay_avg, sigma_t_decay, n_decays

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

end program

subroutine write_info(unit_file, state_name)

  integer, intent(in) :: unit_file
  character(len=*) :: state_name

  write (unit_file,*) "Some info: "
  write (unit_file,*) "Decay times of magnetization at integer multiples (2**n_pow_periods) of the period."
  write (unit_file,*) "Floquet Operator U_F = U_swap e^(-i H) and with multiplication UF -> UF * UF ."
  write (unit_file,*) "Spin-1 chain with hamiltonian H = sum hz * Z + V * ZZ + J * (XX + YY)."
  write (unit_file,*) "Periodic perturbed swap, U_swap = exp(-i(pi/2 + kick) * sum (sigma*sigma)^2 + (sigma*sigma) - 2 )."
  write (unit_file,*) "V = V_{ij}/|i-j|^alpha, with V_{ij} is taken in [-3V/2, -V/2] (OBC); &
    &  hz is taken in [-hz, hz];  J is uniform."
  write (unit_file,*) "Exact diagonalization of the dense matrix has been used to compute U_F."
  write (unit_file,*) "Initial state is "//trim(state_name)//trim(".")

end subroutine

pure function ones(n) result(arr)
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  integer (ip), intent(in) :: n
  real (dp) :: arr(n)

  arr = 1.0

end function
