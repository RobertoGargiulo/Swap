program test_LR

  use functions, only : init_random_seed, dimSpinHalf_Sz, &
    & project => projectState_FullHS_to_Sz, norm_V => normalization_power_law
  use exponentiate, only: diagSYM, expSYM
  use observables, only: sigmaz_Sz
  use matrices, only: buildHSwap => buildSz_HSwap, buildHMBL => buildSz_HMBL_XX_ZZ_LR, &
    & print_hamiltonian_Sz, print_unitary_Sz
  use printing, only: take_time, printmat
  use states, only: buildstate => buildNeelState, &
    & printstate_Sz, printstate
  use omp_lib
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  implicit none

  complex (dcp), parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), parameter :: C_ONE = dcmplx(1._dp, 0._dp)
  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  real (dp), parameter :: pi = 4._dp * atan(1._dp)
  real (dp), parameter :: tol = 1.0e-8
  character(len=*), parameter :: name_initial_state = "Neel"

  integer (ip)     ::  nspin, dim, n_disorder, steps, n_periods
  integer (ip)     ::  i, j, l, r, k, p, q, dim_Sz, Sz
  real (dp) :: norm

  real (dp), dimension(:), allocatable :: hz
  real (dp), dimension(:,:), allocatable :: Jxy, Vzz
  real (dp) :: T0, T1, Jxy_coupling, Vzz_coupling, hz_coupling, alpha
  complex (dcp) :: beta, delta

  real (dp), dimension(:), allocatable :: E
  real (dp), dimension(:,:), allocatable :: H, W_r, H2
  complex (dcp), allocatable :: psi(:), state(:), psi_swap(:), psi_Sz(:), psi_swap2(:)
  complex (dcp), dimension(:,:), allocatable :: U, U2, USwap, Utot

  integer(ip) :: count_beginning, count_end, count_rate

  real (dp), allocatable, dimension(:) :: sigmaz_previous, sigmaz_current, &
    & sigmaz_initial
  integer (ip), allocatable :: swap(:)

  real (dp), allocatable, dimension(:) :: sigmaz_previous2, sigmaz_current2, &
    & sigmaz_initial2

  real (dp), dimension(:), allocatable :: hz_sw
  real (dp), dimension(:,:), allocatable :: Jxy_sw, Vzz_sw

  logical :: SELECT
  EXTERNAL SELECT
  interface
    pure function ones(n) result(arr)
      use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
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
  dim = 2**nspin

  write (*,*) "Number of Disorder Realizations"
  read (*,*) n_disorder
  print*,""

  write (*,*) "Number of Periods per Step"
  read (*,*) n_periods
  print*,""
  !if(mod(n_periods,2) == 0) stop "Error: n_periods must be odd"

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  !write (*,*) "Number of Periods at each step"
  !read (*,*) n_periods
  !print*,""

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

  write (*,*) "Power-Law Coefficient Vzz = V_{ij} / |i-j|^alpha"
  read (*,*) alpha
  print*,""
 
  T1 = pi/4

  call system_clock(count_beginning, count_rate)

  !------------- Initial State ------------------!
  allocate(psi(dim))
  call buildstate(nspin, dim, psi)
  call printstate(nspin, dim, psi, trim(name_initial_state))
  call project(nspin, dim, psi, dim_Sz, Sz, psi_Sz)
  call printstate_Sz(nspin, dim_Sz, Sz, psi_Sz, trim(name_initial_state))
  allocate(sigmaz_initial(nspin))
  sigmaz_initial = sigmaz_Sz(nspin, dim_Sz, Sz, psi_Sz)
  deallocate(psi)

 
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
  print *, "USwap = "
  call print_unitary_Sz(nspin, dim_Sz, Sz, USwap)
  deallocate(E, W_r)

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jxy(nspin, nspin), Vzz(nspin,nspin), hz(nspin))
  allocate( Jxy_sw(nspin, nspin), Vzz_sw(nspin,nspin), hz_sw(nspin))

  !Allocate Floquet and MBL Operators
  allocate(H(dim_Sz,dim_Sz), H2(dim_Sz,dim_Sz), E(dim_Sz), W_r(dim_Sz,dim_Sz))
  allocate(U(dim_Sz,dim_Sz), U2(dim_Sz,dim_Sz), Utot(dim_Sz,dim_Sz))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(psi_swap(dim_Sz), psi_swap2(dim_Sz))

  !Allocate for Dynamics
  allocate( sigmaz_current(nspin) )
  allocate( sigmaz_current2(nspin) )
  allocate(swap(nspin))
  swap = (/ ( merge(i-1, i+1, mod(i,2) == 0 ) , i=1,nspin ) /)

  !$OMP PARALLEL
  call init_random_seed()
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "Current code segment is in parallel: ", omp_in_parallel()
  !$OMP do private(i, j, hz, Vzz, norm, &
  !$OMP & psi_swap, psi_swap2, H, H2, E, W_r, U, U2, Utot, &
  !$OMP & sigmaz_current, sigmaz_current2)
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
 
    !Jxy = -Jxy_coupling
    call random_number(Vzz)
    Vzz = -Vzz_coupling + Vzz_coupling*(Vzz - 0.5) !Vzz in [-V-V/2,-V+V/2]
    do k = 1, nspin-1
      Vzz(k,1:k) = 0
      Jxy(k,1:k) = 0
      do q = k+1, nspin
        Vzz(k,q) = Vzz(k,q) / ( abs(k-q)**alpha )
        Jxy(k,q) = -Jxy_coupling * merge(1, 0, q==k+1)
      enddo
    enddo
    Vzz(nspin,:) = 0
    Jxy(nspin,:) = 0
    Vzz = Vzz / norm_V(alpha, nspin)
 
    call random_number(hz)
    hz = 2*hz_coupling*(hz-0.5) !hz in [-hz_coupling, hz_coupling]

    call buildHMBL( nspin, dim_Sz, Sz, Jxy, Vzz, hz, H2 )
    H = H2
    print *, "H_original = "
    call print_hamiltonian_Sz(nspin, dim_Sz, Sz, H2)
    call diagSYM( 'V', dim_Sz, H2, E, W_r )
    call expSYM( dim_Sz, -C_UNIT*T0, E, W_r, U2 )
    U2 = matmul(USwap,U2)
    Utot = U2
    do r = 1, n_periods-1
      Utot = matmul(U2,Utot)
    enddo
    print *, "U_swap * H_original * U_swap = "
    call print_unitary_Sz(nspin, dim_Sz, Sz, matmul(matmul(USwap,H2),USwap))

    !Swapped parameters
    do k = 1, nspin
      do q = 1, nspin
        Vzz_sw(k,q) = Vzz(swap(k), swap(q)) + Vzz(swap(q), swap(k))
        Jxy_sw(k,q) = Jxy(swap(k), swap(q)) + Jxy(swap(q), swap(k))
      enddo
      Vzz_sw(k,1:k) = 0
      Jxy_sw(k,1:k) = 0
    enddo
    hz_sw = hz(swap)
    call buildHMBL( nspin, dim_Sz, Sz, Jxy_sw, Vzz_sw, hz_sw, H2 )
    print*, "H_swapped = "
    call print_hamiltonian_Sz(nspin, dim_Sz, Sz, H2)
    print *, "H_swapped - U_swap * H_original * U_swap = "
    call print_unitary_Sz(nspin, dim_Sz, Sz, H2 - matmul(matmul(USwap,H),USwap))
    print *, "[ H_swapped, H_original] = "
    call print_hamiltonian_Sz(nspin, dim_Sz, Sz, matmul(H2,H) - matmul(H,H2) )

  
    !Printing parameters
    print *, "Original Parameters: "
    do k = 1, nspin
      write (*,*) "Jxy = ", Jxy(k,:)
    enddo
    do k = 1, nspin
      write (*,*) "Vzz = ", Vzz(k,:)
    enddo
    write (*,*) "hz  = ", hz(:)
    print *, ""
    print *, "Swapped parameters: "
    do k = 1, nspin 
      write (*,*) "Jxy = ", Jxy_sw(k,:)
    enddo
    do k = 1, nspin
      write (*,*) "Vzz = ", Vzz_sw(k,:)
    enddo
    write (*,*) "hz  = ", hz_sw
    print *, ""

    !Skipping of n_periods at kick=0
    Jxy = ceiling(n_periods/2.0) * Jxy + floor(n_periods/2.0) * Jxy_sw
    Vzz = ceiling(n_periods/2.0) * Vzz + floor(n_periods/2.0) * Vzz_sw
    hz = ceiling(n_periods/2.0) * hz + floor(n_periods/2.0) * hz_sw

    print *, "ceil(L/2) = ", ceiling(n_periods/2.0), "floor(L/2) = ", floor(n_periods/2.0)
    print *, "Final Parameters: "
    do k = 1, nspin
      write (*,*) "Jxy = ", Jxy(k,:)
    enddo
    do k = 1, nspin
      write (*,*) "Vzz = ", Vzz(k,:)
    enddo
    write (*,*) "hz  = ", hz(:)
    print *, ""
 
 
    !---------------------------------------------------
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildHMBL( nspin, dim_Sz, Sz, Jxy, Vzz, hz, H )
    !print *, "HMBL = "
    !call print_hamiltonian_Sz(nspin, dim_Sz, Sz, H)
    call diagSYM( 'V', dim_Sz, H, E, W_r )
    call expSYM( dim_Sz, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    print *, "U_single   = "
    call print_unitary_Sz(nspin, dim_Sz, Sz, U2)
    print *, "U_multiple = "
    call print_unitary_Sz(nspin, dim_Sz, Sz, Utot)
    print *, "U_skipping = "
    call print_unitary_Sz(nspin, dim_Sz, Sz, U)
    print *, "U_multiple - U_skipping = "
    call print_unitary_Sz(nspin, dim_Sz, Sz, Utot-U)

    psi_swap = psi_Sz
    psi_swap2 = psi_Sz
    j=1
    sigmaz_current = sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap)
    sigmaz_current2 = sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap2)
    print *, "sigmaz_skipping = ", j, sigmaz_current
    print *, "sigmaz_multiple = ", j, sigmaz_current2
    print *, "sigmaz difference = ", j, sigmaz_current - sigmaz_current2
    print *, ""
    do j = 2, steps

      if(mod(j,10)==0) then
        norm = real(dot_product(psi_swap,psi_swap))
        psi_swap = psi_swap/sqrt(norm)
        norm = real(dot_product(psi_swap2,psi_swap2))
        psi_swap2 = psi_swap2/sqrt(norm)
      endif
      psi_swap = matmul(U, psi_swap)
      do r = 1, n_periods
        psi_swap2 = matmul(U2,psi_swap2)
      enddo
      sigmaz_current = sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap)
      sigmaz_current2 = sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap2)
      print *, "sigmaz_skipping = ", j, sigmaz_current
      print *, "sigmaz_multiple = ", j, sigmaz_current2
      print *, "sigmaz difference = ", j, sigmaz_current - sigmaz_current2
      print *, ""

    enddo
    !print *, ""

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end

pure function ones(n) result(arr)
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  integer (ip), intent(in) :: n
  real (dp) :: arr(n)

  arr = 1.0

end function
  


