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

  integer (c_int)     ::  nspin, dim, iteration, n_iterations, nspin_A, dim_A
  integer (c_int)     ::  i, j, k, p, l, dim_Sz0
  integer (c_int)     ::  unit_ovrlp

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r, LI, OVRLP, QE, IPR_arr
  complex(c_double_complex), dimension(:), allocatable :: PH, init_state, psi_Sz0
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap
  complex (c_double_complex) :: alpha, beta

  real (c_double), allocatable :: MI(:,:,:)
  integer (c_int) :: n_lim

  

  integer (c_int), allocatable :: idx(:), config(:), states(:)

  integer (c_int), allocatable :: idxSz0(:)
  integer (c_int) :: idx_state
  real (c_double) :: IMB_min, LI_min

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, time_min, count1, count2
  character(len=300) :: filestring, string, filestring_Neel, filestring_Nayak

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

  !---Coefficients for the reference state |psi> = \otimes_{k=1}^L P^{k-1}(alpha|0> + beta|1>)
  ! where P|0> = |1>, P|1> = |0>
  alpha = cos(pi/8)
  beta = sin(pi/8)

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
  
  write(filestring_Neel,93) "data/eigen/Sz0_DENSE_SWAP_hz_V_Disorder_Neel_Overlap_nspin", &
    & nspin, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_Nayak,93) "data/eigen/Sz0_DENSE_SWAP_hz_V_Disorder_Nayak_Overlap_nspin", &
    & nspin, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  open(newunit=unit_ovrlp,file=filestring_Neel)

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
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))
  allocate(U(dim_Sz0,dim_Sz0))
  !allocate(idx(dim_Sz0))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))

  !Allocate for Entanglement
  n_lim = min(nspin/2,2)
  allocate(LI(n_iterations,dim_Sz0)) 
  allocate(OVRLP(n_iterations,dim_Sz0)) 
  !allocate(MI(n_lim,n_iterations,dim_Sz0))
  allocate(QE(n_iterations,dim_Sz0))
  allocate(IPR_arr(n_iterations,dim_Sz0))

  !MI = 0
  LI = 0

  !Allocate reference state
  allocate(init_state(dim_Sz0))

  allocate(idxSz0(dim_Sz0))
  call zero_mag_states(nspin, dim_Sz0, idxSz0)

  alpha = cos(pi/8)
  beta = sin(pi/8)
  !call buildNayakState_Sz0(nspin, dim_Sz0, alpha, beta, init_state)
  call buildNeelState_Sz0(nspin, dim_Sz0, init_state)

  allocate(psi_Sz0(dim_Sz0))

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(iteration, h_z, Vint, H, E, W_r, &
  !$OMP & U, W, PH, i, psi_Sz0, nspin_A, dim_A)
  !!! <- Put this back after looking at single iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iteration = 1, n_iterations
    
    if (n_iterations < 10) then
     print *, "iteration = ", iteration
    else
      if (mod(iteration,n_iterations/10)==0) then 
        print *, "iteration = ", iteration
      endif
    endif

    !-------------------------------------------------
    !PARAMETERS
  
    Jint = -J_coupling
    !call random_number(Jint)
    !Jint = 2*J_coupling*(Jint - 0.5) !Jint in [-J,J]

    Vint = -V_coupling
    call random_number(Vint)
    Vint = -V_coupling + V_coupling*(Vint - 0.5) !Vint in [-3V/2,-V/2]

    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz, hz]

    !print *, "J = ", Jint(:)
    !print *, "Vint = ", Vint(:)
    !print *, "hz = ", h_z(:)
  
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    !print *, "UF diagonalized"
    E = real(C_UNIT*log(PH))
    
    !write(string,"(A12)") "MI"
    !print "(A,*(4X,A8))", repeat(trim(string),n_lim), "LI", "Overlap"
    do i = 1, dim_Sz0

      psi_Sz0 = W(1:dim_Sz0,i)

      QE(iteration,i) = E(i)
      LI(iteration,i) = local_imbalance_Sz0(nspin, dim_Sz0, psi_Sz0)
      OVRLP(iteration,i) = abs(dot_product(psi_Sz0,init_state))**2
      IPR_arr(iteration,i) = IPR(psi_Sz0)

      !do nspin_A = 1, n_lim
      !  dim_A = 2**nspin_A
      !  MI(nspin_A,iteration,i) = mutual_information_Sz0(nspin, nspin_A, dim, dim_Sz0, dim_A, psi_Sz0)
      !enddo

      !print "(*(4X,F8.5))", MI(:,iteration,i), LI(iteration,i), OVRLP(iteration,i) 
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  write(unit_ovrlp, *) "dim_Sz0 = ", dim_Sz0
  write(string,"(A26)") "MI"
  !write(unit_ovrlp, "(A12,A26,A,*(A24,2X))") "Iteration", "QE", repeat(trim(string),n_lim), "LI", "Overlap", "IPR"
  write(unit_ovrlp, "(A12,*(A24,2X))") "Iteration", "Overlap", "QE", "LI", "IPR"
  do iteration = 1, n_iterations
    do i = 1, dim_Sz0
      write(unit_ovrlp, *) iteration, OVRLP(iteration,i), QE(iteration,i), LI(iteration,i), IPR_arr(iteration,i)
      !write(unit_ovrlp, *) iteration, QE(iteration,i), MI(:,iteration,i), LI(iteration,i), OVRLP(iteration,i), IPR_arr(iteration,i)
      !if(OVRLP(iteration,i) > 0.05) then
      !  print *, iteration, QE(iteration,i), OVRLP(iteration,i)
      !endif
    enddo
  enddo

  

  deallocate(Jint, Vint, h_z)
  deallocate(H,E,W_r)
  deallocate(USwap)
  deallocate(U, PH, W)
  deallocate(QE, IPR_arr, LI, OVRLP)

  call take_time(count_rate, count_beginning, count1, 'T', "Program")

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end