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
  integer (c_int)     ::  i, j, k, p, l, dim_Sz0
  integer (c_int)     ::  unit_ent

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 

  real (c_double), dimension(:), allocatable :: E, QE
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH, psi_Sz0
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap

  real (c_double), allocatable :: MI(:,:,:), MI_avg(:,:), MI_dis_avg(:)
  integer (c_int) :: nspin_A, dim_A, n_lim
  real (c_double), allocatable, dimension(:,:) :: MBE, IPR_arr, LI, IMBsq, CORR_Z

  integer (c_int), allocatable :: idx(:), config(:), states(:)

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, time_min, count1, count2
  character(len=300) :: filestring, string

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
  
  write(filestring,93) "data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Entanglement_nspin", &
    & nspin, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"
  open(newunit=unit_ent,file=filestring)

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
  allocate(idx(dim_Sz0))
  allocate(QE(dim_Sz0))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))

  !Allocate for Entanglement
  n_lim = min(nspin/2,2)
  allocate(MI(n_lim,n_iterations,dim_Sz0), MI_avg(n_lim,n_iterations), MI_dis_avg(n_lim))
  allocate(LI(n_iterations,dim_Sz0), IMBsq(n_iterations,dim_Sz0), MBE(n_iterations,dim_Sz0), IPR_arr(n_iterations,dim_Sz0) )
  allocate(CORR_Z(n_iterations,dim_Sz0))


  allocate(config(nspin), states(dim_Sz0))
  call zero_mag_states(nspin, dim_Sz0, states)

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(iteration, h_z, H, W_r, &
  !$OMP & U, W, PH, nspin_A, dim_A, i, psi_Sz0) 
  !E <- Put this back after looking at single iterations !!!!!!!!!
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
  
    !call random_number(Jint)
    !Jint = 2*J_coupling*(Jint - 0.5) !Jint in [-J,J]
    Jint = -J_coupling

    call random_number(Vint)
    Vint = V_coupling + V_coupling*(Vint - 0.5) !Jint in [-V,V]
    !Vint = -V_coupling
  
    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]

    !print *, "J = ", Jint(:)
    !print *, "Vint = ", Vint(:)
    !print *, "hz = ", h_z(:)
  
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    !do i = 1, dim_Sz0
    !  print *, E(i)
    !enddo
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    E = real(C_UNIT*log(PH))
    !call dpquicksort(E)
    !print *, E(:)
    !call exact_quasi_energies_Sz0(nspin, dim_Sz0, V_coupling, h_z, QE)

    !do i = 1, dim_Sz0
    !  print *, E(i), QE(i), C_UNIT*log(exp(-C_UNIT*QE(i)))
    !enddo

    
    !write(string,"(A12)") "MI"
    !print "(A,*(4X,A8))", repeat(trim(string),n_lim), "LI", "IPR", "MBE", "I^2", "CORR", "E"
    do i = 1, dim_Sz0

      psi_Sz0 = W(1:dim_Sz0,i)

      !print *, "psi = "
      !call printvec(dim_Sz0, psi_Sz0, 'A')
      !print *, "   abs(c_i)^2     I^2     IPR       l    config"
      !do k = 1, dim_Sz0
      !    l = states(k)
      !    call decode(l, nspin, config)
      !    print 1, abs(psi_Sz0(k))**2, imbalance_sq_basis(nspin, l), IPR(psi_Sz0), l, config(:)
      !enddo
      !1 format (10X,F4.2, 4X,F5.2, 4X,F4.2, 4X,I3, 4X,*(I0))

      IMBsq(iteration,i) = imbalance_sq_Sz0(nspin,dim_Sz0,psi_Sz0)
      LI(iteration,i) = local_imbalance_Sz0(nspin,dim_Sz0,psi_Sz0)
      IPR_arr(iteration,i) = IPR(psi_Sz0)
      MBE(iteration,i) = max_bipartite_entropy_Sz0(nspin, dim, dim_Sz0, psi_Sz0)
      CORR_Z(iteration,i) = sigmaz_corr_c_Sz0(nspin, dim_Sz0, 1, nspin, psi_Sz0)

      do nspin_A = 1, n_lim
        dim_A = 2**nspin_A
        MI(nspin_A,iteration,i) = mutual_information_Sz0(nspin, nspin_A, dim, dim_Sz0, dim_A, psi_Sz0)
      enddo
      !print "(*(4X,F8.4) )", MI(:,iteration,i), LI(iteration,i), &
      !  & IPR_arr(iteration,i), MBE(iteration,i), IMBsq(iteration,i), &
      !  & CORR_Z(iteration,i), E(i)
      !print *, ""

      !do k = 1, nspin
        !do j = k, nspin
          !print *, k, j, sigmaz_corr_c_Sz0(nspin, dim_Sz0, k, j, psi_Sz0)
       ! enddo
      !enddo

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  !print *, E(:)
  write(string,"(A12)") "MI"
  write(unit_ent, "(A,*(4X,A8))") repeat(trim(string),n_lim), "LI", "IPR", "MBE", "I^2", "E"
  do iteration = 1, n_iterations
    do i = 1, dim_Sz0
      write(unit_ent, "(*(4X,F8.4) )") MI(:,iteration,i), LI(iteration,i), &
        & IPR_arr(iteration,i), MBE(iteration,i), IMBsq(iteration,i), & 
        & CORR_Z(iteration,i), E(i)
    enddo
  enddo

  MI_avg = 0
  MI_dis_avg = 0
  !print *, "Mutual Information (MI, nspin_A, iteration)"
  do nspin_A = 1, n_lim
    do iteration = 1, n_iterations
      do i = 1, dim_Sz0
        MI_avg(nspin_A,iteration) = MI_avg(nspin_A,iteration) + MI(nspin_A,iteration,i)
      enddo
      MI_avg = MI_avg/dim_Sz0
      !print *, MI_avg(nspin_A,iteration), nspin_A, iteration
      MI_dis_avg(nspin_A) = MI_dis_avg(nspin_A) + MI_avg(nspin_A,iteration)
    enddo
    MI_dis_avg = MI_dis_avg/n_iterations
    !print *, MI_dis_avg(nspin_A), nspin_A
    !print *, ""
  enddo

  deallocate(Jint, Vint, h_z)
  deallocate(H,E,W_r)
  deallocate(USwap)
  deallocate(U, PH, W)

  call take_time(count_rate, count_beginning, count1, 'T', "Program")

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
