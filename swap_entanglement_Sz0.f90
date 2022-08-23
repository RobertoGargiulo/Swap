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

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH, psi_Sz0
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap

  real (c_double), allocatable :: MI(:,:,:), MI_avg(:,:), MI_dis_avg(:)
  integer (c_int) :: nspin_A, dim_A, n_lim
  real (c_double) :: MBE

  integer (c_int), allocatable :: idx(:), config(:), states(:)

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, time_min, count1, count2
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
  
  write(filestring,93) "data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Entanglement_nspin", &
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
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))
  allocate(U(dim_Sz0,dim_Sz0))
  allocate(idx(dim_Sz0))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))

  !Allocate for Entanglement
  allocate(MI(nspin/2,n_iterations,dim_Sz0), MI_avg(nspin/2,n_iterations), MI_dis_avg(nspin/2))


  allocate(config(nspin), states(dim_Sz0))
  call zero_mag_states(nspin, dim_Sz0, states)

  n_lim = min(nspin/2,1)

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(iteration, h_z, H, E, W_r, &
  !$OMP & U, W, PH, nspin_A, dim_A, i)
  do iteration = 1, n_iterations
    
    if (mod(iteration,10)==0) then 
      print *, "iteration = ", iteration
    endif

    !-------------------------------------------------
    !PARAMETERS
  
    !call random_number(Jint)
    !Jint = 2*J_coupling*(Jint - 0.5) !Jint in [-J,J]
    Jint = -J_coupling

    !call random_number(Vint)
    !Vint = 2*V_coupling*(Vint - 0.5) !Jint in [-V,V]
    Vint = -V_coupling
  
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
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    E = real(C_UNIT*log(PH))
    !call sort_index(E, idx)
    call dpquicksort(E)
    !print *, E(:)

    print *, "       I       I^2        LI       IPR       MBE        MI"
    do i = 1, dim_Sz0
      psi_Sz0 = W(1:dim_Sz0,i)
      !print *, "psi = "
      !call printvec(dim_Sz0, psi_Sz0, 'A')
      MBE = max_bipartite_entropy_Sz0(nspin, dim, dim_Sz0, psi_Sz0)
      !if ( abs(imbalance_sq_Sz0(nspin, dim_Sz0, psi_Sz0)) > 1.0e-10 ) then
      !  print *, "          I             I^2                     LI                      IPR                MBE"
      !  print *, imbalance_Sz0(nspin, dim_Sz0, psi_Sz0), imbalance_sq_Sz0(nspin, dim_Sz0, psi_Sz0), &
      !  & local_imbalance(nspin, dim, psi_Sz0), IPR(psi_Sz0), MBE
      !endif

      !print *, "   abs(c_i)^2     I^2     IPR       l    config"
      do k = 1, dim_Sz0
          l = states(k)
          call decode(l, nspin, config)
          !print 1, abs(psi_Sz0(k))**2, imbalance_sq_basis(nspin, l), IPR(psi_Sz0), l, config(:)
      enddo

      1 format (10X,F4.2, 4X,F5.2, 4X,F4.2, 4X,I3, 4X,*(I0))

      do nspin_A = 1, n_lim
        !print *, "   nspin_A  eigenvector        MI"
        dim_A = 2**nspin_A
        MI(nspin_A,iteration,i) = mutual_information_Sz0(nspin, nspin_A, dim, dim_Sz0, dim_A, W(1:dim_Sz0,i))
        !print *, "       I       I^2        LI       IPR       MBE        MI"
        print "(*(4X,F15.7) )", imbalance_Sz0(nspin, dim_Sz0, psi_Sz0), imbalance_sq_Sz0(nspin, dim_Sz0, psi_Sz0), &
        & local_imbalance(nspin, dim, psi_Sz0), IPR(psi_Sz0), MBE, MI(nspin_A,iteration,i)
        !print *, nspin_A, i, MI(nspin_A,iteration,i)
      enddo
      !print *, ""
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  !print *, E(:)
  !do j = 1, dim_Sz0
  !  write(unit_ph,*) E(j), pair_gaps(j), near_gaps(j), log_gap(j)
  !enddo

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
