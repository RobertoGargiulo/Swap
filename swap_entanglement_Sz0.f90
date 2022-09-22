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

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH, psi_Sz0, init_state
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap
  complex (c_double_complex) :: alpha, beta

  real (c_double), allocatable :: MI(:,:,:)
  integer (c_int) :: nspin_A, dim_A, n_lim, n_min, n_max
  real (c_double), allocatable, dimension(:,:) :: LI, QE, CEE, IPR_arr, CORR_Z, MBE, IMBsq
  real (c_double) :: IMB_min, LI_min

  integer (c_int), allocatable :: idx(:), config(:), states(:), idx_Sz0(:)
  integer (c_int) :: idx_state

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate
  character(len=300) :: filestring, filestring_disV, string

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
  
  write(filestring,93) "data/eigen/Sz0_DENSE_SWAP_hz_Disorder_Entanglement_nspin", &
    & nspin, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_disV,93) "data/eigen/Sz0_DENSE_SWAP_hz_V_Disorder_Entanglement_nspin", &
    & nspin, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  open(newunit=unit_ent,file=filestring_disV)

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
  !allocate(QE(dim_Sz0), QE_new(dim_Sz0), QE_old(dim_Sz0), Es(dim_Sz0), QE_alt(dim_Sz0))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))

  !Allocate for Entanglement
  n_max = min(nspin/2,4)
  n_min = min(nspin/2,2)
  n_lim = n_max - n_min + 1
  allocate(MI(n_lim,n_iterations,dim_Sz0))
  allocate(LI(n_iterations,dim_Sz0),IPR_arr(n_iterations,dim_Sz0))
  allocate(CEE(n_iterations,dim_Sz0))
  !allocate( IMBsq(n_iterations,dim_Sz0), MBE(n_iterations,dim_Sz0) )
  !allocate(CORR_Z(n_iterations,dim_Sz0))

  MI = 0
  IPR_arr = 0
  CEE = 0
  !CORR_Z = 0
  !IMBsq = 0
  !MBE = 0


  allocate(config(nspin), states(dim_Sz0))
  call zero_mag_states(nspin, dim_Sz0, states)

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(iteration, h_z, Vint, H, E, W_r, &
  !$OMP & U, W, PH, nspin_A, dim_A, i, k, psi_Sz0)
  !E <- Put this back after looking at single iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    Vint = -V_coupling + V_coupling*(Vint - 0.5) !Vint in [V/2,3V/2]

    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]

    !print *, "J = ", Jint(:)
    !print *, "Vint = ", Vint(:)
    !print *, "hz = ", h_z(:)
  
    !---------------------------------------------------
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    !E = real(C_UNIT*log(PH))
    !call dpquicksort(E)

    !do i = 1, dim_Sz0
    !  print *, E(i), QE(i)
    !enddo

    
    write(string,"(A12)") "MI"
    !print "(A,*(4X,A8))", repeat(trim(string),n_lim), "LI", "IPR", "CORR"
    do i = 1, dim_Sz0

      psi_Sz0 = W(1:dim_Sz0,i)

      LI(iteration,i) = local_imbalance_Sz0(nspin, dim_Sz0, psi_Sz0)
      IPR_arr(iteration,i) = IPR(psi_Sz0)
      CEE(iteration,i) = comb_entropy_Sz0(nspin, dim, dim_Sz0, psi_Sz0)

      !do nspin_A = n_min, n_max
      !  dim_A = 2**nspin_A
      !  MI(nspin_A,iteration,i) = mutual_information_Sz0(nspin, nspin_A, dim, dim_Sz0, dim_A, psi_Sz0)
      !enddo


      !do k = 1, nspin
        !do j = k, nspin
          !print *, k, j, sigmaz_corr_c_Sz0(nspin, dim_Sz0, k, j, psi_Sz0)
       ! enddo
      !enddo
      !CORR_Z(iteration,i) = sigmaz_corr_c_Sz0(nspin, dim_Sz0, 1, nspin, psi_Sz0)


      !!---- Print to output -------!!!!
      !print *, "     abs(c_i)^2     l    config"
      !do k = 1, dim_Sz0
      !  if( abs(psi_Sz0(k))**2 > 1.0e-4) then
      !    l = states(k)
      !    call decode(l, nspin, config)
      !    print 1, abs(psi_Sz0(k))**2, l, config(:)
      !    !print 2, E(i), exact_energy(nspin, Vint, h_z, l), &
      !    !  & exact_quasi_energy(nspin, Vint, h_z, l), &
      !    !  & psi_Sz0(k), &
      !    !  & l, config(:)
      !  endif
      !enddo
      !1 format (4X,F6.3, 4X,I3, 4X,*(I0))
      !2 format (4X,F6.3, 4X,F8.3, 4X,F8.3, 4X, 1(F8.3,F8.3X'i':X), 4X,I3, 4X,*(I0))

      !write(string,"(A26)") "MI"
      !print "(A12,A,*(A24,2X))", "Iteration" , repeat(trim(string),n_lim), "CEE", "LI", "IPR"!, "MBE", "I^2", "CORR_Z", "E"
      !print *, iteration, MI(:,iteration,i), CEE(iteration,i), &
      !  &  LI(iteration,i), IPR_arr(iteration,i)!, CORR_Z(iteration,i), E(i)
      !---------------------------------------------------------------------
      !print*, ""

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  !print *, E(:)
  write(string,"(A26)") "MI"
  write(unit_ent, "(A12,A,*(A24,2X))") "Iteration" , repeat(trim(string),n_lim), "CEE", "LI", "IPR"!, "MBE", "I^2", "CORR_Z", "E"
  do iteration = 1, n_iterations
    do i = 1, dim_Sz0
      write(unit_ent, *) iteration, MI(:,iteration,i), CEE(iteration,i), &
        &  LI(iteration,i), IPR_arr(iteration,i)!, CORR_Z(iteration,i), E(i)
    enddo
  enddo

  !MI_avg = 0
  !MI_dis_avg = 0
  !!print *, "Mutual Information (MI, nspin_A, iteration)"
  !do nspin_A = 1, n_lim
  !  do iteration = 1, n_iterations
  !    do i = 1, dim_Sz0
  !      MI_avg(nspin_A,iteration) = MI_avg(nspin_A,iteration) + MI(nspin_A,iteration,i)
  !    enddo
  !    MI_avg = MI_avg/dim_Sz0
  !    !print *, MI_avg(nspin_A,iteration), nspin_A, iteration
  !    MI_dis_avg(nspin_A) = MI_dis_avg(nspin_A) + MI_avg(nspin_A,iteration)
  !  enddo
  !  MI_dis_avg = MI_dis_avg/n_iterations
  !  !print *, MI_dis_avg(nspin_A), nspin_A
  !  !print *, ""
  !enddo

  deallocate(Jint, Vint, h_z)
  deallocate(H,E,W_r)
  deallocate(USwap)
  deallocate(U, PH, W)

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
