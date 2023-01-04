program swap

  use exponentiate
  use genmat
  use printing
  use MBL
  use omp_lib
  !use sorts
  use iso_c_binding
  !use general
  implicit none

  !complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  !complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, steps, n_iterations
  integer (c_int)     ::  j, l, i, dim_Sz0
  integer (c_int)     ::  unit_avg, unit_zmag
  integer (c_int)     ::  start

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 

  real(c_double) :: imb_t_avg, imb_t_sigma, imb_t_avg2, imb_t_sigma2
 
  real (c_double) :: norm
  real (c_double), dimension(:), allocatable :: imb_avg, imb_sq!, imb_avg2, imb_sq2
  !real (c_double), dimension(:), allocatable :: li_avg, li_sq, lo_avg, lo_sq
  real (c_double), allocatable :: zmag(:), zmag_avg(:,:), zmag_sq(:,:)

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap, UF

  complex(c_double_complex), dimension(:), allocatable :: init_state, state!, state2

  integer(c_int) :: count_beginning, count_end, count_rate!, time_min, count1, count2
  character(len=200) :: filestring, filestring_Neel, filestring_Nayak, &
    & filestring_zmag, filestring_imb_LI1, filestring_zmag_I0, filestring_imb_I0LI1
  

  real (c_double) :: IMB, LI, theta, fr_theta
  integer (c_int) :: idx_state
  integer (c_int), allocatable :: idxSz0(:), config(:)
  complex (c_double_complex) :: alpha, beta

  !Parametri Modello: J, V, h_z, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine, Steps di Evoluzione, Stato Iniziale

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin
  dim_Sz0 = binom(nspin,nspin/2)

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
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

  write (*,*) "Fraction fr_theta for angle theta = fr_theta * pi for alpha = cos(theta), beta = sin(theta)"
  read (*,*) fr_theta
  print *, ""

  !---Read below for distributions of J, V, hz
  
  T1 = pi/4 + kick


  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES



  write(filestring_Neel,93) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_Neel_AVG_FLUCT_Imbalance_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_Nayak,93) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_Nayak_AVG_FLUCT_Imbalance_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_imb_LI1,94) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_LI1ProdState_frtheta", &
    & fr_theta, "_AVG_FLUCT_Imbalance_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_zmag,94) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_LI1ProdState_frtheta", &
    & fr_theta, "_AVG_FLUCT_ZMag_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_imb_I0LI1,94) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_I0LI1ProdState_frtheta", &
    & fr_theta, "_AVG_FLUCT_Imbalance_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  write(filestring_zmag_I0,94) "data/magnetizations/Sz0_DENSE_SWAP_hz_V_Disorder_I0LI1ProdState_frtheta", &
    & fr_theta, "_AVG_FLUCT_ZMag_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_iterations", n_iterations, &
    & "_J", J_coupling, "_V", int(V_coupling), V_coupling-int(V_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, ".txt"

  open(newunit=unit_zmag,file=filestring_zmag_I0)
  open(newunit=unit_avg,file=filestring_imb_I0LI1)

  !91  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  !92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)
  93  format(A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)
  94  format(A,F4.2, A,I0, A,I0, A,F4.2, A,I0, A,F4.2, A,I0,F0.2, A,I0,F0.2, A,F4.2, A)
!
 
  !------------------------------------------------

  !BUILD INITIAL STATE (of type staggered)
  allocate(init_state(dim_Sz0))

  theta = pi * fr_theta
  alpha = cos(theta)
  beta = sin(theta)
  !call buildNayakState_Sz0(nspin, dim_Sz0, alpha, beta, init_state)
  !call buildNeelState_Sz0(nspin, dim_Sz0, init_state)
  call buildI0LI1ProdState_Sz0(nspin, dim_Sz0, alpha, beta, init_state)
  !call buildLI1ProdState_Sz0(nspin, dim_Sz0, alpha, beta, init_state)

  !call printvec(dim_Sz0, init_state, 'A')

  allocate(idxSz0(dim_Sz0), config(nspin))
  call zero_mag_states(nspin, dim_Sz0, idxSz0)
  do i = 1, dim_Sz0
    if (abs(init_state(i))**2 > 1.0e-6) then
      l = idxSz0(i)
      call decode(l, nspin, config)
      print "( 4X,F8.4, 4X,I6, 4X,*(I0) )", abs(init_state(i))**2, l, config(:)
    endif
  enddo


  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  call buildSz0_HSwap(nspin, dim_Sz0, H)
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim_Sz0, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)
  
  !print *, "USwap = "
  !call printmat(dim_Sz0, USwap, 'C')
  !print *, "USwap*USwap^dagger = "
  !call printmat(dim_Sz0, matmul(USwap,USwap), 'C')

  !---------------------------------------------------
  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_z(nspin))

  !Allocate Floquet and MBL Operators
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(U(dim_Sz0,dim_Sz0), H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))
  !allocate(UF(dim_Sz0,dim_Sz0))

  !Allocate initial and generic state
  allocate( state(dim_Sz0) )
  !allocate( state2(dim_Sz0) )

  !Allocate observables and averages
  allocate( imb_avg(steps), imb_sq(steps))
  !allocate( imb_avg2(steps), imb_sq2(steps))
  allocate( zmag(nspin), zmag_avg(steps,nspin), zmag_sq(steps,nspin))
  !allocate( li_avg(steps), li_sq(steps))
  !allocate( lo_avg(steps), lo_sq(steps))

  !Allocate for Eigenvalues/Eigenvectors
  !allocate(PH(dim_Sz0))
  !allocate(W(dim_Sz0,dim_Sz0))

  imb_avg = 0
  imb_sq = 0
  !imb_avg2 = 0
  !imb_sq2 = 0
  zmag_avg = 0
  zmag_sq = 0
  !li_avg = 0
  !li_sq = 0
  !lo_avg = 0
  !lo_sq = 0
  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+:imb_avg, imb_sq, zmag_avg, zmag_sq) & 
  !$OMP & private(iteration, h_z, Vint, j, norm, &
  !$OMP & state, H, E, W_r, U, zmag )
  !, UF, state2), imb_avg2, imb_sq2), li_avg, li_sq, lo_avg, lo_sq
  !PH, W)
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

    Vint = -V_coupling
    call random_number(Vint)
    Vint = -V_coupling + V_coupling*(Vint - 0.5) !Jint in [-V,V]
 
    call random_number(h_z)
    h_z = 2*hz_coupling*(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]
 
!    write (*,*) "Jint = ", Jint(:)
!    write (*,*) "Vint = ", Vint(:)
!    write (*,*) "h_z = ", h_z(:)
!    print *, ""
 
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    !print *, "Building Dense Operator"
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)

    j = 1
    state = init_state
    norm = real(dot_product(state,state))
    state = state/sqrt(norm)
    imb_avg(j) = imb_avg(j) + imbalance_Sz0(nspin, dim_Sz0, state)
    imb_sq(j) = imb_sq(j) + imbalance_Sz0(nspin, dim_Sz0, state)**2
    call local_zmag_Sz0(nspin, dim_Sz0, state, zmag )
    zmag_avg(j,:) = zmag_avg(j,:) + zmag(:)
    zmag_sq(j,:) = zmag_sq(j,:) + zmag(:)**2
    !li_avg(j) = li_avg(j) + local_imbalance_Sz0(nspin, dim_Sz0, state)
    !li_sq(j) = li_sq(j) + local_imbalance_Sz0(nspin, dim_Sz0, state)**2
    !lo_avg(j) = lo_avg(j) + local_overlap_Sz0(nspin, dim_Sz0, state)
    !lo_sq(j) = lo_sq(j) + local_overlap_Sz0(nspin, dim_Sz0, state)**2


    !state2 = init_state
    !norm = real(dot_product(state2,state2))
    !state2 = state2/sqrt(norm)
    !imb_avg2(j) = imb_avg2(j) + imbalance_Sz0(nspin, dim_Sz0, state2)
    !imb_sq2(j) = imb_sq2(j) + imbalance_Sz0(nspin, dim_Sz0, state2)**2

    do j = 2, steps
      state = matmul(U,state) !U -> UF
      norm = real(dot_product(state,state))
      state = state / sqrt(norm)
      imb_avg(j) = imb_avg(j) + imbalance_Sz0(nspin, dim_Sz0, state) 
      imb_sq(j) = imb_sq(j) + imbalance_Sz0(nspin, dim_Sz0, state)**2
      call local_zmag_Sz0(nspin, dim_Sz0, state, zmag )
      zmag_avg(j,:) = zmag_avg(j,:) + zmag(:)
      zmag_sq(j,:) = zmag_sq(j,:) + zmag(:)**2
      !li_avg(j) = li_avg(j) + local_imbalance_Sz0(nspin, dim_Sz0, state)
      !li_sq(j) = li_sq(j) + local_imbalance_Sz0(nspin, dim_Sz0, state)**2
      !lo_avg(j) = lo_avg(j) + local_overlap_Sz0(nspin, dim_Sz0, state)
      !lo_sq(j) = lo_sq(j) + local_overlap_Sz0(nspin, dim_Sz0, state)**2

      !state2 = matmul(U,state2)
      !norm = real(dot_product(state2,state2))
      !state2 = state2 / sqrt(norm)
      !imb_avg2(j) = imb_avg2(j) + imbalance_Sz0(nspin, dim_Sz0, state2)
      !imb_sq2(j) = imb_sq2(j) + imbalance_Sz0(nspin, dim_Sz0, state2)**2
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 


  imb_avg = imb_avg/n_iterations
  imb_sq = sqrt( (imb_sq/n_iterations - imb_avg**2) / n_iterations )
  zmag_avg = zmag_avg/n_iterations
  zmag_sq = sqrt( (zmag_sq/n_iterations - zmag_avg**2) / n_iterations )
  !li_avg = li_avg/n_iterations
  !li_sq = sqrt( (li_sq/n_iterations - li_avg**2) / n_iterations )
  !lo_avg = lo_avg/n_iterations
  !lo_sq = sqrt( (lo_sq/n_iterations - lo_avg**2) / n_iterations )
  !imb_avg2 = imb_avg2/n_iterations
  !imb_sq2 = sqrt( (imb_sq2/n_iterations - imb_avg2**2) / n_iterations )

  do j = 1, steps
    write(unit_avg,*) j*T0, imb_avg(j), imb_sq(j)!, li_avg(j), li_sq(j), lo_avg(j), lo_sq(j)!, imb_avg2(j), imb_sq2(j)
    !print *, j*T0, imb_avg(j), imb_sq(j), imb_avg2(j), imb_sq2(j)
  enddo
  
  do j = 1, steps
    write(unit_zmag,*) j*T0, zmag_avg(j,:), zmag_sq(j,:)
  enddo

  start = int(100/T0) !The average starts from the step for which 100 = start*T0
  call time_avg('T', steps, start, imb_avg, imb_sq, imb_t_avg, imb_t_sigma)
  !call time_avg('F', steps, start, imb_avg2, imb_sq2, imb_t_avg2, imb_t_sigma2)


  !write(unit_avg,*) "Imbalance Time Averages and Errors of Imbalance"
  !write(unit_avg,*) imb_t_avg, imb_t_sigma, imb_t_avg2, imb_t_sigma2

  print *,"Imbalance Time Averages and Errors"
  print *, imb_t_avg, imb_t_sigma!, imb_t_avg2, imb_t_sigma2

  deallocate(Jint, Vint, h_z)
  deallocate(H,E,W_r,U)
  deallocate(USwap)
  deallocate(state, imb_avg, imb_sq)
  !deallocate(state2, imb_avg2, imb_sq2)
  !deallocate(USwap, UF)

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

  !close(unit_avg)

  !call take_time(count_rate, count_beginning, count_end)

end program swap
