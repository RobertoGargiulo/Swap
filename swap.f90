program swap

  use exponentiate
  use genmat
  use printing
  use iso_c_binding
  !use general
  implicit none

  complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, steps, n_iterations
  integer (c_int)     ::  i, j, k, p
  integer (c_int)     ::  unit_mag, unit_ph, unit_w
  integer (c_int), dimension(:), allocatable  :: base_state, config

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_x, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, h_coupling, hz_coupling, kick 
  
  real (c_double) :: mag, norm
  complex (c_double_complex) :: alpha, beta

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r

  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:), allocatable :: state, init_state
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap, U_MBL

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, day, month, year, date(8), time_min
  real (c_double) :: time_s
  character(len=100) :: filestring
  character(len=8) :: time_string


  !Parametri Modello: J, h_x, h_z, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine, Steps di Evoluzione, Stato Iniziale

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""

  !Standard Values
  T0 = 1
  J_coupling = 1
  V_coupling = J_coupling
  hz_coupling = J_coupling
  h_coupling = 0.3
  kick = 0.1
  T1 = pi/4 + kick
  
  write (*,*) "Longitudinal Interaction Constant J * ZZ"
  read (*,*) J_coupling
  print*,""

  write (*,*) "Transverse Interaction Constant V * (XX + YY)"
  read (*,*) V_coupling
  print*,""

  write (*,*) "Transverse Field h_x * X"
  read (*,*) h_coupling
  print*,""

  write (*,*) "Longitudinal Field h_z * Z"
  read (*,*) hz_coupling
  print*,""
  !---Read below for distribution of J, V, hx, hz
  
  write (*,*) "Perturbation on Kick, epsilon = T1 - pi/4"
  read (*,*) kick
  print*,""

  T1 = pi/4 + kick

    !Coefficienti dello stato iniziale |psi> = (alpha|up>+beta|down>)^L
  alpha = 1
  beta = 0

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !Data Files
  write(filestring,92) "data/magnetizations/Swap_Sz_nspin", nspin, "_steps", steps, &
  &  "_iterations", n_iterations, "_J", J_coupling, "_h", h_coupling, "_kick", kick, ".txt"
  !open(newunit=unit_mag,file=filestring)

  write(filestring,92) "data/eigenvalues/Swap_PH_nspin", nspin, "_steps", steps, &
  &  "_iterations", n_iterations, "_J", J_coupling, "_h", h_coupling, "_kick", kick, ".txt"
  !open(newunit=unit_ph, file=filestring)
  
  write(filestring,92) "data/eigenvalues/Swap_W_nspin", nspin, "_steps", steps, &
  &  "_iterations", n_iterations, "_J", J_coupling, "_h", h_coupling, "_kick", kick, ".txt"
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A,F4.2, A)

  !open(newunit=unit_w, file=filestring)
 
  !------------------------------------------------

  !BUILD INITIAL STATE (of type staggered)
  call buildStaggState(nspin, dim, alpha, beta, init_state)

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim,dim), E(dim), W_r(dim,dim), USwap(dim,dim))
  call buildHSwap(nspin, dim, H)
  call diagSYM( 'V', dim, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim, -C_UNIT*T1, E, W_r, USwap) 
  deallocate(E, W_r)
!  print *, "USwap = "
!  call printmat(dim, USwap, 'R')

  !Allocate local interactions and fields
  allocate( Jint(nspin-1), Vint(nspin-1), h_x(nspin), h_z(nspin))

  !Allocate Floquet and MBL Operators
  allocate(U(dim,dim), H_MBL(dim,dim), U_MBL(dim,dim), E(dim), W_r(dim,dim))

  !Allocate initial and generic state
  allocate(state(dim), init_state(dim))

  !Allocate for Eigenvalues/Eigenvectors
  !allocate(PH(dim), W(dim,dim))


  do iteration = 1, n_iterations
    
    if (mod(iteration,10)==0) then 
      print *, "iteration = ", iteration
    endif

    !-------------------------------------------------
    !PARAMETERS
  
    call random_number(Jint)
    Jint = 2*J_coupling*(Jint - 0.5) !Jint in [-J,J]
    !Jint = 2*pi

    call random_number(Vint)
    Vint = 2*V_coupling*(Vint - 0.5) !Jint in [-J,J]
    !Vint = Jint
  
    call random_number(h_x)
    h_x = 2*h_coupling*(h_x - 0.5) !h_x in [-h_coupling, h_coupling]
  
    call random_number(h_z)
    h_z = 2*hz_coupling(h_z-0.5) !h_z in [-hz_coupling, hz_coupling]
  
    !write (*,*) "Jint = ", Jint(:)
    !write (*,*) "h_x = ", h_x(:)
    !write (*,*) "h_z = ", h_z(:)
    !print *, ""
  
    !---------------------------------------------------
  
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildHMBL( nspin, dim, Jint, Vint, h_z, h_z, H_MBL )
    call diagSYM( 'V', dim, H_MBL, E, W_r )
    call expSYM( dim, -C_UNIT*T0, E, W_r, U_MBL  )

    U = matmul(U_MBL, USwap)
  
!    print *, "UF = "
!    call printmat(dim, U,'C')
    !-------------------------------------------------

    !DIAGONALIZE FLOQUET OPERATOR
!    call diagUN( SELECT, dim, U, PH, W)
  
!    print *, "Eigenvalues of U_F "
!    print "(*(/f15.10spf15.10' i'))", PH(:)
!    print *,""
!  
!    print *, "Eigenvectors of U_F "
!    do i = 1,dim
!      write (*,"(*('|',f5.2spf5.2' i'))") W(i,:)
!    enddo
!    print *,""
    !------------------------------------------------
  
    !PRINT Eigenvalues/Eigenvectors to file
    !call writevec(unit_ph,dim,PH,'C')
    !call writemat(unit_w,dim,W,'C')
  
    !-----------------------------------------------
  
    !EVOLUTION OF INITIAL STATE and COMPUTATION OF MAGNETIZATION 
 
    state = init_state
    norm = dot_product(state,state)
    !call printvec(dim, state, 'R')
    j = 1
    
    !write(unit_mag,*) "iteration = ", iteration
    !write(unit_mag,*) mag_stag_z(nspin, dim, state), j
    
    print *, mag_stag_z(nspin, dim, state), j, norm
  
    do j = 2, steps
      state = matmul(U,state)
      norm = dot_product(state,state)
      state = state / sqrt(norm)
      !write(unit_mag,*) mag_stag_z(nspin, dim, state), j
      print *, mag_stag_z(nspin, dim, state), j, norm
    enddo
    print *, ""
 

  enddo 
  deallocate(state, U, PH, W, Jint, h_x, h_z)

  close(unit_mag)
  close(unit_ph)
  close(unit_w)

  call system_clock(count_end)

  time_s = real(count_end - count_beginning) / real(count_rate)
  time_min = int(time_s/60)

  print "(A,1X,I4,A,2X,F15.10,A)", "Elapsed Time: ", time_min, "min", time_s - 60*time_min, "sec"
  print *, ""

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
