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

  integer (c_int)     ::  nspin, dim, i, j, k, p, iteration, steps, n_iterations
  integer (c_int)     :: unit_mag, unit_ph, unit_w
  integer (c_int), dimension(:), allocatable  :: base_state, config

  real(c_double), dimension(:), allocatable :: Jint, h_x, h_z
  real(c_double) :: T0, T1, h_coupling, kick 
  
  real (c_double) :: mag, mag_avg, norm
  complex (c_double_complex) :: alpha, beta

  real (c_double), dimension(:), allocatable :: ESwap
  real (c_double), dimension(:,:), allocatable :: HSwap, WSwap

  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:), allocatable :: state, init_state
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap

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
  p = nspin/2

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

  write (*,*) "Number of Iterations"
  read (*,*) n_iterations
  print*,""

  !Standard Values (w/ J = h_z)
  T0 = 1
  h_coupling = 0.3
  kick = 0.1

  write (*,*) "Coupling Constant, h = h_x / J"
  read (*,*) h_coupling
  print*,""
  
  write (*,*) "Perturbation on Kick, epsilon = T1 - pi/2"
  read (*,*) kick
  print*,""

  T1 = pi/4+kick

    !Coefficienti dello stato iniziale |psi> = (alpha|up>+beta|down>)^L
  alpha = cos(pi/8)
  beta = sin(pi/8)

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !Data Files

  !filestring = 'data.txt'

  !Standard Parameters
!  write(filestring,91) "data/magnetizations/Swap_Sz_nspin", nspin, "_steps", steps, &
!  &  "_iterations", n_iterations, ".txt"
!  open(newunit=unit_mag,file=filestring)
!
!  write(filestring,91) "data/eigenvalues/Swap_PH_nspin", nspin, "_steps", steps, &
!  &  "_iterations", n_iterations, ".txt"
!  open(newunit=unit_ph, file=filestring)
!  
!  write(filestring,91) "data/eigenvalues/Swap_W_nspin", nspin, "_steps", steps, &
!  &  "_iterations", n_iterations, ".txt"
!  91  format(A,I0,A,I0,A,I0,A)

  !Input Parameters
  write(filestring,92) "data/magnetizations/Swap_Sz_nspin", nspin, "_steps", steps, &
  &  "_iterations", n_iterations, "_h", h_coupling, "_kick", kick, ".txt"
  open(newunit=unit_mag,file=filestring)

  write(filestring,92) "data/eigenvalues/Swap_PH_nspin", nspin, "_steps", steps, &
  &  "_iterations", n_iterations, "_h", h_coupling, "_kick", kick, ".txt"
  !open(newunit=unit_ph, file=filestring)
  
  write(filestring,92) "data/eigenvalues/Swap_W_nspin", nspin, "_steps", steps, &
  &  "_iterations", n_iterations, "_h", h_coupling, "_kick", kick, ".txt"
  92  format(A,I0, A,I0, A,I0, A,F4.2, A,F4.2, A)

  !open(newunit=unit_w, file=filestring)
 
  !------------------------------------------------

  allocate( Jint(nspin-1), h_x(nspin), h_z(nspin))
  allocate(U(dim,dim), PH(dim), W(dim,dim))
  allocate(state(dim), init_state(dim))
  allocate(HSwap(dim,dim), ESwap(dim), WSwap(dim,dim), USwap(dim,dim))
  !allocate(PH(dim), W(dim,dim))
  !allocate(mag_avg(steps))


  call buildStaggState(nspin, dim, alpha, beta, init_state)

  do iteration = 1, n_iterations
    
    if (mod(iteration,10)==0) then 
      write (*,*) "iteration = ", iteration
    endif

    !Parametri
  
    call random_number(Jint)
    Jint = (Jint + 0.5) !Jint in [1/2,3/2], J == 1
    !Jint = 2*pi
  
    call random_number(h_x)
    h_x = h_coupling*h_x !h_x in [0, 0.3*J]
  
    call random_number(h_z)
    !h_z = 0
  
    !write (*,*) "Jint = ", Jint(:)
    !write (*,*) "h_x = ", h_x(:)
    !write (*,*) "h_z = ", h_z(:)
    !print *, ""
  
    !----------------------------
  
    !Costruzione Hamiltoniane

!    call buildHSwap(nspin, dim, HSwap)
!
!    print *, "HSwap = "
!    call printmat(dim, HSwap, 'R')
!
!    call diagSYM( 'V', dim, HSwap, ESwap, WSwap )
!    call expSYM( dim, -C_UNIT*T1, ESwap, WSwap, USwap )
!
!    print *, "USwap = "
!    call printmat(dim, USwap, 'R')

  
    call buildUFSwap(nspin, dim, Jint, h_x, h_z, T0, T1, U )
!    print *, "UF = "
!    call printmat(dim, U,'C')

!    call diagUN( SELECT, dim, U, PH, W)
  
    !Verify WW^\dagger = 1
!    call printmat(dim, matmul(W,transpose(conjg(W))),'R')
!  
!    print *, "Eigenvalues of U_F "
!    print "(*(/f15.10spf15.10' i'))", PH(:)
!    print *,""
!  
!    print *, "Eigenvectors of U_F "
!    do i = 1,dim
!      write (*,"(*('|',f5.2spf5.2' i'))") W(i,:)
!    enddo
!    print *,""
    !--------------------------------
  
    !Print data to file
    !call writevec(unit_ph,dim,PH,'C')
    !call writemat(unit_w,dim,W,'C')
  
    !-------------------------------------
  
    !Evolution, Magnetization
   
 
    state = init_state
    norm = dot_product(state,state)
    !call printvec(dim, state, 'R')
    j = 1
    
    write(unit_mag,*) "iteration = ", iteration
    write(unit_mag,*) mag_stag_z(nspin, dim, state), j
    
    !print *, mag_stag_z(nspin, dim, state), j, norm
  
    do j = 2, steps
      state = matmul(U,state)
      norm = dot_product(state,state)
      state = state / sqrt(norm)
      !call print_mag(nspin,dim,state,mag)
      write(unit_mag,*) mag_stag_z(nspin, dim, state), j
      !print *, mag_stag_z(nspin, dim, state), j, norm
    enddo
    !print *, ""
 

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
