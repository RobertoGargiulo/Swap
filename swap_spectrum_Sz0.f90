program swap

  use exp_sparse
  use exponentiate
  use genmat
  use printing
  use MBL
  !use omp_lib
  use sorts
  use iso_c_binding
  !use general
  implicit none

  complex (c_double_complex), parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int)     ::  nspin, dim, iteration, steps, n_iterations, dim_eff
  integer (c_int)     ::  i, j, k, p, nz_dim, krylov_dim, dim_Sz0, nz_Sz0_dim
  integer (c_int)     ::  unit_mag, unit_ph, unit_w, unit_avg

  real(c_double), dimension(:), allocatable :: Jint, Vint, h_z
  real(c_double) :: T0, T1, J_coupling, V_coupling, hz_coupling, kick 
  real(c_double) :: t_avg, t_sigma, r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2
  
  real (c_double) :: norm
  real (c_double), dimension(:), allocatable :: avg, sigma, avg2, sigma2, r_avg, r_sigma, r_avg2, r_sigma2
  complex (c_double_complex) :: alpha, beta

  real (c_double), dimension(:), allocatable :: E
  real (c_double), dimension(:,:), allocatable :: H, W_r
  complex(c_double_complex), dimension(:), allocatable :: PH
  complex(c_double_complex), dimension(:,:), allocatable :: U, W, USwap

  real (c_double), dimension(:,:,:), allocatable :: MI
  real (c_double), dimension(:,:), allocatable :: MI_avg
  integer (c_int) :: nspin_A, dim_A

  logical :: SELECT
  EXTERNAL SELECT

  integer(c_int) :: count_beginning, count_end, count_rate, time_min, count1, count2
  real (c_double) :: time_s
  character(len=200) :: filestring
  character(len=8) :: time_string


  !Parametri Modello: J, V, h_z, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine, Steps di Evoluzione, Stato Iniziale

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""
  dim = 2**nspin
  dim_Sz0 = binom(nspin,nspin/2)
  nz_Sz0_dim = non_zero_HMBL_Sz0(nspin)

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
  alpha = 1
  beta = 0

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !DATA FILES
 
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
  !allocate(H_sparse(nz_Sz0_dim), ROWS(nz_Sz0_dim), COLS(nz_Sz0_dim))
  allocate(H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))
  allocate(U(dim_Sz0,dim_Sz0))

  !Allocate observables and averages
  allocate( r_avg(n_iterations), r_sigma(n_iterations))
  allocate( r_avg2(n_iterations), r_sigma2(n_iterations))

  !Allocate for Eigenvalues/Eigenvectors
  allocate(PH(dim_Sz0))
  allocate(W(dim_Sz0,dim_Sz0))
  allocate(MI(nspin/2,n_iterations,dim_Sz0), MI_avg(nspin/2,n_iterations))

  r_avg = 0
  r_sigma = 0
  r_avg2 = 0
  r_sigma2 = 0

  !$OMP PARALLEL
  call init_random_seed() 
  !print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP DO private(iteration, h_z, H, E, W_r, &
  !$OMP & U, W, PH)
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
  
!    write (*,*) "Jint = ", Jint(:)
!    write (*,*) "Vint = ", Vint(:)
!    write (*,*) "h_z = ", h_z(:)
!    print *, ""
  
    !---------------------------------------------------
    !call take_time(count_rate, count_beginning, count1, 'F', filestring)
    !BUILD FLOQUET (EVOLUTION) OPERATOR
    call buildSz0_HMBL( nspin, dim_Sz0, Jint, Vint, h_z, H )
    call diagSYM( 'V', dim_Sz0, H, E, W_r )
    call gap_ratio(dim_Sz0, E, r_avg2(iteration), r_sigma2(iteration))
    call expSYM( dim_Sz0, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)
    call diagUN( SELECT, dim_Sz0, U, PH, W)
    E = real(C_UNIT*log(PH))
    call dpquicksort(E)
    call gap_ratio(dim_Sz0, E, r_avg(iteration), r_sigma(iteration))

    do nspin_A = 1, nspin/2
      dim_A = 2**nspin_A
      do i = 1, dim_Sz0
        MI(nspin_A,iteration,i) = mutual_information_Sz0(nspin, nspin_A, dim, dim_Sz0, dim_A, W(1:dim_Sz0,i))
      enddo
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 


  call time_avg('F', n_iterations, 1, r_avg, r_sigma, r_dis_avg, r_dis_sigma)
  call time_avg('F', n_iterations, 1, r_avg2, r_sigma2, r_dis_avg2, r_dis_sigma2)

  print *, "Average and Variance of Gap Ratio (over the spectrum and then disorder)"
  print *, r_dis_avg, r_dis_sigma, r_dis_avg2, r_dis_sigma2

  MI_avg = 0
  print *, "Mutual Information (MI, log 2, nspin_A, iteration)"
  do nspin_A = 1, nspin/2
    do iteration = 1, n_iterations
      do i = 1, dim_Sz0
        MI_avg(nspin_A,iteration) = MI_avg(nspin_A,iteration) + MI(nspin_A,iteration,i)
      enddo
      MI_avg = MI_avg/dim_Sz0
      print *, MI_avg(nspin_A,iteration), log(2.d0), nspin_A, iteration
    enddo
  enddo

  call take_time(count_rate, count_beginning, count1, 'T', "Program")

  deallocate(Jint, Vint, h_z)
  deallocate(H,E,W_r)
  !deallocate(USwap)
  !deallocate(U, PH, W)

end program swap

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end

subroutine take_time(count_rate, count_start, count_end, opt, filestring)
  implicit none
  integer, intent(in) :: count_rate, count_start
  integer, intent(out) :: count_end
  character :: opt*1
  character (*) :: filestring

  real :: time_s
  integer :: time_min

  call system_clock(count_end)

  time_s = real(count_end - count_start) / real(count_rate)
  time_min = int(time_s/60)
  if(opt == 'T') then
    print "(A,A,A,1X,I4,A,2X,F15.10,A)", "Elapsed Time for ", filestring, ": ", time_min, "min", time_s - 60*time_min, "sec"
  else if(opt == 'F') then
    print *, ""
  endif
  !print *, ""
end subroutine
