program test_LR

  use functions, only: init_random_seed, dimSpin1_Sz, &
    & norm_V => normalization_power_law, &
    & project => projectState_FullHS_to_Sz
  use exponentiate, only: diagSYM, expSYM
  use observables, only: sigmaz_Sz
  use matrices, only: buildHSwap => buildSz_HSwap, buildHMBL => buildSz_HMBL_LR
  use states, only: buildstate => buildUpZeroState, &
    & printstate_sz
  use printing, only: take_time, printmat
  use sorts, only: sort => dpquicksort
  use omp_lib
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  implicit none

  complex (dcp), parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), parameter :: C_ONE = dcmplx(1._dp, 0._dp)
  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  integer (ip), parameter :: dimSpin1 = 3
  character(len=*), parameter :: name_initial_state = "UpZero"
  real (dp), parameter :: pi = 4._dp * atan(1._dp)

  integer (ip)     ::  nspin, dim, n_disorder, steps
  integer (ip)     ::  i, j, l, k, p, q, m, dim_Sz, Sz
  integer (ip)     ::  unit_sigmaz

  real(dp), dimension(:), allocatable :: Jxy, hz
  real(dp), dimension(:,:), allocatable :: Vzz
  real(dp) :: T0, T1, Jxy_coupling, Vzz_coupling, hz_coupling, kick, alpha

  real (dp), dimension(:), allocatable :: E
  real (dp), dimension(:,:), allocatable :: H, W_r
  complex(dcp), dimension(:,:), allocatable :: U, W, USwap

  real (dp) :: norm
  complex (dcp), allocatable :: psi(:), state(:), psi_swap(:), psi_Sz(:)
  real (dp), allocatable :: sigmaz_avg(:,:), sigmaz_sq(:,:)

  logical :: SELECT
  EXTERNAL SELECT

  integer(ip) :: count_beginning, count_end, count_rate
  character(len=200) :: filestring
  character(len=:), allocatable :: columns

  !--------------- (Explicit) Interfaces ------------!
  interface
    function column_titles(nspin) result(columns)
      integer, intent(in) :: nspin
      character(26*(4*nspin+3)) :: columns
    end function
  end interface

  !Parametri Modello: J, V, hz, T0, T1/epsilon, nspin/L
  !Parametri Simulazione: Iterazioni di Disordine

  !Parametri Iniziali

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print*,""

  write (*,*) "Number of Disorder Realizations"
  read (*,*) n_disorder
  print*,""

  write (*,*) "Number of Steps"
  read (*,*) steps
  print*,""

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

  write (*,*) "Perturbation on Kick T1 = pi/2 + kick"
  read (*,*) kick
  print*,""

  write (*,*) "Power-Law Coefficient Vzz = V_{ij} / |i-j|^alpha"
  read (*,*) alpha
  print*,""
 
  dim = dimSpin1**nspin
  T1 = pi/2 + kick
  !Sz = 0 !!!!!!! <------- Magnetization of Subspace
  !dim_Sz = dimSpin1_Sz(nspin, Sz)

  call system_clock(count_beginning, count_rate)
  !---------------------------------------------

  !------------------- File Manager ------------------------!

  write(filestring,93) "data/dynamics/sigmaz_Swap_LR_" // trim(name_initial_state) // "_nspin", &
    & nspin, "_steps", steps, "_period", T0, "_n_disorder", n_disorder, &
    & "_J", Jxy_coupling, "_V", int(Vzz_coupling), Vzz_coupling-int(Vzz_coupling), &
    & "_hz", int(hz_coupling), hz_coupling-int(hz_coupling), "_kick", kick, &
    & "_alpha", int(alpha), alpha-int(alpha), ".txt"
  print *, "Output: ", trim(filestring)

  93  format(A,I0, AI0, A,F4.2, A,I0, A,F7.5, A,I0,F0.2, A,I0,F0.2, A,F5.3, A,I0,F0.2, A)

  open(newunit=unit_sigmaz,file=filestring)
  call write_info(unit_sigmaz, trim(name_initial_state))


  !------------- Initial State ------------------!
  allocate(psi(dim))
  call buildstate(nspin, dim, psi)
  call project(nspin, dim, psi, dim_Sz, Sz, psi_Sz)
  call printstate_Sz(nspin, dim_Sz, Sz, psi_Sz, trim(name_initial_state))
 
  !------------------------------------------------

  !BUILD DRIVING PROTOCOL (NO DISORDER) USwap = exp(-i*(pi/4 + eps)*HSwap)
  allocate(H(dim_Sz,dim_Sz), E(dim_Sz), W_r(dim_Sz,dim_Sz), USwap(dim_Sz,dim_Sz))
  call buildHSwap(nspin, dim_Sz, Sz, H)
  call diagSYM( 'V', dim_Sz, H, E, W_r)
!  print *, "HSwap = "
!  call printmat(dim, H, 'R')
  deallocate(H)
  call expSYM( dim_Sz, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)
  !print *, "USwap = "
  !call printmat(dim_Sz, USwap, 'C')

  !------ Allocate state, evolution operators, observables -----!
  allocate(psi_swap(dim_Sz))

  allocate( Jxy(nspin-1), Vzz(nspin-1,nspin), hz(nspin))
  allocate(U(dim_Sz,dim_Sz), H(dim_Sz,dim_Sz), E(dim_Sz), W_r(dim_Sz,dim_Sz))

  allocate(sigmaz_avg(steps,nspin), sigmaz_sq(steps,nspin))
  sigmaz_avg = 0
  sigmaz_sq = 0

  !$OMP PARALLEL
  call init_random_seed() 
  print *, "Size of Thread team: ", omp_get_num_threads()
  !print *, "Verify if current code segment is in parallel: ", omp_in_parallel()
  !$OMP do reduction(+: sigmaz_avg, sigmaz_sq) private(i, hz, Vzz, norm, j, &
  !$OMP & psi_swap, H, E, W_r, U )
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
 
    !call random_number(Jxy)
    !Jxy = 2*Jxy_coupling*(Jxy - 0.5) !Jxy in [-J,J]
    Jxy = -Jxy_coupling

    !Vzz = -Vzz_coupling
    call random_number(Vzz)
    Vzz = -Vzz_coupling + Vzz_coupling*(Vzz - 0.5) !Vzz in [-V-V/2,-V+V/2]
    !Vzz = 1
    do k = 1, nspin-1
      Vzz(k,1:k) = 0
      do q = k+1, nspin
        Vzz(k,q) = Vzz(k,q) / ( abs(k-q)**alpha )
      enddo
      !write (*,*) Vzz(k,:)
    enddo
    Vzz = Vzz / norm_V(alpha, nspin)
 
    call random_number(hz)
    hz = 2*hz_coupling*(hz-0.5) !hz in [-hz_coupling, hz_coupling]
 
    !write (*,*) "Jxy = ", Jxy(:)
    !write (*,*) "Vzz = ", Vzz(:)
    !write (*,*) "hz = ", hz(:)
    !print *, ""
 
    !------------ Floquet Operator(s) --------------!
    call buildHMBL(nspin, dim_Sz, Sz, Jxy, Vzz, hz, H)
    call diagSYM( 'V', dim_Sz, H, E, W_r)
    call expSYM( dim_Sz, -C_UNIT*T0, E, W_r, U )
    U = matmul(USwap,U)

    psi_swap = psi_Sz
    j = 1
    sigmaz_avg(j,:) = sigmaz_avg(j,:) + sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap)
    sigmaz_sq(j,:) = sigmaz_sq(j,:) + sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap) ** 2
    do j = 2, steps

      if(mod(j,10)==0) then
        norm = real(dot_product(psi_swap,psi_swap))
        psi_swap = psi_swap/sqrt(norm)
      endif
      psi_swap = matmul(U, psi_swap)
      sigmaz_avg(j,:) = sigmaz_avg(j,:) + sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap)
      sigmaz_sq(j,:) = sigmaz_sq(j,:) + sigmaz_Sz(nspin, dim_Sz, Sz, psi_swap) ** 2

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL 

  sigmaz_avg = sigmaz_avg / n_disorder
  sigmaz_sq = sqrt( (sigmaz_sq/n_disorder - sigmaz_avg**2) / n_disorder )

  columns = column_titles(nspin)
  write (*,*) columns
  do j = 1, steps, max(1,steps/100)
    write (*,*) j*T0, sigmaz_avg(j,:), sigmaz_sq(j,:)
  enddo

  write (unit_sigmaz,*) columns
  do j = 1, steps
    write (unit_sigmaz,*) j*T0, sigmaz_avg(j,:), sigmaz_sq(j,:)
  enddo

  deallocate(sigmaz_avg, sigmaz_sq)
  deallocate(psi_swap, psi_Sz)
  deallocate(USwap, H, E, W_r, U)
  deallocate(Vzz, hz)

  close(unit_sigmaz)

  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end

subroutine write_info(unit_file, state_name)

  integer, intent(in) :: unit_file
  character(len=*) :: state_name

  write (unit_file,*) "Some info: "
  write (unit_file,*) "Dynamics of disorder average of magnetization at integer multiples of the period."
  write (unit_file,*) "Spin-1 chain with hamiltonian H = sum hz * Z + V * ZZ + J * (XX + YY)."
  write (unit_file,*) "Periodic perturbed swap, U_swap = exp(-i(pi/2 + kick) * sum (sigma*sigma)^2 + (sigma*sigma) - 2 )."
  write (unit_file,*) "V = V_{ij}/|i-j|^alpha, with V_{ij} is taken in [-3V/2, -V/2] (OBC); &
  &  hz is taken in [-hz, hz];  J is uniform."
  write (unit_file,*) "Exact diagonalization of the dense matrix has been used to compute U_F and diagonalize it."
  write (unit_file,*) "Initial state is "//trim(state_name)//trim(".")

end subroutine

function column_titles(nspin) result(columns)

  use iso_c_binding
  implicit none

  integer (c_int), intent(in) :: nspin
  character (26*(4*nspin+3)) :: columns
  character (nspin) :: s_index
  integer (c_int) :: i, j1, j2, j3, j4


  columns = ""
  columns(2:26+2) = "j*T0"
  do i = 1, nspin
    write(s_index, "(I0)") i
    j1 = 26*i+2
    j2 = 26*(i+1)+2
    j3 = 26*(i+nspin)+3
    j4 = 26*(i+nspin+1)+3
    !print *, j1, j2, j3, j4
    columns(j1:j2) = trim("sigma_z^")//trim(s_index)
    columns(j3:j4) = trim("err(sigma_z^")//trim(s_index)//trim(")")
    !print *, columns
  enddo
  !print *, columns
  !write(columns,"(3X,A4,20X,A)") "j*T0", trim(columns)
  !print *, columns

end function
