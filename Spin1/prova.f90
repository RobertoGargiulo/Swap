program prova

  use functions, only : dimSpin1_Sz0
  use matrices, only: buildSz0_HSwap, buildSz0_HMBL, &
    & print_hamiltonian_Sz0, print_unitary_Sz0
  use states, only: buildNeelState_Sz0, printstate_Sz0
  use printing
  use exponentiate, only: diagSYM, expSYM, diagUN
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)
  integer (c_int), parameter :: dimSpin1 = 3

  integer (c_int) :: nspin, dim, dim_Sz0, Sz
  integer (c_int), allocatable :: config(:), idxSz0(:)
  integer (c_int) :: i, j, k, n 

  complex (c_double_complex) :: alpha(dimSpin1)

  complex (c_double_complex), allocatable :: init_state(:), state(:)
  integer (c_int) :: n_disorder, steps
  real (c_double) :: hz_coupling, J_coupling, V_coupling, T0, T1, kick
  real (c_double), allocatable :: H(:,:), hz(:), Jxy(:), Vz(:)
  real (c_double), allocatable :: E(:), W_r(:,:)
  complex (c_double_complex), allocatable :: USwap(:,:), PH(:), W(:,:), U(:,:)

  real (c_double), allocatable :: imb_avg(:), imb_sq(:)
  real (c_double), allocatable :: sigmaz_avg(:,:), sigmaz_sq(:,:)

  integer(c_int) :: count_beginning, count_end, count_rate

  character(len=200) :: filestring_Neel, filestring_imb, filestring_sigmaz

  logical :: SELECT
  EXTERNAL SELECT


  !------------ Input ----------------!

  write (*,*) "Number of Spins"
  read (*,*) nspin
  print *, ""

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

  call system_clock(count_beginning, count_rate)

  dim = dimSpin1**nspin
  dim_Sz0 = dimSpin1_Sz0(nspin)
  T1 = pi/2 + kick

  !-------- Swap Operator -------------!
  allocate(H(dim_Sz0,dim_Sz0))

  call buildSz0_HSwap(nspin, dim_Sz0, H)
  call print_hamiltonian_Sz0(nspin, dim_Sz0, H)

  allocate(E(dim_Sz0), W_r(dim_Sz0,dim_Sz0), USwap(dim_Sz0,dim_Sz0))
  call diagSYM( 'V', dim_Sz0, H, E, W_r)
  deallocate(H)

  call expSYM( dim_Sz0, -C_UNIT*pi/2.d0, E, W_r, USwap )
  deallocate(E, W_r)
  call print_unitary_Sz0(nspin, dim_Sz0, USwap)

  !------------- Initial State --------------!
  allocate(init_state(dim_Sz0))
  call buildNeelState_Sz0(nspin, dim_Sz0, init_state)
  call printstate_Sz0(nspin, dim_Sz0, init_state, "Neel state:")


!  !Allocate state, evolution operators, observables
!  allocate(state(dim_Sz0))
!
!  allocate(hz(nspin), Jxy(nspin-1), Vz(nspin-1))
!  allocate(U(dim_Sz0,dim_Sz0), H(dim_Sz0,dim_Sz0), E(dim_Sz0), W_r(dim_Sz0,dim_Sz0))
!
!  allocate(imb_avg(steps), imb_sq(steps))
!  allocate(sigmaz(nspin), sigmaz_avg(steps,nspin), sigmaz_sq(steps,nspin))
!  
!  do i = 1, n_disorder
!
!    !call random_number(Jxy)
!    call random_number(Vz)
!    call random_number(hz)
!
!    hz = h * 2*(hz-0.5) !hz in [-h, h]
!    Jxy = -J
!    Vz = -V + V * (Vz - 0.5) !Vz in [-3/2V, -V/2] <=> -V \pm V/2
!
!    call buildSz0_HMBL(nspin, dim_Sz0, Jxy, Vz, hz, H)
!
!    !allocate(PH(dim_Sz0),W(dim_Sz0,dim_Sz0))
!    !call diagUN( SELECT, dim_Sz0, U, PH, W)
!
!    state = init_state
!    j = 1
!    sigmaz_avg(j) = sigmaz_avg(j) + sigmaz_Sz0(nspin, dim_Sz0, state)
!    sigmaz_sq(j) = sigmaz_sq(j) + sigmaz_Sz0(nspin, dim_Sz0, state) ** 2
!    do j = 2, steps
!
!      state = matmul(U, state)
!      sigmaz_avg(j) = sigmaz_avg(j) + sigmaz_Sz0(nspin, dim_Sz0, state)
!      sigmaz_sq(j) = sigmaz_sq(j) + sigmaz_Sz0(nspin, dim_Sz0, state) ** 2
!
!
!
!    enddo
!
!  enddo
!  sigmaz_avg = sigmaz_avg / n_disorder



  call take_time(count_rate, count_beginning, count_end, 'T', "Program")

  print *, "End Program"



end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
