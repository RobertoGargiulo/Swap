program test_Sz

  use functions, only: basis_Sz, dimSpin1_Sz, project => projectState_FullHS_to_Sz
  use states, only: buildNeelState, printstate_Sz, printstate, buildUpZeroState
  use matrices, only: buildSz_HSwap, buildSz_HMBL, &
    & print_hamiltonian_Sz, print_unitary_Sz
  use exponentiate, only: diagSYM, expSYM, diagUN
  use observables, only: exact_energies_Sz, exact_quasi_energies_Sz
  use sorts, only: sort => dpquicksort
  use iso_c_binding, only: ip => c_int, dp => c_double, cp => c_double_complex
  implicit none

  integer (ip), parameter :: dimSpin1 = 3
  complex (cp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)
  real (dp), parameter :: pi = 4.d0 * datan(1.d0)
  integer (ip) :: nspin, dim_Sz, Sz, dim
  integer (ip), allocatable :: idxSz(:)

  complex (cp), allocatable :: psi(:), psi_Sz(:)
  complex (cp), allocatable :: USwap(:,:), U(:,:), PH(:), W(:,:)
  real (dp), allocatable :: H(:,:), E(:), W_r(:,:)
  real (dp), allocatable :: E_MBL(:), QE(:), QE_ex(:)
  real (dp), allocatable :: Jxy(:), Vz(:), hz(:)
  real (dp) :: T1, T0
  integer :: i

  logical :: SELECT
  EXTERNAL SELECT

  write (*,*) "nspin = "
  read (*,*) nspin
  print *, ""
  dim = dimSpin1 ** nspin

  !dim = 0
  !do Sz = -nspin, nspin
  !  dim_Sz = dimSpin1_Sz(nspin, Sz)
  !  dim = dim + dim_Sz
  !  print *, "Sz = ", Sz, ", dim_Sz = ", dim_Sz
  !  allocate(idxSz(dim_Sz))
  !  call basis_Sz(nspin, dim_Sz, Sz, idxSz)
  !  deallocate(idxSz)
  !  print *, ""
  !  print *, "(16 * dim_Sz)^2 / 10^9 = ", (16_dp*dim_Sz)**2 / 10_dp**9
  !  print *, ""
  !enddo
  !print *, "dim = ", dim, ", dimSpin1**nspin = ", dimSpin1**nspin
  !print *, ""

  !------------ Initial State --------------!

  allocate(psi(dim))
!  call buildNeelState(nspin, dim, psi)
!  call printstate(nspin, dim, psi, "Neel state:")
!
!  call project(nspin, dim, psi, dim_Sz, Sz, psi_Sz)
!
!  call printstate_Sz(nspin, dim_Sz, Sz, psi_Sz, "Projected Neel state:")

  call buildUpZeroState(nspin, dim, psi)
  call printstate(nspin, dim, psi, "UpZero state:")

  call project(nspin, dim, psi, dim_Sz, Sz, psi_Sz)
  call printstate_Sz(nspin, dim_Sz, Sz, psi_Sz, "Projected UpZero state:")


  !-------- Swap --------!

  allocate(H(dim_Sz, dim_Sz))
  call buildSz_HSwap(nspin, dim_Sz, Sz, H)
  print *, "H_Swap = "
  call print_hamiltonian_Sz(nspin, dim_Sz, Sz, H)

  allocate(E(dim_Sz), W_r(dim_Sz,dim_Sz), USwap(dim_Sz,dim_Sz))
  call diagSYM( 'V', dim_Sz, H, E, W_r)
  deallocate(H)

  T1 = pi/2
  call expSYM( dim_Sz, -C_UNIT*T1, E, W_r, USwap )
  deallocate(E, W_r)
  print *, "U_Swap = "
  call print_unitary_Sz(nspin, dim_Sz, Sz, USwap)
  
  !-------- Hamiltonian ----------!
  allocate(Jxy(nspin-1), Vz(nspin-1), hz(nspin))
  Jxy = 0
  call random_number(Vz)
  call random_number(hz)
  Vz = 2 * Vz
  hz = 5 * 2 * (hz - 0.5)
  T0 = 1

  allocate(H(dim_Sz, dim_Sz))
  call buildSz_HMBL(nspin, dim_Sz, Sz, Jxy, Vz, hz, H)
  print *, "H_MBL = "
  call print_hamiltonian_Sz(nspin, dim_Sz, Sz, H)

  allocate(E(dim_Sz), W_r(dim_Sz,dim_Sz), E_MBL(dim_Sz))
  call diagSYM( 'V', dim_Sz, H, E, W_r)
  E_MBL = exact_energies_Sz(nspin, dim_Sz, Sz, Vz, hz)
  call sort(E_MBL)
  print *, "E = "
  do i = 1, dim_Sz
    print *, E_MBL(i), E(i)
  enddo
  print *, ""

  allocate(U(dim_Sz,dim_Sz),PH(dim_Sz),W(dim_Sz,dim_Sz))
  call expSYM( dim_Sz, -C_UNIT*T0, E, W_r, U )
  print *, "U_MBL = "
  call print_unitary_Sz(nspin, dim_Sz, Sz, U)
  U = matmul(USwap,U)
  print *, "U_F = U_Swap * U_MBL = "
  call print_unitary_Sz(nspin, dim_Sz, Sz, U)
  deallocate(E, W_r)
  call diagUN( SELECT, dim_Sz, U, PH, W)

  allocate(QE(dim_Sz), QE_ex(dim_Sz))
  QE = real(C_UNIT*log(PH), dp)
  QE_ex = exact_quasi_energies_Sz(nspin, dim_Sz, Sz, Vz, hz)
  call sort(QE)
  call sort(QE_ex)
  print*, "QE = "
  do i = 1, dim_Sz
    print *, QE_ex(i), QE(i), QE(i) + pi, QE(i) - pi
  enddo
  print *, ""

  U = matmul(U,U)
  print *, "U_F^2 = "
  call print_unitary_Sz(nspin, dim_Sz, Sz, U)


end program

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end
