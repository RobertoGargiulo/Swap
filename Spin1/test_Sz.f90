program test_Sz

  use functions, only: basis_Sz, dimSpin1_Sz, project => projectState_FullHS_to_Sz
  use states, only: buildNeelState, printstate_Sz, printstate, buildUpZeroState
  use iso_c_binding, only: ip => c_int, dp => c_double, cp => c_double_complex
  implicit none

  integer (ip), parameter :: dimSpin1 = 3
  integer (ip) :: nspin, dim_Sz, Sz, dim
  integer (ip), allocatable :: idxSz(:)

  complex (cp), allocatable :: psi(:), psi_Sz(:)

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

  allocate(psi(dim))
  call buildNeelState(nspin, dim, psi)
  call printstate(nspin, dim, psi, "Neel state:")

  call project(nspin, dim, psi, Sz, psi_Sz)
  dim_Sz = dimSpin1_Sz(nspin, Sz)

  call printstate_Sz(nspin, dim_Sz, Sz, psi_Sz, "Projected Neel state:")

  call buildUpZeroState(nspin, dim, psi)
  call printstate(nspin, dim, psi, "UpZero state:")

  call project(nspin, dim, psi, Sz, psi_Sz)
  dim_Sz = dimSpin1_Sz(nspin, Sz)

  call printstate_Sz(nspin, dim_Sz, Sz, psi_Sz, "Projected UpZero state:")


end program
