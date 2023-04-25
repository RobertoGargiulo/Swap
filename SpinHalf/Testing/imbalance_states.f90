program prova

  use iso_c_binding
  use printing
  use sorts
  use genmat
  use MBL
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)
  real (c_double), parameter :: pi = 4.d0 * datan(1.d0)

  integer (c_int) :: nspin, dim, dim_Sz0, i, j, k, l
  integer (c_int) :: count_rate, count_start, count_end
  real (c_double), dimension(:), allocatable :: E, IMB_arr
  complex (c_double_complex), dimension(:), allocatable :: PH, psi
  real (c_double) :: IMB, LI, tol
  integer (c_int), allocatable :: states(:), idxSz0(:), config(:)
  integer (c_int) :: int_1d, int_2d(2), idx

do nspin = 8, 12, 2
  print *, "nspin = ", nspin
  dim = 2**nspin
  dim_Sz0 = binom(nspin,nspin/2)

  call system_clock(count_start, count_rate)

  allocate(states(dim_Sz0))

  tol = 1.0e-10

  allocate(IMB_arr(nspin/2 + 1))
  allocate( idxSz0(dim_Sz0), config(nspin) )
  call zero_mag_states(nspin, dim_Sz0, idxSz0)
  IMB_arr = 0
  k = 0
  idx = 0
  print "(4X,A4, 4X,A6, 4X,A6, 2X,A)", "l", "I", "LI", "config"
  do i = 1, dim_Sz0

    
    l = idxSz0(i)
    call decode(l, nspin, config)
    IMB = imbalance_basis(nspin, l)
    LI = local_imbalance_basis(nspin, l)
    !print "(4X,I4, 4X,F6.3, 4X,F6.3, 2X,*(I0))", l, IMB, LI, config(:)

    if (idx == 0 .AND. abs(IMB) < tol) then
      idx = 1
      k = k+1
      IMB_arr(k) = IMB
    else if ( ANY( abs(IMB_arr(1:k)-IMB) < tol ) ) then
    else
      k = k+1
      IMB_arr(k) = IMB
    endif

  enddo
  
  call dpquicksort(IMB_arr(1:k))
  do i = 1, k
    !print *, "Imbalance = ", IMB_arr(i)
    !call finite_imbalance_states_Sz0(nspin, dim_Sz0, IMB_arr(i), states)
  enddo
  print *, ""



  IMB = 0.5
  LI = 0.5
  call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, states)

  allocate(psi(dim_Sz0))
  call buildRndLarge_IMB_LI_State_basis_Sz0(nspin, dim_Sz0, IMB, LI, psi)
  call buildRndLarge_IMB_LI_State_Sz0(nspin, dim_Sz0, IMB, LI, psi)
  call buildLarge_IMB_LI_State_basis_Sz0(nspin, dim_Sz0, 3, IMB, LI, psi)
  deallocate(psi)
  !deallocate(states)
  deallocate(IMB_arr, idxSz0, config, states )

  call take_time(count_rate, count_start, count_end, 'T', "Program")
  print *, ""

enddo


end program prova



