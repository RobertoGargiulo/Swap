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
  complex (c_double_complex), dimension(:), allocatable :: PH
  real (c_double) :: IMB, LI, tol
  integer (c_int), allocatable :: states(:), idxSz0(:), config(:)
  integer (c_int), allocatable :: n_imb(:), IMB_tot(:)
  integer (c_int) :: int_1d, int_2d(2), idx

do nspin = 2, 4, 2
  print *, "nspin = ", nspin
  !read (*,*) nspin
  dim = 2**nspin
  dim_Sz0 = binom(nspin,nspin/2)

  !print *, "Imbalance = "
  !read (*,*) IMB

  call system_clock(count_start, count_rate)


  !allocate(states(dim))
  !print *, "states allocated"

  tol = 1.0e-10

  !IMB_arr = 0
  !k = 0
  !do i = 1, (nspin/2 + 1)**2

  !  int_2d = int_1dto2d(i-1)
  !  !print *, int_2d
  !  IMB = real(int_2d(2) - int_2d(1),c_double) / real(nspin - int_2d(2) - int_2d(1), c_double)
  !  if(int_2d(2) == int_2d(1)) IMB = 0
  !  !print *, "Imbalance = ", IMB

  !  idx = 0
  !  if ( ANY( abs(IMB_arr-IMB) < tol ) ) then
  !    idx = 0
  !  else
  !    idx = 1
  !  endif

  !  if (idx == 1) then
  !    k = k+1
  !    !print *, k, IMB
  !    IMB_arr(k) = IMB
  !  endif

  !  !call finite_imbalance_states(nspin, dim, IMB, states)
  !enddo
  !print *, "# of Imbalances = ", k+1
  !print *, ""

  allocate(IMB_arr(nspin/2 + 1))
  allocate( idxSz0(dim_Sz0), config(nspin) )
  allocate(n_imb(dim_Sz0), IMB_tot(dim_Sz0))
  call zero_mag_states(nspin, dim_Sz0, idxSz0)
  IMB_arr = 0
  n_imb = 0
  k = 0
  j = 0
  idx = 0
  do i = 1, dim_Sz0

    !print *, int_2d
    
    l = idxSz0(i)
    call decode(l, nspin, config)
    IMB = imbalance_basis(nspin, l)
    LI = local_imbalance_basis(nspin, l)
    print "(4X,I0, 4X,F6.3, 4X,F6.3, 2X,*(I0))", l, IMB, LI, config(:)

    !IMB_tot(i) = 
    !print *, l, "Imbalance = ", IMB

    if (idx == 0 .AND. abs(IMB) < tol) then
      idx = 1
      k = k+1
      IMB_arr(k) = IMB
      !print *, "Zero IMB ", i, k, IMB
    else if ( ANY( abs(IMB_arr(1:k)-IMB) < tol ) ) then
      !idx = 0
      n_imb(k+1) = n_imb(k+1) + 1
      !print *, "Same IMB ", i, n_imb(k+1), IMB
    else
      !idx = 1
      k = k+1
      !print *, "Different IMB ", i, k, IMB
      IMB_arr(k) = IMB
      !print "(4X,I0, 4X,F6.3, 2X,*(I0))", l, IMB, config(:)
    endif

    !if (idx == 1) then
    !endif

    !call finite_imbalance_states(nspin, dim, IMB, states)
  enddo
  
  call dpquicksort(IMB_arr(1:k))
  do i = 1, k
    print *, "Imbalance = ", IMB_arr(i)!, k==nspin/2+1
  enddo

  !do i = 1, dim_Sz0


  call take_time(count_rate, count_start, count_end, 'T', "Program")

  !deallocate(states)
  deallocate(IMB_arr, idxSz0, config )
  deallocate(n_imb, IMB_tot)

enddo

end program prova
