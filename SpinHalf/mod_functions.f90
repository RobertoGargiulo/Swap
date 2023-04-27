!This module contains various functions and procedures used throughout the project

!Specifically it includes all decoding/encoding procedures to go 
!   from an integer to a basis (full Hilbert space, Sz=0 subspace)

!Binomial function, which is used to compute Hilbert space dimensions

!Some procedures which compute the imbalance, local imbalance of the spin basis

!A procedure which generates the full Hilbert-space state corresponding to
!   the Sz=0 subspace

!Other utilities: binary search algorithm, conversion from 1D grid to 2D grid

module functions
  
  !use ifport
  use iso_c_binding, only: dp => c_double, ip => c_int, cp => c_double_complex
  use printing
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  integer (c_int), private, parameter :: dimSpinHalf = 2
  real (c_double), private, parameter :: tol = 1.0e-6
  !integer (c_int) :: i, j, s
  !real (c_double), parameter, private :: SigmaZ(3,3) = reshape((/ (( (2-i)*merge(1,0,i==j), j=1,3),i=1,3) /), shape(SigmaZ)) 

    !SigmaZ = 0
    !do i = 1, 3
    !  s = 2 - i
    !  SigmaZ(i,i) = s
    !enddo

contains

  subroutine decode(decimal, nbits, bitstring)

    integer (c_int), intent(in) :: decimal, nbits
    integer (c_int), intent(out) :: bitstring(nbits)
    integer (c_int) :: i

    bitstring = 0
    do i = 1, nbits
      if (btest(decimal, i - 1)) bitstring(i) = 1
    enddo
  end subroutine decode

  subroutine encode(bitstring, nbits, intgr)

    integer (c_int), intent(in) :: nbits
    integer (c_int), intent(in) :: bitstring(nbits)
    integer (c_int), intent(out) :: intgr
    integer (c_int) :: i

    intgr = 0
    do i = 1, nbits
      intgr = intgr !+ bistring(i) * 2**(i-1)
    enddo
  end subroutine encode

  function swap_config(nspin, i) result(j)

    integer (c_int), intent(in) :: nspin, i
    integer (c_int) :: config(nspin), k, j

    call decode(i, nspin, config)
    j = 0
    do k = 1, nspin/2
      j = j + config(2*k) * 2**(2*k-2)  + config(2*k-1) * 2**(2*k-1)
    enddo

  end function 

  subroutine init_random_seed()
   use iso_fortran_env, only: int64
   implicit none
   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid
   integer(int64) :: t
 
   call random_seed(size = n)
   allocate(seed(n))
   ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
      read(un) seed
      close(un)
   else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(t)
      if (t == 0) then
         call date_and_time(values=dt)
         t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
              + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
              + dt(3) * 24_int64 * 60 * 60 * 1000 &
              + dt(5) * 60 * 60 * 1000 &
              + dt(6) * 60 * 1000 + dt(7) * 1000 &
              + dt(8)
      end if
      pid = getpid()
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
         seed(i) = lcg(t)
      end do
   end if
   call random_seed(put=seed)
 contains
   ! This simple PRNG might not be good enough for real work, but is
   ! sufficient for seeding a better PRNG.
   function lcg(s)
     integer :: lcg
     integer(int64) :: s
     if (s == 0) then
        s = 104729
     else
        s = mod(s, 4294967296_int64)
     end if
     s = mod(s * 279470273_int64, 4294967291_int64)
     lcg = int(mod(s, int(huge(0), int64)), kind(0))
   end function lcg
 end subroutine init_random_seed


  integer function binom(n,k)
    integer(c_int), intent(in) :: n,k
    integer, parameter :: i8 = selected_int_kind(18)
    integer, parameter :: dp = selected_real_kind(15)

    if (k == n) then
        binom = 1
    else if (k == 1) then
        binom = n
    else if ((k /= 1) .and. (k /= n)) then
      binom = nint(exp(log_gamma(n+1.0_dp)-log_gamma(n-k+1.0_dp)-log_gamma(k+1.0_dp)),kind=i8)
    end if 
  end function




  subroutine buildState_Sz0_to_FullHS(nspin, dim, dim_Sz0, psi_Sz0, psi)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (c_int), intent(in) :: nspin, dim, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    complex (c_double_complex), intent(out) :: psi(dim)

    integer :: i, l, states(dim_Sz0)

    call zero_mag_states(nspin, dim_Sz0, states)

    psi = 0
    !print *, "psi initialized"
    do i = 1, dim_Sz0
      l = states(i) + 1
      psi(l) = psi_Sz0(i)
      !print *,l, psi(l)
      !print *, l, dim, i, dim_Sz0
    enddo

  end subroutine buildState_Sz0_to_FullHS



  integer function non_zero_HMBL_Sz0(L)

    integer, intent(in) :: L
    integer :: nz

    if (L==2) then
      nz = 4
    else if (L==4) then
      nz = 18
    else if (L==6) then
      nz = 80
    else if (L==8) then
      nz = 350
    else if (L==10) then
      nz = 1512
    else if (L==12) then
      nz = 6468
    else if (L==14) then
      nz = 27456
    else if (L==16) then
      nz = 115830
    else if (L==18) then
      nz = 486200
    else if (L==20) then 
      nz = 2032316
    else if (L==22) then
      nz = 8465184
    else
      print *, "Chain Length is odd or too large"
      nz = 0
    endif
    
    non_zero_HMBL_Sz0 = nz

  end function non_zero_HMBL_Sz0



  subroutine zero_mag_states(nspin, dim_Sz0, states)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (c_int), intent(in) :: nspin, dim_Sz0
    integer (c_int), intent(out) :: states(dim_Sz0)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    k = 0
    states = 0
    do i = 0, 2**nspin-1

      call decode(i, nspin, config)

      if (sum(config)==nspin/2) then
        k = k+1
        states(k) = i
        !print *, i, k, states(k)
      endif
    enddo

  end subroutine zero_mag_states

  subroutine basis_Sz0_inv(nspin, dim_Sz0, states, inverse)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (c_int), intent(in) :: nspin, dim_Sz0
    integer (c_int), intent(out) :: states(dim_Sz0), inverse(dimSpinHalf**nspin)

    integer :: config(nspin), spin(nspin)
    integer (c_int) :: i, j, l, dim

    dim = dimSpinHalf**nspin
    l = 0
    states = 0
    inverse = -1
    do i = 0, dim-1

      call decode(i, nspin, config)
      spin = 1 - 2*config

      if (sum(spin)==0) then
        l = l+1
        states(l) = i
        inverse(i+1) = l
        !print "(2(I4,4X),*(I0))", i, l, config(:)
      endif
    enddo
    !print *, "dim_Sz0 = ", dim_Sz0

  end subroutine





  integer (c_int) function binsearch(val, array)
  
  
    implicit none
    integer (c_int), intent(in) :: val, array(:)
    integer (c_int) :: mid, start, finish, range
    
    binsearch = -1
    start = 1
    finish = size(array)
    
    range = finish - start
    mid = (finish + start) / 2
    
    do while (array(mid) /= val .and. range > 0)
    
      if (val > array(mid)) then
        start = mid + 1
      else
        finish = mid - 1
      end if
      
      range = finish - start
      mid = (start + finish) / 2
    
    end do
    
    if (array(mid) == val) then
      binsearch = mid
    end if
  
  end function binsearch

  function binsearch_closest_from_above(val, array) result(indx)

    implicit none
    real (c_double), intent(in) :: val, array(:)
    integer (c_int) :: indx, mid, left, right, range

    indx = -1
    left = 1
    right = size(array) + 1

    range = right - left

    !print "(3(A12),2(A26))", "left", "right", "mid", "array(mid)", "val"

    do while (range > 0)

      mid = (right + left) / 2

      if (val < array(mid)) then
        right = mid
      else
        left = mid + 1
      end if
      range = right - left
      
      !print *, left, right, mid, array(mid), val

    end do

    indx = right! - 1 

    !if (array(mid) == val) then
    !  binsearch = mid
    !end if
 
  end function

  integer function int_2dto1d(int_2d)

    implicit none
    integer, intent(in) :: int_2d(2)
    integer :: int_1d
    

    if (int_2d(1) < int_2d(2)) then
      int_1d = int_2d(2)**2 + int_2d(1)
    else if (int_2d(1) >= int_2d(2)) then
      int_1d = int_2d(1)**2 + int_2d(1) + int_2d(2)
    endif

    int_2dto1d = int_1d

  end function int_2dto1d


  function int_1dto2d(int_1d)

    implicit none
    integer, intent(in) :: int_1d
    integer :: int_2d(2)
    integer :: int_1dto2d(2)

    if (int_1d - floor(sqrt(real(int_1d)))**2 < floor(sqrt(real(int_1d))) ) then
      int_2d(1) = int_1d - floor(sqrt(real(int_1d)))**2
      int_2d(2) = floor(sqrt(real(int_1d)))
    else
      int_2d(1) = floor(sqrt(real(int_1d)))
      int_2d(2) = int_1d - floor(sqrt(real(int_1d)))**2 - floor(sqrt(real(int_1d)))
    endif

    int_1dto2d = int_2d

  end function int_1dto2d

  real (dp) function norm(alpha, nspin)

    real (dp) :: alpha
    integer (ip) :: nspin

    if (alpha > 1) then
      norm = 1
    else if (alpha == 1) then
      norm = log(real(nspin, kind(dp)))
    else
      norm = nspin**(1-alpha)
    endif

  end function



end module functions
