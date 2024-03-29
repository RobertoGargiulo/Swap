!This module contains various functions and procedures used throughout the project

!Specifically it includes all decoding/encoding procedures to go 
!   from an integer to a basis (full Hilbert space, Sz=0 subspace, generic Sz subspace)

!Binomial and Multinomial function, which are used to compute Hilbert space dimensions

!Some procedures which compute the imbalance, local imbalance of the spin basis

!A procedure which generates the full Hilbert-space state corresponding to
!   the Sz=0 subspace and generic Sz subspace 

!Other utilities: binary search algorithm, conversion from 1D grid to 2D grid

module functions
  
#ifdef __INTEL_COMPILER
    USE IFPORT
#endif
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  use printing
  implicit none

  complex (dcp), private, parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), private, parameter :: C_ONE = dcmplx(1._dp, 0._dp) 
  complex (dcp), private, parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  integer (ip), private, parameter :: dimSpin1 = 3
  real (dp), private, parameter :: tol = 1.0e-6
  real (dp), parameter, private :: pi = 4.d0 * datan(1.d0)
  !integer (ip) :: i, j, s
  !real (dp), parameter, private :: SigmaZ(3,3) = reshape((/ (( (2-i)*merge(1,0,i==j), j=1,3),i=1,3) /), shape(SigmaZ)) 

    !SigmaZ = 0
    !do i = 1, 3
    !  s = 2 - i
    !  SigmaZ(i,i) = s
    !enddo

contains

  subroutine decode(decimal, ntrits, tritstring)

    !Decode integer label for spin basis state |i> as a list of integers:
    ! i = sum_{k=1}^L a_k 3^{k-1} such that i + 1 = 1, ..., dim is an array index <->
    ! a_k = mod( floor(i/3^{k-1}) , 3) = 0, 1, 2 such that a_k + 1 = 1, 2, 3 is an array index (in tensor product notation, on a
    ! single Hilbert space)
    integer (ip), intent(in) :: decimal, ntrits
    integer (ip), intent(out) :: tritstring(ntrits)
    integer (ip) :: k, num

    tritstring = 0
    num = decimal
    do k = 1, ntrits
        tritstring(k) = mod(num,dimSpin1)
        num = num / dimSpin1
    enddo

  end subroutine decode

  function trit_to_spin(trit)

    !Correspondence between integer configuration a_k and spin configuration:
    ! s_k = 1 - a_k = 1, 0, -1 <-> a_k = 1 - s_k = 0, 1, 2
    integer (ip), intent(in) :: trit
    integer (ip) :: trit_to_spin
    
    trit_to_spin = 1 - trit

  end function

  function spin_to_trit(spin)

    !Correspondence between integer configuration a_k and spin configuration:
    ! s_k = 1 - a_k = 1, 0, -1 <-> a_k = 1 - s_k = 0, 1, 2
    integer (ip), intent(in) :: spin
    integer (ip) :: spin_to_trit
    
    spin_to_trit = 1 - spin

  end function

  subroutine init_random_seed()
    use iso_fortran_env, only: int64
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/random", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat) !urandom vs random
    if (istat == 1) then
      !print *, "reading from /dev/random"
      read(un) seed
      close(un)
    else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.

      !print *, "using clock time and PID"
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
#ifdef _OPENMP
      pid = pid + omp_get_thread_num() * 100
      !print *, "Working with openmp: adding thread number to pid to ensure independence &
      !  & even with thread race. "
      !print *, " 'New PID' to XOR with (should be different for each thread): ", pid
#endif
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    print *, "seed = ", seed(:)
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

    !integer function omp_get_thread_num
    !end function
  end subroutine init_random_seed


  integer function binom(n,k)
    integer(ip), intent(in) :: n,k
    integer, parameter :: i8 = selected_int_kind(18)
    integer, parameter :: dp = selected_real_kind(15)

    !if (k == n) then
    !    binom = 1
    !else if (k == 1) then
    !    binom = n
    !else if ((k /= 1) .and. (k /= n)) then
      binom = nint(exp(log_gamma(n+1.0_dp)-log_gamma(n-k+1.0_dp)-log_gamma(k+1.0_dp)),kind=i8)
    !end if 
  end function

  integer function multinom(n,k_vec)
    integer (ip), intent(in) :: n
    integer (ip), allocatable, intent(in) :: k_vec(:)
    

    integer, parameter :: i8 = selected_int_kind(18)
    integer, parameter :: dp = selected_real_kind(15)
    !integer :: i

    !print *, "k_vec = ", k_vec(:), "n = ", n
    !print *, "sum = ", sum(k_vec)
    if(sum(k_vec) /= n .OR. any(k_vec<0) ) stop "Error in multinom(n,k_vec): sum(k_vec) != n"

    multinom = nint(exp( log_gamma(n+1.0_dp)-sum(log_gamma(k_vec+1.0_dp)) ),kind=i8)

  end function


!
!
!
!
!
!
!  integer function non_zero_HMBL_Sz0(L)
!
!    integer, intent(in) :: L
!    integer :: nz
!
!    if (L==2) then
!      nz = 4
!    else if (L==4) then
!      nz = 18
!    else if (L==6) then
!      nz = 80
!    else if (L==8) then
!      nz = 350
!    else if (L==10) then
!      nz = 1512
!    else if (L==12) then
!      nz = 6468
!    else if (L==14) then
!      nz = 27456
!    else if (L==16) then
!      nz = 115830
!    else if (L==18) then
!      nz = 486200
!    else if (L==20) then 
!      nz = 2032316
!    else if (L==22) then
!      nz = 8465184
!    else
!      print *, "Chain Length is odd or too large"
!      nz = 0
!    endif
!    
!    non_zero_HMBL_Sz0 = nz
!
!  end function non_zero_HMBL_Sz0
!
!
!
  function dimSpin1_Sz0(nspin)

    integer (ip), intent(in) :: nspin
    integer (ip) :: dimSpin1_Sz0
    integer (ip) :: n_up, ntot
    integer (ip), allocatable :: k_vec(:)

    allocate(k_vec(dimSpin1))

    !In S_z = 0 configurations, N_up = N_down with the constraint N_up + N_down + N_zero = nspin
    ntot = 0
    do n_up = 0, nspin/2
      k_vec(1) = n_up
      k_vec(2) = n_up
      k_vec(3) = nspin - 2*n_up
      ntot = ntot + multinom(nspin,k_vec)
    enddo

    dimSpin1_Sz0 = ntot

  end function

  function dimSpin1_Sz(nspin, Sz)

    integer (ip), intent(in) :: nspin, Sz
    integer (ip) :: dimSpin1_Sz
    integer (ip) :: n_up, ntot
    integer (ip), allocatable :: k_vec(:)

    allocate(k_vec(dimSpin1))

    !In S_z = 0 configurations, N_up = N_down with the constraint N_up + N_down + N_zero = nspin
    if( Sz < -nspin .OR. Sz > nspin ) then
      print *, "Error: The value of Sz is invalid. Sz = ", Sz
      stop
    endif
    ntot = 0
    do n_up = max(0,Sz), (nspin+Sz)/2
      k_vec(1) = n_up
      k_vec(2) = n_up - Sz
      k_vec(3) = nspin + Sz - 2*n_up
      ntot = ntot + multinom(nspin,k_vec)
    enddo

    dimSpin1_Sz = ntot

  end function



  subroutine zero_mag_states(nspin, dim_Sz0, states)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (ip), intent(in) :: nspin, dim_Sz0
    integer (ip), intent(out) :: states(dim_Sz0)

    integer :: config(nspin)
    integer (ip) :: i, k, dim

    dim = dimSpin1**nspin
    k = 0
    states = 0
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(1-config)==0) then
        k = k+1
        states(k) = i
        !print "(2(I4,4X),*(I0))", i, k, config(:)
      endif
    enddo

  end subroutine zero_mag_states

  subroutine zero_mag_states_inv(nspin, dim_Sz0, states, inverse)
 
     !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
     integer (ip), intent(in) :: nspin, dim_Sz0
     integer (ip), intent(out) :: states(dim_Sz0), inverse(dimSpin1**nspin)
 
     integer :: config(nspin)
     integer (ip) :: i, l, dim
 
     dim = dimSpin1**nspin
     l = 0
     states = 0
     inverse = -1
 
     !print "(2(A4,4X),A4)", "i", "l", "conf"
     do i = 0, dim-1
 
       call decode(i, nspin, config)
 
       if (sum(1-config)==0) then
         l = l+1
         states(l) = i
         inverse(i+1) = l
         !print "(2(I4,4X),*(I0))", i, l, config(:)
       endif
     enddo
     !print *, ""
 
   end subroutine

  subroutine basis_Sz(nspin, dim_Sz, Sz, states)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    integer (ip), intent(out) :: states(dim_Sz)

    integer :: config(nspin)
    integer (ip) :: i, k, dim

    dim = dimSpin1**nspin
    k = 0
    states = 0
    !print "(2(A4,4X),A4)", "k", "i", "conf"
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(1-config)==Sz) then
        k = k+1
        states(k) = i
        !print "(2(I4,4X),*(I0))", k, i, config(:)
      endif
    enddo

  end subroutine

  subroutine basis_Sz_inv(nspin, dim_Sz, Sz, states, inverse)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    integer (ip), intent(out) :: states(dim_Sz), inverse(dimSpin1**nspin)

    integer :: config(nspin)
    integer (ip) :: i, l, dim

    dim = dimSpin1**nspin
    l = 0
    states = 0
    inverse = -1
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(1-config)==Sz) then
        l = l+1
        states(l) = i
        inverse(i+1) = l
        !print "(2(I4,4X),*(I0))", i, l, config(:)
      endif
    enddo

  end subroutine

  subroutine buildState_Sz0_to_FullHS(nspin, dim, dim_Sz0, psi_Sz0, psi)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (ip), intent(in) :: nspin, dim, dim_Sz0
    complex (dcp), intent(in) :: psi_Sz0(dim_Sz0)
    complex (dcp), intent(out) :: psi(dim)

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

  subroutine buildState_Sz_to_FullHS(nspin, dim, dim_Sz, Sz, psi_Sz, psi)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (ip), intent(in) :: nspin, dim, dim_Sz, Sz
    complex (dcp), intent(in) :: psi_Sz(dim_Sz)
    complex (dcp), intent(out) :: psi(dim)

    integer :: i, l, states(dim_Sz)

    call basis_Sz(nspin, dim_Sz, Sz, states)

    psi = 0
    !print *, "psi initialized"
    do l = 1, dim_Sz
      i = states(l) + 1
      psi(i) = psi_Sz(l)
      !print *,l, psi(l)
      !print *, l, dim, i, dim_Sz0
    enddo

  end subroutine buildState_Sz_to_FullHS

  subroutine projectState_FullHS_to_Sz(nspin, dim, psi, dim_Sz, Sz, psi_Sz)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (ip), intent(in) :: nspin, dim
    complex (dcp), intent(in) :: psi(dim)
    complex (dcp), allocatable, intent(out) :: psi_Sz(:)
    integer (ip), intent(out) :: dim_Sz, Sz

    integer :: i, l, config(nspin)
    real :: mag_psi, mag_s
    integer, allocatable :: states(:)
    logical :: flag


    !--------- Check Sz ----------!
    mag_psi = 0
    do i = 1, dim
      call decode(i-1, nspin, config)
      mag_s = sum(1-config)
      mag_psi = mag_psi + abs(psi(i))**2 * mag_s
    enddo
    print *, "Magnetization of State:", mag_psi


    flag = .True.
    do i = 1, dim
      call decode(i-1, nspin, config)
      mag_s = sum(1-config)
      if ( abs( (mag_s - mag_psi) * psi(i) ) > tol ) flag = .False.
    enddo


    Sz = int(mag_psi)
    dim_Sz = dimSpin1_Sz(nspin, Sz)
    if (flag) then
      print *, "The State has a definite magnetization Sz = ", Sz
    else
      stop "The State is not an eigenstate of Sz"
    endif

    allocate(psi_Sz(dim_Sz), states(dim_Sz))

    psi_Sz = 0
    call basis_Sz(nspin, dim_Sz, Sz, states)
    !print *, "psi initialized"
    do l = 1, dim_Sz
      i = states(l) + 1
      psi_Sz(l) = psi(i)
      !print *,l, psi(l)
      !print *, l, dim, i, dim_Sz0
    enddo

    if(abs( dot_product(psi_Sz,psi_Sz) - dot_product(psi,psi) ) > tol) & 
      & stop "Error with projection of state to Sz subspace"

  end subroutine


  integer (ip) function binsearch(val, array)
  
  
    implicit none
    integer (ip), intent(in) :: val, array(:)
    integer (ip) :: mid, start, finish, range
    
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

  function binsearch_closest_in_circle(theta, array) result(indx)

    !Binary search for closest element to given value where the array contains angles in [-pi, pi]

    implicit none
    real (dp), intent(in) :: theta, array(:)
    integer (ip) :: indx, mid, left, right!, range
    integer (ip) :: n!, prev, next

    indx = -1
    left = 1
    right = size(array) + 1
    n = size(array)

    !range = right - left

    if (theta < array(1)) then
      indx = merge(1, size(array), abs(array(1)-theta)<= abs( (theta+2*pi)-array(size(array)) ) )
      !print *, "indx(1) = ", indx
      return
    else if (theta > array(size(array))) then
      indx = merge(size(array), 1, abs( theta - array(size(array)) )<= abs(array(1)-(theta-2*pi)) )
      !print *, "indx(n) = ", indx
      return
    endif

    !print "(6X,4(A12),2(A26))", "indx", "left", "right", "mid", "array(mid)", "val"

    mid = (right + left) / 2
    do while (right > left)

      !prev = merge(mid-1,n, mid.ne.1)
      !next = merge(mid+1,1, mid.ne.n)
      !print *, "before", indx, left, right, mid, array(mid), theta
      !print *, array(prev), array(mid), array(next)

      if (array(mid) == theta) then
        indx = mid
        !print *, "indx = ", indx
        return
      else if (array(mid) > theta) then
        if (mid>1 .and. array(mid-1) < theta) then
          indx = merge(mid, mid-1, theta-array(mid-1) >= array(mid)-theta)
          !print *, "indx(mid-1:mid) = ", indx
          return
        endif
        right = mid-1
      else ! array(mid) < val) 
        if (mid<size(array) .and. array(mid+1) > theta) then
          indx = merge(mid, mid+1, array(mid+1)-theta >= theta-array(mid))
          !print *, "indx(mid:mid+1) = ", indx
          return
        endif
        left = mid+1
      end if
      !range = right - left
      mid = (right + left) / 2
      
      !print *, "after ", indx, left, right, mid, array(mid), theta
      !print *, array(prev), array(mid), array(next)

    end do

    !indx = right! - 1 
    indx = mid
    !print *, "indx(end) = ", indx
 
  end function

  subroutine disorder_average(local_avg, local_sq, global_avg, global_sigma)
    
    real (dp), intent(in) :: local_avg(:), local_sq(:)
    real (dp), intent(out) :: global_avg, global_sigma
  
    real (dp) :: global_sq
    integer (ip) :: n
  
    n = size(local_avg)
  
    global_avg = sum(local_avg) / n
    global_sq = sum(local_sq) / n
    global_sigma = sqrt( (global_sq - global_avg**2) / n )
  
  end subroutine

  function normalization_power_law(alpha, nspin) result(norm)

    real (dp) :: alpha, norm
    integer (ip) :: nspin

    if (alpha > 1) then
      norm = 1
    else if (abs(alpha-1)<tol) then
      norm = log(real(nspin, kind=dp))
    else
      norm = nspin**(1-alpha)
    endif

  end function

end module functions
