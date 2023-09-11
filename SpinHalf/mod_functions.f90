!This module contains various functions and procedures used throughout the project

!Specifically it includes all decoding/encoding procedures to go 
!   from an integer to a basis (full Hilbert space, Sz=0 subspace, generic Sz)

!Binomial function, which is used to compute Hilbert space dimensions

!Some procedures which compute the imbalance, local imbalance of the spin basis

!A procedure which generates the full Hilbert-space state corresponding to
!   the Sz=0 subspace and generic Sz subspace

!Other utilities: binary search algorithm(s), conversion from 1D grid to 2D grid

module functions
  
  !use ifport
#ifdef __INTEL_COMPILER
    USE IFPORT
#endif

  use iso_c_binding, only: dp => c_double, ip => c_int, cp => c_double_complex
  use printing
  implicit none

  complex (cp), private, parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (cp), private, parameter :: C_ONE = dcmplx(1._dp, 0._dp) 
  complex (cp), private, parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  integer (ip), private, parameter :: dimSpinHalf = 2
  real (dp), parameter, private :: pi = 4._dp * atan(1._dp)
  real (dp), private, parameter :: tol = EPSILON(1._dp)
  !integer (ip) :: i, j, s
  !real (dp), parameter, private :: SigmaZ(3,3) = reshape((/ (( (2-i)*merge(1,0,i==j), j=1,3),i=1,3) /), shape(SigmaZ)) 

    !SigmaZ = 0
    !do i = 1, 3
    !  s = 2 - i
    !  SigmaZ(i,i) = s
    !enddo

contains

  subroutine decode(decimal, nbits, bitstring)

    integer (ip), intent(in) :: decimal, nbits
    integer (ip), intent(out) :: bitstring(nbits)
    integer (ip) :: i

    bitstring = 0
    do i = 1, nbits
      if (btest(decimal, i - 1)) bitstring(i) = 1
    enddo
  end subroutine decode

  function swap_config(nspin, i) result(j)

    integer (ip), intent(in) :: nspin, i
    integer (ip) :: config(nspin), k, j

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


  integer (ip) function binom(n,k)
    use iso_c_binding, only: lp => c_float128
    integer (ip), intent(in) :: n,k

    if (k == n) then
        binom = 1
    else if (k == 1) then
        binom = n
    else if ((k /= 1) .and. (k /= n)) then
      binom = nint(exp(log_gamma(n+1.0_lp)-log_gamma(n-k+1.0_lp)-log_gamma(k+1.0_lp)),kind=ip)
    end if 
    !print *, "c_float_result binom =", exp(log_gamma(n+1.0_lp)-log_gamma(n-k+1.0_lp)-log_gamma(k+1.0_lp))
    !print *, "integer result = ", binom
  end function




  subroutine buildState_Sz0_to_FullHS(nspin, dim, dim_Sz0, psi_Sz0, psi)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (ip), intent(in) :: nspin, dim, dim_Sz0
    complex (cp), intent(in) :: psi_Sz0(dim_Sz0)
    complex (cp), intent(out) :: psi(dim)

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
    integer (ip), intent(in) :: nspin, dim_Sz0
    integer (ip), intent(out) :: states(dim_Sz0)

    integer :: config(nspin)
    integer (ip) :: i, k

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
    integer (ip), intent(in) :: nspin, dim_Sz0
    integer (ip), intent(out) :: states(dim_Sz0), inverse(dimSpinHalf**nspin)

    integer :: config(nspin), spin(nspin)
    integer (ip) :: i, l, dim

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





  function binsearch(val, array) result(indx)
  
  
    implicit none
    integer (ip), intent(in) :: val, array(:)
    integer (ip) :: indx
    integer (ip) :: mid, start, finish, range
    
    indx = -1
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
      indx = mid
    end if
  
  end function

  function binsearch_closest_from_above(val, array) result(indx)

    implicit none
    real (dp), intent(in) :: val, array(:)
    integer (ip) :: indx, mid, left, right, range

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
    !print *, "indx = ", indx
 
  end function

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


  subroutine find_degeneracies( n, E, idxuE, deg )
    !Finds all degeneracies of a real ordered 1D array 'E' of length n (up to a fixed tolerance)
    !idxuE contains the indices of the unique values of E
    !deg contains the corresponding degeneracy
  
    use iso_c_binding
    implicit none
  
    integer (ip), intent(in) :: n
    real (dp), intent(in) :: E(n) !'E' has to be sorted already
    !real (dp), intent(out)  :: idxuE(n)
    integer (ip), intent(out)  :: deg(n), idxuE(n)
  
    integer :: i, j
    real (dp) :: tol = EPSILON(1.0_dp) ! 1.0e-10

    idxuE = 0
    deg = 1
  
    j = 1
    idxuE(j) = 1 !E(j)
    do i = 2, n
      if ( abs(E(i) - E(i-1)) > tol ) then
        j = j + 1
        !idxuE(j) = E(i)
        idxuE(j) = i
      else
        deg(j) = deg(j) + 1
      endif
    enddo
    deg(j+1:n) = 0
  
    !print *, "Degenerate E: "
    !print "(*(A15))", "i", "idxunique(i)", "unique(idx)", "deg(i)"
    !do i = 1, n
    !  if (deg(i)>0) then
    !    !print *, i, idxuE(i), E(idxuE(i)), deg(i)
    !  endif
    !enddo
  
  end subroutine

  subroutine printmat_as_list_C(nspin, dim, M, t)

    integer, intent(in) :: nspin, dim
    complex(cp), intent(in), dimension(dim,dim) :: M
    character, intent(in) :: t*1

    integer :: i, j, config1(nspin), config2(nspin)

    if (t == 'C') then
      do i = 1, dim
        do j = 1, dim
          call decode(i-1, nspin, config1)
          call decode(j-1, nspin, config2)
          print "(*(I0))", config1
          write (*, 97) M(i,j), config2
          97 format('|',(1(sf26.20spf26.20x'i':x),*(I0)))
        enddo
      enddo
    else if (t == 'R') then
      do i = 1, dim
        do j = 1, dim
          print *, i, j, real(M(i,j), kind=kind(M))
        enddo
      enddo
    else if (t == 'I') then
      do i = 1, dim
        do j = 1, dim
          print *, i, j, int(M(i,j))
        enddo
      enddo
    else if (t == 'A') then
      do i = 1, dim
        do j = 1, dim
          print *, i, j, abs(M(i,j))
        enddo
      enddo
    end if

  end subroutine


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

  function dimSpinHalf_Sz(nspin, Sz) result(dim_Sz)

    integer (ip), intent(in) :: nspin, Sz
    integer (ip) :: dim_Sz

    dim_Sz = binom(nspin, (nspin + Sz)/2)

  end function

  subroutine basis_Sz(nspin, dim_Sz, Sz, states)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    integer (ip), intent(out) :: states(dim_Sz)

    integer :: config(nspin)
    integer (ip) :: i, l, dim

    dim = dimSpinHalf**nspin
    l = 0
    states = 0
    print *, "dim_Sz = ", dim_Sz, "binom(L, L/2) =", binom(nspin, nspin/2)
    !print "(2(A4,4X),A4)", "k", "i", "conf"
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(1-2*config)==Sz) then
        l = l+1
        if(l>dim_Sz) then
          print*, "Error:"
          print *, config(:), l
        endif
        states(l) = i
        !print "(2(I4,4X),*(I0))", k, i, config(:)
      endif
    enddo

  end subroutine

  subroutine basis_Sz_inv(nspin, dim_Sz, Sz, states, inverse)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    integer (ip), intent(out) :: states(dim_Sz), inverse(dimSpinHalf**nspin)

    integer :: config(nspin)
    integer (ip) :: i, l, dim

    dim = dimSpinHalf**nspin
    l = 0
    states = 0
    inverse = -1
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(1-2*config)==Sz) then
        l = l+1
        states(l) = i
        inverse(i+1) = l
        !print "(2(I4,4X),*(I0))", i, l, config(:)
      endif
    enddo

  end subroutine

  subroutine buildState_Sz_to_FullHS(nspin, dim, dim_Sz, Sz, psi_Sz, psi)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (ip), intent(in) :: nspin, dim, dim_Sz, Sz
    complex (cp), intent(in) :: psi_Sz(dim_Sz)
    complex (cp), intent(out) :: psi(dim)

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
    complex (cp), intent(in) :: psi(dim)
    complex (cp), allocatable, intent(out) :: psi_Sz(:)
    integer (ip), intent(out) :: dim_Sz, Sz

    integer :: i, l, config(nspin)
    real (dp) :: mag_psi, mag_s
    integer, allocatable :: states(:)
    logical :: flag


    !--------- Check Sz ----------!
    mag_psi = 0
    do i = 1, dim
      call decode(i-1, nspin, config)
      mag_s = sum(1-2*config)
      mag_psi = mag_psi + abs(psi(i))**2 * mag_s
    enddo
    print *, "Magnetization of State:", mag_psi


    flag = .True.
    do i = 1, dim
      call decode(i-1, nspin, config)
      mag_s = sum(1-2*config)
      if ( abs( (mag_s - mag_psi) * psi(i) ) > tol ) flag = .False.
    enddo


    Sz = nint(mag_psi)
    dim_Sz = dimSpinHalf_Sz(nspin, Sz)
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
      i = states(l)
      !call decode(i, nspin, config)
      psi_Sz(l) = psi(i+1)
      !print *, l, psi(l)
      !print *, l, dim, i, dim_Sz, psi(l)
      !print "(*(I0))", config(:)
    enddo

    if(abs( dot_product(psi_Sz,psi_Sz) - dot_product(psi,psi) ) > tol) & 
      & stop "Error with projection of state to Sz subspace"

  end subroutine


end module functions
