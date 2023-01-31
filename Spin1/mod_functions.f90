!This module contains various functions and procedures used throughout the project

!Specifically it includes all decoding/encoding procedures to go 
!   from an integer to a basis (full Hilbert space, Sz=0 subspace, generic Sz subspace)

!It also includes a simple binary search algorithm

!Binomial and Multinomial function, which are used to compute Hilbert space dimensions

!Some procedures which compute the imbalance, local imbalance of the spin basis

!A procedure which generates the full Hilbert space vector corresponding to
!   the Sz=0 subspace and generic Sz subspace 


module functions
  
  !use ifport
  use iso_c_binding
  use printing
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  integer (c_int), private, parameter :: dimSpin1 = 3
  real (c_double), private, parameter :: tol = 1.0e-6
  !integer (c_int) :: i, j, s
  !real (c_double), parameter, private :: SigmaZ(3,3) = reshape((/ (( (2-i)*merge(1,0,i==j), j=1,3),i=1,3) /), shape(SigmaZ)) 

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
    integer (c_int), intent(in) :: decimal, ntrits
    integer (c_int), intent(out) :: tritstring(ntrits)
    integer (c_int) :: k, num

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
    integer (c_int), intent(in) :: trit
    integer (c_int) :: trit_to_spin
    
    trit_to_spin = 1 - trit

  end function

  function spin_to_trit(spin)

    !Correspondence between integer configuration a_k and spin configuration:
    ! s_k = 1 - a_k = 1, 0, -1 <-> a_k = 1 - s_k = 0, 1, 2
    integer (c_int), intent(in) :: spin
    integer (c_int) :: spin_to_trit
    
    spin_to_trit = 1 - spin

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

    !if (k == n) then
    !    binom = 1
    !else if (k == 1) then
    !    binom = n
    !else if ((k /= 1) .and. (k /= n)) then
      binom = nint(exp(log_gamma(n+1.0_dp)-log_gamma(n-k+1.0_dp)-log_gamma(k+1.0_dp)),kind=i8)
    !end if 
  end function

  integer function multinom(n,k_vec)
    integer (c_int), intent(in) :: n
    integer (c_int), allocatable, intent(in) :: k_vec(:)
    

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

    integer (c_int), intent(in) :: nspin
    integer (c_int) :: dimSpin1_Sz0
    integer (c_int) :: n_up, ntot
    integer (c_int), allocatable :: k_vec(:)

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

    integer (c_int), intent(in) :: nspin, Sz
    integer (c_int) :: dimSpin1_Sz
    integer (c_int) :: n_up, ntot
    integer (c_int), allocatable :: k_vec(:)

    allocate(k_vec(dimSpin1))

    !In S_z = 0 configurations, N_up = N_down with the constraint N_up + N_down + N_zero = nspin
    if( Sz < -nspin .OR. Sz > nspin ) stop "Error: The value of Sz is invalid."
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
    integer (c_int), intent(in) :: nspin, dim_Sz0
    integer (c_int), intent(out) :: states(dim_Sz0)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, dim

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
     integer (c_int), intent(in) :: nspin, dim_Sz0
     integer (c_int), intent(out) :: states(dim_Sz0), inverse(dimSpin1**nspin)
 
     integer :: config(nspin)
     integer (c_int) :: i, j, l, dim
 
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
    integer (c_int), intent(in) :: nspin, dim_Sz, Sz
    integer (c_int), intent(out) :: states(dim_Sz)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, dim

    dim = dimSpin1**nspin
    k = 0
    states = 0
    print "(2(A4,4X),A4)", "k", "i", "conf"
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(1-config)==Sz) then
        k = k+1
        states(k) = i
        print "(2(I4,4X),*(I0))", k, i, config(:)
      endif
    enddo

  end subroutine

  subroutine basis_Sz_inv(nspin, dim_Sz, Sz, states, inverse)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (c_int), intent(in) :: nspin, dim_Sz, Sz
    integer (c_int), intent(out) :: states(dim_Sz), inverse(dimSpin1**nspin)

    integer :: config(nspin)
    integer (c_int) :: i, j, l, dim

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
!
!  subroutine finite_imbalance_states(nspin, dim, IMB, states)
!    integer (c_int), intent(in) :: nspin, dim
!    real (c_double), intent(in) :: IMB
!    integer (c_int), intent(out) :: states(dim)
!
!    integer (c_int) :: i, k, l, config(nspin), sum_even, sum_odd
!    real (c_double) :: loc_imb
!
!    l = 0
!    print *, "       i     l  config"
!    do i = 0, dim-1
!      call decode(i, nspin, config)
!      sum_odd = 0
!      sum_even = 0
!      !print *, i, sum_even, sum_odd
!      do k = 1, nspin/2
!        sum_odd = sum_odd + config(2*k-1)
!        sum_even = sum_even + config(2*k)
!        !print *, k, sum_even, sum_odd
!      enddo
!      loc_imb = real(sum_odd - sum_even, c_double) / real(nspin - sum_odd - sum_even, c_double)
!      !print *, loc_imb
!      if(abs(IMB) < tol .AND. sum_odd == sum_even) then
!        l = l+1
!        states(l) = i
!        print "(2(8X,I0),4X,*(I0))", i, l, config(:)
!      else if(abs(IMB) > tol .AND. abs(loc_imb - IMB) < tol) then
!        l = l+1
!        states(l) = i
!        print "(2(8X,I0),4X,*(I0))", i, l, config(:)
!      endif
!    enddo
!
!  end subroutine finite_imbalance_states
!
!  subroutine large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, states)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: IMB, LI
!    integer (c_int), intent(out) :: states(dim_Sz0)
!
!    integer :: i, l, m, idxSz0(dim_Sz0)
!    real (c_double) :: IMBc, LIc
!
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    m = 0
!    states = 0
!    !print "( 4X,A4, 4X,A4, 4X,A6, 4X,A6 )", "m", "l", "IMB", "LI"
!    do i = 1, dim_Sz0
!
!      l = idxSz0(i)
!      IMBc = imbalance_basis(nspin, l)
!      LIc = local_imbalance_basis(nspin, l)
!      if (IMB >= 0) then
!        if (IMBc >= IMB .AND. LIc >= LI) then
!          m = m+1
!          states(m) = l
!          !print "(4X,I4, 4X,I4, 4X,F6.3, 4X,F6.3, 4X,F6.3)", m, l, IMBc, LIc
!        endif
!      else if (IMB < 0) then
!        if (IMBc <= IMB .AND. LIc >= LI) then
!          m = m+1
!          states(m) = l
!          !print "(4X,I4, 4X,I4, 4X,F6.3, 4X,F6.3)", m, l, IMBc, LIc
!        endif
!      endif
!    enddo
!
!    !print *, "states = ", states(1:m)
!
!  end subroutine
!
!  subroutine finite_imbalance_states_Sz0(nspin, dim_Sz0, IMB, states)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: IMB
!    integer (c_int), intent(out) :: states(dim_Sz0)
!
!    integer (c_int) :: i, k, m, l, config(nspin), sum_even, sum_odd, idxSz0(dim_Sz0)
!    real (c_double) :: loc_imb
!
!    m = 0
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    print *, "       i     l  config"
!    do i = 1, dim_Sz0
!      l = idxSz0(i)
!      call decode(l, nspin, config)
!      sum_odd = 0
!      sum_even = 0
!      !print *, i, sum_even, sum_odd
!      do k = 1, nspin/2
!        sum_odd = sum_odd + config(2*k-1)
!        sum_even = sum_even + config(2*k)
!        !print *, k, sum_even, sum_odd
!      enddo
!      loc_imb = real(sum_odd - sum_even, c_double) / real(nspin - sum_odd - sum_even, c_double)
!      !print *, loc_imb
!      if(abs(IMB) < tol .AND. sum_odd == sum_even) then
!        m = m+1
!        states(m) = l
!        print "(2(8X,I0),4X,*(I0))", l, m, config(:)
!      else if(abs(IMB) > tol .AND. abs(loc_imb - IMB) < tol) then
!        m = m+1
!        states(m) = l
!        print "(2(8X,I0),4X,*(I0))", l, m, config(:)
!      endif
!    enddo
!
!  end subroutine
!
!
!
!
!
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

  subroutine buildState_Sz_to_FullHS(nspin, dim, dim_Sz, Sz, psi_Sz, psi)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (c_int), intent(in) :: nspin, dim, dim_Sz, Sz
    complex (c_double_complex), intent(in) :: psi_Sz(dim_Sz)
    complex (c_double_complex), intent(out) :: psi(dim)

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

  subroutine projectState_FullHS_to_Sz(nspin, dim, psi, Sz, psi_Sz)
    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
    !simply by constructing a vector which has zero components outside the Sz=0 subspace

    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), allocatable, intent(out) :: psi_Sz(:)
    integer (c_int), intent(out) :: Sz

    integer :: i, l, dim_Sz, sigmaz(nspin), config(nspin)
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
      psi_Sz(l) = psi_Sz(i)
      !print *,l, psi(l)
      !print *, l, dim, i, dim_Sz0
    enddo

    if(abs( dot_product(psi_Sz,psi_Sz) - dot_product(psi,psi) ) > tol) & 
      & stop "Error with projection of state to Sz subspace"

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




end module functions
