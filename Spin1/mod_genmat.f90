module genmat
  
  !use ifport
  use iso_c_binding
  use printing
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  integer (c_int), private, parameter :: dimSpin1 = 3
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

    if (k == n) then
        binom = 1
    else if (k == 1) then
        binom = n
    else if ((k /= 1) .and. (k /= n)) then
      binom = nint(exp(log_gamma(n+1.0_dp)-log_gamma(n-k+1.0_dp)-log_gamma(k+1.0_dp)),kind=i8)
    end if 
  end function

  function KDelta(i, j)

    integer (c_int), intent(in) :: i, j
    integer (c_int) :: KDelta 

    KDelta = merge(1,0,i==j)

  end function

  function identity(N)

    integer (c_int), intent(in) :: N
    real (c_double) :: identity(N,N)

    integer :: i, j

    do i = 1, N
      do j = 1, N
        identity(i,j) = KDelta(i,j)
      enddo
    enddo
    
  end function

  

!  function SigmaZ()
!
!    real (c_double) :: SigmaZ(3,3)
!    integer (c_int) :: i, s
!
!    SigmaZ = 0
!    do i = 1, 3
!      s = 2 - i
!      SigmaZ(i,i) = s
!    enddo
!
!  end function
!
!  function SigmaX()
!
!    real (c_double) :: SigmaX(3,3)
!    integer (c_int) :: i, j, s, sp
!
!    SigmaX = 0
!    do i = 1, 3
!      do j = 1, 3
!        s = 2 - i
!        sp = 2 - j
!        SigmaX(j,i) = (KDelta(sp,s+1) + KDelta(sp,s-1) ) * sqrt(2.d0)
!      enddo
!    enddo
!
!  end function
!
!  function SigmaY()
!
!    complex (c_double_complex) :: SigmaY(3,3)
!    integer (c_int) :: i, j, s, sp
!
!    SigmaY = 0
!    do i = 1, 3
!      do j = 1, 3
!        s = 2 - i
!        sp = 2 - j
!        SigmaY(j,i) = (KDelta(sp,s+1) - KDelta(sp,s-1) ) / ( sqrt(2.d0) * C_UNIT )
!      enddo
!    enddo
!
!  end function
!
!  function SigmaZZ()
!
!    real (c_double) :: SigmaZZ(3,3,3,3)
!    integer (c_int) :: i1, i2, s1, s2
!
!    SigmaZZ = 0
!    do i1 = 1, 3
!      do i2 = 1, 3
!        s1 = 2 - i1
!        s2 = 2 - i2
!        SigmaZZ(i1,i2,i1,i2) = s1 * s2
!      enddo
!    enddo
!
!  end function
!
!  function SigmaXY()
!
!    real (c_double) :: SigmaXY(3,3,3,3)
!    integer (c_int) :: i1, i2, j1, j2, s1, s2, sp1, sp2
!
!    SigmaXY = 0
!    do i1 = 1, 3
!      do i2 = 1, 3
!        do j1 = 1, 3
!          do j2 = 1, 3
!            s1 = 2 - i1
!            s2 = 2 - i2
!            sp1 = 2 - j1
!            sp2 = 2 - j2
!            SigmaXY(j1,j2, i1,i2) = kdelta(sp1, s1 + 1) * kdelta(sp2, s2 - 1) * kdelta(sp1, s1 - 1) * kdelta(sp2, s2 + 1)
!          enddo
!        enddo
!      enddo
!    enddo
!
!  end function

  subroutine buildOneBody(nspin, dim, h, HOneBody)

    !Example of OneBody operator of the form:
    !HOneBody = sum_{k=1}^L h_k OPR_k  where OPR_k is a one-body operator acting on spin 'k'
    ! (defined in the function, although it might be an argument)
    ! In this case OPR = sigma_k^z is a spin-1 operator
    
    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: h(nspin)
    real (c_double), intent(out) :: HOneBody(dim,dim)
    real (c_double) :: OPR(dimSpin1,dimSpin1)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, idx, idxp
    
    do idx = 1, dimSpin1
      do idxp = 1, dimSpin1
        OPR(idxp, idx) = (2-idx)*KDelta(idxp, idx)
      enddo
    enddo

    !print *, "OPR = "
    !call printmat(dimSpin1, OPR, 'R')

    HOneBody = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin

        do idxp = 0, dimSpin1-1

          idx = config(k)
          j = i + (idxp - idx) * dimSpin1**(k-1)
          HOneBody(j+1,i+1) = HOneBody(j+1,i+1) + h(k) * OPR(idxp+1, idx+1)

        enddo

      enddo
    enddo

    !print*, "h = "
    !print "(*(F5.1))", h(:)

    !print*, "H1 = "
    !call printmat(dim,HOneBody, 'R')

  end subroutine buildOneBody

  subroutine buildTwoBody(nspin, dim, Jint, HTwoBody)

    !Example of TwoBody operator of the form:
    !HOneBody = sum_{k=1}^L J_k OPR_{k,k+1}  
    !   where OPR_{k,k+1} is a two-body operator acting on spins 'k' and 'k+1'
    ! (defined in the subroutine, although it might be an argument)
    ! In this case OPR = sigma_k^x sigma_{k+1}^x + sigma_k^y sigma_{k+1}^y =
    !                  = (sigma_k^+ sigma_{k+1}^- + sigma_k^- sigma_{k+1}^+)/2
    ! In the sigma^z basis:
    ! 
    
    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1)
    real (c_double), intent(out) :: HTwoBody(dim,dim)
    real (c_double) :: OPR(dimSpin1,dimSpin1,dimSpin1,dimSpin1)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, idx1, idx2, idxp1, idxp2
    
    do idx1 = 1, dimSpin1
      do idx2 = 1, dimSpin1
        do idxp1 = 1, dimSpin1
          do idxp2 = 1, dimSpin1
            OPR(idxp1, idxp2, idx1, idx2) = ( KDelta(idxp1, idx1+1) * KDelta(idxp2,idx2-1) + &
            &  KDelta(idxp1, idx1-1) * KDelta(idxp2, idx2+1) )
            !print*, OPR(idxp1, idxp2, idx1, idx2), idxp1, idxp2, idx1, idx2
          enddo
        enddo
      enddo
    enddo

!    print *, "OPR = "
!    call printmat(dimSpin1, OPR, 'R')

    HTwoBody = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin-1

        do idxp1 = 0, dimSpin1-1
          do idxp2 = 0, dimSpin1-1

            idx1 = config(k)
            idx2 = config(k+1)
            j = i + (idxp1 - idx1) * dimSpin1**(k-1) + (idxp2 - idx2) * dimSpin1**k
            HTwoBody(j+1,i+1) = HTwoBody(j+1,i+1) + Jint(k) * OPR(idxp1+1, idxp2+1, idx1+1, idx2+1)

          enddo
        enddo

      enddo
    enddo

    !print*, "J = "
    !print "(*(F5.1))", Jint(:)

    !print*, "H2 = "
    !call printmat(dim,HTwoBody, 'R')

  end subroutine buildTwoBody

  subroutine buildTwoBodyZZ(nspin, dim, Vint, HTwoBody)

    !Example of TwoBody operator of the form:
    !HOneBody = sum_{k=1}^L V_k OPR_{k,k+1}  
    !   where OPR_{k,k+1} is a two-body operator acting on spins 'k' and 'k+1'
    ! (defined in the subroutine, although it might be an argument)
    ! In this case OPR = sigma_k^z sigma_{k+1}^z

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Vint(nspin-1)
    real (c_double), intent(out) :: HTwoBody(dim,dim)
    real (c_double) :: OPR(dimSpin1,dimSpin1,dimSpin1,dimSpin1)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, idx1, idx2, idxp1, idxp2, s1, s2
    
    do idx1 = 1, dimSpin1
      do idx2 = 1, dimSpin1
        do idxp1 = 1, dimSpin1
          do idxp2 = 1, dimSpin1
            s1 = 2 - idx1
            s2 = 2 - idx2
            OPR(idxp1, idxp2, idx1, idx2) = s1 * KDelta(idxp1, idx1) * s2 * KDelta(idxp2,idx2)
          enddo
        enddo
      enddo
    enddo

!    print *, "OPR = "
!    call printmat(dimSpin1, OPR, 'R')

    HTwoBody = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin-1

        do idxp1 = 0, dimSpin1-1
          do idxp2 = 0, dimSpin1-1

            idx1 = config(k)
            idx2 = config(k+1)
            j = i + (idxp1 - idx1) * dimSpin1**(k-1) + (idxp2 - idx2) * dimSpin1**k
            HTwoBody(j+1,i+1) = HTwoBody(j+1,i+1) + Vint(k) * OPR(idxp1+1, idxp2+1, idx1+1, idx2+1)

          enddo
        enddo

      enddo
    enddo

    !print*, "V = "
    !print "(*(F5.1))", Vint(:)

    !print*, "H2 = "
    !call printmat(dim,HTwoBody, 'R')

  end subroutine buildTwoBodyZZ

  subroutine buildHSwap(nspin, dim, HSwap)

    !Hamiltonian for Swap Operator U_Swap = e^{-i pi/2 * HSwap):
    !HOneBody = sum_{k=1}^{L/2} OPR_{k,k+1}
    !   where OPR = OPR2 + OPR1,
    !   OPR1 = XX + YY + ZZ = ZZ + (+- + -+)/2
    !   OPR2 = OPR1^2
    
    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(out) :: HSwap(dim,dim)
    real (c_double) :: OPR1(dimSpin1,dimSpin1,dimSpin1,dimSpin1)
    real (c_double) :: OPR2(dimSpin1,dimSpin1,dimSpin1,dimSpin1)
    real (c_double) :: OPR(dimSpin1,dimSpin1,dimSpin1,dimSpin1)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, idx1, idx2, idxp1, idxp2, idxs1, idxs2, s1, s2
    
    OPR1 = 0
    do idx1 = 1, dimSpin1
      do idx2 = 1, dimSpin1
        do idxp1 = 1, dimSpin1
          do idxp2 = 1, dimSpin1
            s1 = 2 - idx1
            s2 = 2 - idx2
            OPR1(idxp1, idxp2, idx1, idx2) = s1 * KDelta(idxp1, idx1) * s2 * KDelta(idxp2,idx2) + &
              & ( KDelta(idxp1, idx1+1) * KDelta(idxp2,idx2-1) + &
              &  KDelta(idxp1, idx1-1) * KDelta(idxp2, idx2+1) )
          enddo
        enddo
      enddo
    enddo

    OPR2 = 0
    do idx1 = 1, dimSpin1
      do idx2 = 1, dimSpin1
        do idxp1 = 1, dimSpin1
          do idxp2 = 1, dimSpin1

            do idxs1 = 1, dimSpin1
              do idxs2 = 1, dimSpin1
                OPR2(idxp1, idxp2, idx1, idx2) = OPR2(idxp1, idxp2, idx1, idx2) + &
                  & OPR1(idxp1, idxp2, idxs1, idxs2) * OPR1(idxs1, idxs2, idx1, idx2)
              enddo
            enddo
            OPR(idxp1, idxp2, idx1, idx2) = OPR2(idxp1, idxp2, idx1, idx2) + OPR1(idxp1, idxp2, idx1, idx2) !- &
            ! & KDelta(idxp1,idx1) * KDelta(idxp2,idx2)


          enddo
        enddo
      enddo
    enddo


!    print *, "OPR = "
!    call printmat(dimSpin1, OPR, 'R')

    HSwap = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin/2

        do idxp1 = 0, dimSpin1-1
          do idxp2 = 0, dimSpin1-1

            idx1 = config(2*k-1)
            idx2 = config(2*k)
            j = i + (idxp1 - idx1) * dimSpin1**(2*k-2) + (idxp2 - idx2) * dimSpin1**(2*k-1)
            HSwap(j+1,i+1) = HSwap(j+1,i+1) + OPR(idxp1+1, idxp2+1, idx1+1, idx2+1)

          enddo
        enddo

      enddo
    enddo


    !print*, "HSwap = "
    !call printmat(dim, HSwap, 'R')

  end subroutine buildHSwap


  subroutine buildHMBL(nspin, dim, Jint, Vint, hz, H)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, idx1, idx2, idxp1, idxp2, s1, s2
    real (c_double) :: OPRz(dimSpin1,dimSpin1), OPRzz(dimSpin1,dimSpin1,dimSpin1,dimSpin1)
    real (c_double) :: OPRxy(dimSpin1,dimSpin1,dimSpin1,dimSpin1)

    do idx1 = 1, dimSpin1
      do idx2 = 1, dimSpin1

        OPRz(idx2, idx1) = (2-idx1)*KDelta(idx2, idx1)
        do idxp1 = 1, dimSpin1
          do idxp2 = 1, dimSpin1
            s1 = 2 - idx1
            s2 = 2 - idx2
            OPRzz(idxp1, idxp2, idx1, idx2) = s1 * KDelta(idxp1, idx1) * s2 * KDelta(idxp2,idx2)
            OPRxy(idxp1, idxp2, idx1, idx2) = ( KDelta(idxp1, idx1+1) * KDelta(idxp2,idx2-1) + &
              &  KDelta(idxp1, idx1-1) * KDelta(idxp2, idx2+1) )
          enddo
        enddo

      enddo
    enddo

    H = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin-1

        do idxp1 = 0, dimSpin1-1
          do idxp2 = 0, dimSpin1-1

            idx1 = config(k)
            idx2 = config(k+1)
            j = i + (idxp1 - idx1) * dimSpin1**(k-1) + (idxp2 - idx2) * dimSpin1**k
            H(j+1,i+1) = H(j+1,i+1) + Vint(k) * OPRzz(idxp1+1, idxp2+1, idx1+1, idx2+1) + &
              & Jint(k) * OPRxy(idxp1+1, idxp2+1, idx1+1, idx2+1)

          enddo
          idx1 = config(k)
          j = i + (idxp1 - idx1) * dimSpin1**(k-1)
          H(j+1,i+1) = H(j+1,i+1) + hz(k) * OPRz(idxp1+1, idx1+1)

        enddo

      enddo

      k = nspin
      do idxp1 = 0, dimSpin1-1
        idx1 = config(k)
        j = i + (idxp1 - idx1) * dimSpin1**(k-1)
        H(j+1,i+1) = H(j+1,i+1) + hz(k) * OPRz(idxp1+1, idx1+1)
      enddo

    enddo

  end subroutine buildHMBL
!
!
!
!  subroutine buildSPARSE_HMBL(nspin, dim, nz_dim, Jint, Vint, hz, H, ROWS, COLS)
!
!    integer (c_int), intent(in) :: nspin, dim, nz_dim
!    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
!    real (c_double), intent(out) :: H(nz_dim)
!    integer (c_int), intent(out) :: ROWS(nz_dim), COLS(nz_dim)
!
!    integer :: config(nspin)
!    integer (c_int) :: i, j, k, m
!    
!    H = 0
!    ROWS = 0
!    COLS = 0
!    m = 0
!    do i = 1, dim
!
!    call decode(i-1,nspin,config)
!
!      m = m+1
!      do k = 1, nspin-1
!        H(m) = H(m) + Vint(k) * (1 - 2 * config(k)) * &
!          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
!      enddo
!      k = nspin
!      H(m) = H(m) + hz(k) * (1 - 2 * config(k))
!
!      ROWS(m) = i
!      COLS(m) = i
!
!      do k = 1, nspin-1
!
!        if (config(k)/=config(k+1)) then
!          m = m+1
!          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
!          H(m) = H(m) + Jint(k) * 2 * (config(k) - config(k+1))**2
!          COLS(m) = i
!          ROWS(m) = j
!        endif
!      enddo
!
!    enddo
!    !print *, m
!    !print *, "H_MBL_SPARSE ="
!    !do m = 1, (nspin+1)*dim/2
!    !  if(abs(H(m))<0e-10)  print *, H(m), ROWS(m), COLS(m), m
!    !enddo
!
!  end subroutine buildSPARSE_HMBL
!
!
!  subroutine buildHESPARSE_HMBL(nspin, dim, Jint, Vint, hz, H, ROWS, COLS)
!
!    integer (c_int), intent(in) :: nspin, dim
!    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
!    real (c_double), intent(out) :: H((nspin+3)*dim/4)
!    integer (c_int), intent(out) :: ROWS((nspin+3)*dim/4), COLS((nspin+3)*dim/4)
!
!    integer :: config(nspin)
!    integer (c_int) :: i, j, k, m
!    
!    H = 0
!    ROWS = 0
!    COLS = 0
!    m = 0
!    do i = 1, dim
!
!    call decode(i-1,nspin,config)
!
!      m = m+1
!      do k = 1, nspin-1
!        H(m) = H(m) + Vint(k) * (1 - 2 * config(k)) * &
!          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
!      enddo
!      k = nspin
!      H(m) = H(m) + hz(k) * (1 - 2 * config(k))
!
!      ROWS(m) = i
!      COLS(m) = i
!
!      do k = 1, nspin-1
!
!        if (config(k)/=config(k+1)) then
!          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
!          if(j>i) then
!            m = m+1
!            H(m) = H(m) + Jint(k) * 2 * (config(k) - config(k+1))**2
!            COLS(m) = i
!            ROWS(m) = j
!          endif
!        endif
!      enddo
!
!    enddo
!    !print *, m
!    !print *, "H_MBL_SPARSE ="
!    !do m = 1, (nspin+1)*dim/2
!      !if(abs(H(m))<0e-10)  print *, H(m), ROWS(m), COLS(m), m
!      !print *, H(m), ROWS(m), COLS(m), m
!    !enddo
!
!  end subroutine buildHESPARSE_HMBL
!
!
!  subroutine buildSz0_SPARSE_HMBL(nspin, dim_Sz0, nz_dim, Jint, Vint, hz, H, ROWS, COLS)
!
!    !
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0, nz_dim
!    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
!    real (c_double), intent(out) :: H(nz_dim)
!    integer (c_int), intent(out) :: ROWS(nz_dim), COLS(nz_dim)
!    
!
!    integer :: config(nspin), states(dim_Sz0)
!    integer (c_int) :: i, j, k, m, n, l, r
!
!    H = 0
!    ROWS = 0
!    COLS = 0
!    m = 0
!    call zero_mag_states(nspin, dim_Sz0, states)
!
!
!    do l = 1, dim_Sz0
!
!      i = states(l)
!      call decode(i,nspin,config)
!
!      m = m+1
!      do k = 1, nspin-1
!        H(m) = H(m) + Vint(k) * (1 - 2 * config(k)) * &
!          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
!      enddo
!      k = nspin
!      H(m) = H(m) + hz(k) * (1 - 2 * config(k))
!
!      ROWS(m) = l
!      COLS(m) = l
!
!      do k = 1, nspin-1
!
!        if (config(k)/=config(k+1)) then
!          m = m+1
!          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
!          H(m) = H(m) + Jint(k) * 2 * (config(k) - config(k+1))**2
!          n = binsearch(j,states)
!          COLS(m) = l
!          ROWS(m) = n
!        endif
!      enddo
!
!    enddo
!    n = m
!
!  end subroutine buildSz0_SPARSE_HMBL
!
!
!  subroutine buildSz0_HMBL(nspin, dim_Sz0, Jint, Vint, hz, H)
!
!    !dim_Sz0 is the dimensione of the Sz=0 subsapce
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
!    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)
!
!    integer :: config(nspin), states(dim_Sz0)
!    integer (c_int) :: i, j, k, m, n, l, r, rflag
!    
!    H = 0
!    call zero_mag_states(nspin, dim_Sz0, states)
!    do l = 1, dim_Sz0
!
!      i = states(l)
!      call decode(i,nspin,config)
!
!      do k = 1, nspin-1
!
!        H(l,l) = H(l,l) + Vint(k) * (1 - 2 * config(k)) * &
!          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
!
!        if (config(k)/=config(k+1)) then
!          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
!          r = binsearch(j,states)
!
!          H(r,l) = H(r,l) + Jint(k) * 2 * (config(k) - config(k+1))**2
!        endif
!      enddo
!      k = nspin
!      H(l,l) = H(l,l) + hz(k) * (1 - 2 * config(k))
!    enddo
!
!  end subroutine buildSz0_HMBL
!
!  subroutine buildSz0_HMBL_NNN(nspin, dim_Sz0, Jint, Vint, Vint2, hz, H)
!
!    !dim_Sz0 is the dimensione of the Sz=0 subsapce
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), Vint2(nspin-2), hz(nspin)
!    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)
!
!    integer :: config(nspin), states(dim_Sz0)
!    integer (c_int) :: i, j, k, m, n, l, r, rflag
!    
!    H = 0
!    call zero_mag_states(nspin, dim_Sz0, states)
!    do l = 1, dim_Sz0
!
!      i = states(l)
!      call decode(i,nspin,config)
!      config = 1 - 2*config
!
!      do k = 1, nspin-2
!
!        H(l,l) = H(l,l) + Vint(k) * config(k) * &
!          & config(k+1) + Vint2(k) * config(k) * config(k+2) + hz(k) * config(k)
!
!        if (config(k)/=config(k+1)) then
!          j = i + config(k)*2**(k-1) + config(k+1)*2**(k)
!          r = binsearch(j,states)
!
!          H(r,l) = H(r,l) + Jint(k) * (config(k) - config(k+1))**2/2
!        endif
!      enddo
!
!      k = nspin - 1
!      H(l,l) = H(l,l) + Vint(k) * config(k) * &
!        & config(k+1) + hz(k) * config(k)
!
!      if (config(k)/=config(k+1)) then
!        j = i + config(k)*2**(k-1) + config(k+1)*2**(k)
!        r = binsearch(j,states)
!
!        H(r,l) = H(r,l) + Jint(k) * (config(k) - config(k+1))**2/2
!      endif
!
!      k = nspin
!      H(l,l) = H(l,l) + hz(k) * config(k)
!    enddo
!
!  end subroutine buildSz0_HMBL_NNN
!
!
!
!  subroutine buildSz0_HSwap(nspin, dim_Sz0, H)
!
!    !dim_Sz0 is the dimensione of the Sz=0 subsapce
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)
!
!    integer :: config(nspin), states(dim_Sz0)
!    integer (c_int) :: i, j, k, m, n, l, r, rflag
! 
!    H = 0
!    call zero_mag_states(nspin, dim_Sz0, states)
!    do l = 1, dim_Sz0
!
!      i = states(l)
!      call decode(i,nspin,config)
!
!      do k = 1, nspin/2
!
!        H(l,l) = H(l,l) + (1 - 2 * config(2*k-1)) * &
!          & (1 - 2 * config(2*k))-1
!
!        if (config(2*k-1)/=config(2*k)) then
!          j = i + (1-2*config(2*k-1))*2**(2*k-2) + (1-2*config(2*k))*2**(2*k-1)
!          r = binsearch(j,states)
!          
!
!          H(r,l) = H(r,l) + 2 * (config(2*k-1) - config(2*k))**2
!        endif
!      enddo
!    enddo
!
!  end subroutine buildSz0_HSwap
!
!
!
!
  subroutine buildProdState(nspin, dim, alpha, state)

    !Builds a product state of the form: |phi> \otimes |phi> \otimes ...
    ! where: |phi> = alpha_1|+1> + alpha_2|0> + alpha_3|-1>
 
    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: alpha(dimSpin1)
    complex (c_double_complex), intent(out) :: state(dim)
    integer :: config(nspin), num(dimSpin1)
    complex (c_double_complex) :: alpha_n(dimSpin1)
    real (c_double) :: norm
    integer (c_int) :: i, k, idx

    alpha_n = alpha/dot_product(alpha,alpha)
    !print *, alpha_n

    !num() contains the number of all spins with a certain orientation
    ! (up, zero, down) <-> (0,1,2) in a given configuration (encoded in an integer i)
    ! It counts the number of 0, 1, 2 in 'config'
    state = 1
    do i = 1, dim
      call decode(i-1,nspin,config)
      !print *, "config = ", config(:)
      num = 0
      do k = 1, nspin
        num(config(k)+1) = num(config(k)+1) + 1
      enddo
      !print *, "num = ", num(:)

      do idx = 1, dimSpin1
        state(i) = state(i) * alpha_n(idx) ** num(idx)
        !print *, "alpha_n(", idx ,") =", alpha_n(idx) ** num(idx)
        !print *, "state(i) = ", state(i)
      enddo
      !print *, "state = ", state(i)

    enddo

  end subroutine buildProdState

!  subroutine buildNayakState(nspin,dim,alpha,beta,state)
!
!    integer (c_int), intent(in) :: nspin, dim
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim)
!    
!    integer :: config(nspin), i, k
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    do i = 1, dim
!
!      call decode(i-1,nspin,config)
!
!      do k = 1, nspin/2
!        state(i) = state(i) * ( (1-config(2*k-1)) * alpha_n + config(2*k-1) * beta_n ) * &
!          & ( config(2*k) * alpha_n + (1-config(2*k)) * beta_n )
!      enddo
!    enddo
!
!
!  end subroutine buildNayakState
!
!  subroutine buildNayakState_Sz0(nspin,dim_Sz0,alpha,beta,state)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim_Sz0)
!    
!    integer :: config(nspin), i, k, l, indx(dim_Sz0)
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    call zero_mag_states(nspin, dim_Sz0, indx)
!    do i = 1, dim_Sz0
!      l = indx(i)
!      call decode(l,nspin,config)
!
!      do k = 1, nspin/2
!        state(i) = state(i) * ( (1-config(2*k-1)) * alpha_n + config(2*k-1) * beta_n ) * &
!          & ( config(2*k) * alpha_n + (1-config(2*k)) * beta_n )
!      enddo
!    enddo
!
!
!  end subroutine buildNayakState_Sz0
!
!
!  subroutine buildI0LI1ProdState_Sz0(nspin, dim_Sz0, alpha, beta, state)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim_Sz0)
!    
!    integer (c_int) :: config(nspin), i, k, l, indx(dim_Sz0)
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    call zero_mag_states(nspin, dim_Sz0, indx)
!    do i = 1, dim_Sz0
!      l = indx(i)
!      call decode(l,nspin,config)
!
!      do k = 1, nspin/4
!        state(i) = state(i) * ( (1-config(4*k-3)) * config(4*k-2) * alpha_n + &
!          &  config(4*k-3) * (1-config(4*k-2)) * beta_n ) * &
!          & ( (1 - config(4*k-1)) * config(4*k) * beta_n + &
!          &  config(4*k-1) * (1-config(4*k)) * alpha_n )
!      enddo
!      if(mod(nspin/2,2) /= 0) then
!        k = nspin/4 + 1
!        state(i) = state(i) * ( (1-config(4*k-3)) * config(4*k-2) * alpha_n + &
!          &  config(4*k-3) * (1-config(4*k-2)) * beta_n )
!      endif
!    enddo
!
!
!  end subroutine
!
!  subroutine buildNeelState_Sz0(nspin, dim_Sz0, psi)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
!    
!    integer :: config(nspin), i, k, indx(dim_Sz0), l
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!    call zero_mag_states(nspin, dim_Sz0, indx)
!
!    psi = 0
!    l = 0
!    do k = 1, nspin/2
!      l = l + 2**(2*k-1)
!    enddo
!    !print *, l
!    
!    do i = 1, dim_Sz0
!      if(indx(i)==l) then
!        psi(i) = 1
!        exit
!      endif
!    enddo
!
!  end subroutine buildNeelState_Sz0
!
!  subroutine buildLI1ProdState_Sz0(nspin, dim_Sz0, alpha, beta, state)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim_Sz0)
!    
!    integer (c_int) :: config(nspin), i, k, l, indx(dim_Sz0)
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    call zero_mag_states(nspin, dim_Sz0, indx)
!    do i = 1, dim_Sz0
!      l = indx(i)
!      call decode(l,nspin,config)
!
!      do k = 1, nspin/2
!        state(i) = state(i) * ( (1-config(2*k-1)) * config(2*k) * alpha_n + &
!          &  config(2*k-1) * (1-config(2*k)) * beta_n )
!      enddo
!
!    enddo
!
!  end subroutine
!
!
!  subroutine buildRndLarge_IMB_LI_State_Sz0(nspin, dim_Sz0, IMB, LI, psi)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: IMB, LI
!    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
!
!    integer (c_int) :: idxSz0(dim_Sz0), nz, l, i, config(nspin)
!    integer (c_int), allocatable :: idx(:)
!    real (c_double) :: rand
!
!    call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, idxSz0)
!
!    do i = 1, dim_Sz0
!      if(idxSz0(i) == 0) then
!        nz = i-1
!        exit
!      endif
!    enddo
!
!    allocate(idx(nz))
!
!    idx = idxSz0(1:nz)
!
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    call init_random_seed()
!    do i = 1, nz
!      l = idx(i)
!      call random_number(rand)
!      psi(binsearch(l,idxSz0)) = rand
!    enddo
!    psi = psi/sqrt(dot_product(psi,psi))
!
!    !print *, imbalance_Sz0(nspin, dim_Sz0, psi), local_imbalance_Sz0(nspin, dim_Sz0, psi)
!
!  end subroutine
!
!
!  subroutine buildRndLarge_IMB_LI_State_basis_Sz0(nspin, dim_Sz0, IMB, LI, psi)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: IMB, LI
!    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
!
!    integer (c_int) :: idxSz0(dim_Sz0), nz, l, i, config(nspin)
!    integer (c_int), allocatable :: idx(:)
!    real (c_double) :: rand
!
!    call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, idxSz0)
!
!    do i = 1, dim_Sz0
!      if(idxSz0(i) == 0) then
!        nz = i-1
!        exit
!      endif
!    enddo
!
!    allocate(idx(nz))
!
!    idx = idxSz0(1:nz)
!
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    call init_random_seed()
!    call random_number(rand)
!    l = idx(floor(nz*rand + 1))
!    i = binsearch(l,idxSz0)
!
!    psi = 0
!    psi(i) = 1
!
!    call decode(l, nspin, config)
!    !print *, IMB, LI
!    !print *, imbalance_basis(nspin,l), local_imbalance_basis(nspin,l), l
!    1 format ( 4X,F6.3, 4X,F6.3, 4X,I4, 4X,I0 )
!
!
!  end subroutine
!    
!  subroutine buildLarge_IMB_LI_State_basis_Sz0(nspin, dim_Sz0, i, IMB, LI, psi)
!    integer (c_int), intent(in) :: nspin, dim_Sz0, i
!    real (c_double), intent(in) :: IMB, LI
!    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
!
!    integer (c_int) :: idxSz0(dim_Sz0), nz, j, l, config(nspin)
!    integer (c_int), allocatable :: idx(:)
!    real (c_double) :: rand
!
!    call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, idxSz0)
!
!    do j = 1, dim_Sz0
!      if(idxSz0(j) == 0) then
!        nz = j-1
!        exit
!      endif
!    enddo
!    print *, "nz = ", nz, "i = ", i
!    if (i>nz) stop "Error: Not enough number of states."
!
!    l = idxSz0(i)
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    j = binsearch(l, idxSz0)
!    psi = 0
!    psi(j) = 1
!
!    call decode(l, nspin, config)
!    print *, "build Large IMB, LI basis state"
!    print *, IMB, LI
!    print *, imbalance_basis(nspin,l), local_imbalance_basis(nspin,l), l
!    1 format ( 4X,F6.3, 4X,F6.3, 4X,I4, 4X,I0 )
!
!  end subroutine
!
!
!  subroutine buildState_Sz0_to_FullHS(nspin, dim, dim_Sz0, psi_Sz0, psi)
!    !Goes from the state in the subspace Sz=0 to the state in the full Hilbert 
!    !simply by constructing a vector which has zero components outside the Sz=0 subspace
!
!    integer (c_int), intent(in) :: nspin, dim, dim_Sz0
!    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
!    complex (c_double_complex), intent(out) :: psi(dim)
!
!    integer :: i, l, states(dim_Sz0)
!
!    call zero_mag_states(nspin, dim_Sz0, states)
!
!    psi = 0
!    !print *, "psi initialized"
!    do i = 1, dim_Sz0
!      l = states(i) + 1
!      psi(l) = psi_Sz0(i)
!      !print *,l, psi(l)
!      !print *, l, dim, i, dim_Sz0
!    enddo
!
!  end subroutine buildState_Sz0_to_FullHS
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
  subroutine zero_mag_states(nspin, dim_Sz0, states)

    !dim_Sz0 is the dimension of the Sz=0 subspace, the length of 'states'
    integer (c_int), intent(in) :: nspin, dim_Sz0
    integer (c_int), intent(out) :: states(dim_Sz0)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m, dim

    dim = dimSpin1**nspin
    k = 0
    states = 0
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(config)==nspin/2) then
        k = k+1
        states(k) = i
        !print *, i, k, states(k)
      endif
    enddo

  end subroutine zero_mag_states
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
!      if(abs(IMB) < 1.0e-10 .AND. sum_odd == sum_even) then
!        l = l+1
!        states(l) = i
!        print "(2(8X,I0),4X,*(I0))", i, l, config(:)
!      else if(abs(IMB) > 1.0e-10 .AND. abs(loc_imb - IMB) < 1.0e-10) then
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
!      if(abs(IMB) < 1.0e-10 .AND. sum_odd == sum_even) then
!        m = m+1
!        states(m) = l
!        print "(2(8X,I0),4X,*(I0))", l, m, config(:)
!      else if(abs(IMB) > 1.0e-10 .AND. abs(loc_imb - IMB) < 1.0e-10) then
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
!
!  subroutine magntz(i, nspin, mag)
!
!    integer(c_int), intent(in) :: i, nspin
!    integer(c_int) :: config(nspin)
!    real(c_double) :: mag
!    integer (c_int) :: j, k, m
!
!    call decode(i, nspin, config)
!
!    mag = 0
!    do k = 1, nspin
!      mag = mag + (1 - 2 * config(k))
!    enddo
!    
!  end subroutine magntz
!
!
!  real function mag_z(nspin, dim, state)
!
!    integer(c_int), intent(in) :: nspin, dim
!    complex(c_double_complex), intent(in) :: state(dim)
!    real(c_double) :: mag, magaux
!    integer :: config(nspin)
!    integer (c_int) :: i, j, k, m
!
!    mag = 0
!    do i = 1, dim
!      call decode(i-1,nspin,config)
!      magaux = 0
!      do k = 1, nspin
!        magaux = magaux + (1._c_double - 2._c_double * config(k))
!      enddo
!      magaux = magaux * abs(state(i))**2
!      mag = mag + magaux
!    enddo
!    mag = mag/nspin
!    mag_z = mag
!
!  end function mag_z
!
!  real function mag_z_p(nspin, dim, state, p)
!
!    integer(c_int), intent(in) :: nspin, dim, p
!    complex(c_double_complex), intent(in) :: state(dim)
!    real(c_double) :: mag
!    integer :: config(nspin)
!    integer (c_int) :: i, j, k, m
!
!    mag = 0
!    do i = 1, dim
!      call decode(i-1,nspin,config)
!      mag = mag + abs(state(i))**2 * (1 - 2*config(p))
!    enddo
!    mag_z_p = mag
!
!  end function mag_z_p
!
!  real function mag_stag_z(nspin, dim, state)
!
!    integer(c_int), intent(in) :: nspin, dim
!    complex(c_double_complex), intent(in) :: state(dim)
!    real(c_double) :: mag
!    real(c_double) :: magaux
!    integer :: config(nspin)
!    integer (c_int) :: i, j, k, m
!
!    mag = 0
!    do i = 1, dim
!      call decode(i-1,nspin,config)
!      magaux = 0
!      do k = 1, nspin
!        magaux = magaux + (-1)**k * (1._c_double - 2._c_double * config(k))
!      enddo
!      magaux = magaux * abs(state(i))**2
!      mag = mag + magaux
!    enddo
!    mag = mag/nspin
!    mag_stag_z = mag
!
!  end function mag_stag_z
!  
!  real function imbalance(nspin, dim, state)
!
!    integer(c_int), intent(in) :: nspin, dim
!    complex(c_double_complex), intent(in) :: state(dim)
!    real(c_double) :: imb, imbaux, mag, magaux
!    integer :: config(nspin)
!    integer (c_int) :: i, j, k, m
!
!    mag = 0
!    imb = 0
!    do i = 1, dim
!      call decode(i-1,nspin,config)
!      imbaux = 0
!      magaux = 0
!      do k = 1, nspin
!        imbaux = imbaux + (-1)**k * (1 - 2 * config(k))
!        magaux = magaux + (1 - 2*config(k))
!      enddo
!      magaux = magaux * abs(state(i))**2
!      imbaux = imbaux * abs(state(i))**2
!      mag = mag + magaux
!      imb = imb + imbaux
!    enddo
!    mag = mag/nspin
!    imb = imb/nspin
!    imb = imb/(1+mag)
!    imbalance = imb
!
!
!  end function imbalance
!
!  real function imbalance_Sz0(nspin, dim_Sz0, state)
!
!    integer(c_int), intent(in) :: nspin, dim_Sz0
!    complex(c_double_complex), intent(in) :: state(dim_Sz0)
!    real(c_double) :: imb, imbaux
!    integer :: config(nspin)
!    integer (c_int) :: i, j, k, m, l, indx(dim_Sz0)
!
!    imb = 0
!    call zero_mag_states(nspin, dim_Sz0, indx)
!
!    do i = 1, dim_Sz0
!
!      l = indx(i)
!      call decode(l,nspin,config)
!      imbaux = 0
!      do k = 1, nspin
!        imbaux = imbaux + (-1)**k * (1 - 2 * config(k))
!      enddo
!      imbaux = imbaux * abs(state(i))**2
!      imb = imb + imbaux
!    enddo
!    imb = imb/nspin
!    imbalance_Sz0 = imb
!
!  end function imbalance_Sz0
!
!  real function imbalance_sq_Sz0(nspin, dim_Sz0, state)
!
!    integer(c_int), intent(in) :: nspin, dim_Sz0
!    complex(c_double_complex), intent(in) :: state(dim_Sz0)
!    real(c_double) :: imb, imbaux
!    integer :: config(nspin)
!    integer (c_int) :: i, k1, k2, l, indx(dim_Sz0)
!
!    imb = 0
!    call zero_mag_states(nspin, dim_Sz0, indx)
!
!    do i = 1, dim_Sz0
!
!      l = indx(i)
!      call decode(l,nspin,config)
!      imbaux = 0
!      do k1 = 1, nspin
!        do k2 = 1, nspin
!          imbaux = imbaux + (-1)**k1 * (1 - 2 * config(k1)) * (-1)**k2 * (1 - 2 * config(k2))
!        enddo
!      enddo
!      imbaux = imbaux * abs(state(i))**2
!      imb = imb + imbaux
!    enddo
!    imb = imb/nspin**2
!    imbalance_sq_Sz0 = imb
!
!  end function imbalance_sq_Sz0
!
!  function imbalance_basis(nspin, i)
!
!    integer(c_int), intent(in) :: nspin, i
!    real (c_double) :: imb, imbalance_basis
!    integer :: config(nspin), k
!
!    call decode(i,nspin,config)
!    imb = 0
!    do k = 1, nspin
!        imb = imb + (-1)**k * (1 - 2 * config(k))
!        !print *, k, imb
!    enddo
!    imb = imb/(2*nspin - 2*sum(config))
!    imbalance_basis = imb
!
!  end function imbalance_basis
!
!  function imbalance_sq_basis(nspin, i)
!
!    integer(c_int), intent(in) :: nspin, i
!    real (c_double) :: imb, imbalance_sq_basis
!    integer :: config(nspin), k
!
!    call decode(i,nspin,config)
!    imb = 0
!    do k = 1, nspin
!        imb = imb + (-1)**k * (1 - 2 * config(k))
!        !print *, k, imb
!    enddo
!    imb = imb/(2*nspin - 2*sum(config))
!    imbalance_sq_basis = imb**2
!
!  end function imbalance_sq_basis
!
!  function local_imbalance_basis(nspin, i)
!    integer (c_int), intent(in) :: nspin, i
!    real (c_double) :: local_imbalance_basis
!
!    integer (c_int) :: k, config(nspin)
!    real (c_double) :: LI
!
!    LI = 0
!    call decode(i, nspin, config)
!    do k = 1, nspin/2
!
!      LI = LI + (config(2*k) - config(2*k-1))**2
!    enddo
!    LI = 2 * LI / nspin
!    local_imbalance_basis = LI
!
!  end function
!
!  function local_imbalance(nspin, dim, psi)
!    integer (c_int), intent(in) :: nspin, dim
!    complex (c_double_complex), intent(in) :: psi(dim)
!    real (c_double) :: local_imbalance
!
!    integer (c_int) :: i, k, config(nspin)
!    real (c_double) :: LI, LI_part
!
!    LI = 0
!    do i = 1, dim
!    call decode(i-1, nspin, config)
!      LI_part = 0
!      do k = 1, nspin/2
!        LI_part = LI_part + (config(2*k) - config(2*k-1))**2
!      enddo
!      LI = LI + abs(psi(i))**2 * 2 * LI_part / nspin
!    enddo
!    local_imbalance = LI
!
!  end function
!
!  function local_imbalance_Sz0(nspin, dim_Sz0, psi_Sz0)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
!    real (c_double) :: local_imbalance_Sz0
!
!    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)
!    real (c_double) :: LI, LI_part
!
!    call zero_mag_states(nspin, dim_Sz0, states)
!    LI = 0
!    do i = 1, dim_Sz0
!      l = states(i)
!      call decode(l, nspin, config)
!      LI_part = 0
!      do k = 1, nspin/2
!        LI_part = LI_part + (config(2*k) - config(2*k-1))**2
!      enddo
!      LI = LI + abs(psi_Sz0(i))**2 * 2 * LI_part / nspin
!    enddo
!    local_imbalance_Sz0 = LI
!
!  end function
!
!  function local_overlap(nspin, dim, psi)
!    integer (c_int), intent(in) :: nspin, dim
!    complex (c_double_complex), intent(in) :: psi(dim)
!    real (c_double) :: local_overlap
!
!    integer (c_int) :: i, j, k, config(nspin)
!    real (c_double) :: LO
!
!    LO = 0
!    do i = 1, dim
!      call decode(i-1, nspin, config)
!      do k = 1, nspin/2
!        if (config(2*k)/=config(2*k-1)) then
!          j = i + (config(2*k-1) - config(2*k)) * (2**(2*k-1) - 2**(2*k-2))
!          LO = LO + abs(psi(i))**2 - psi(i) * dconjg(psi(j))
!        endif
!      enddo
!    enddo
!    LO = LO/(nspin/2)
!    local_overlap = LO
!
!  end function
!
!  function local_overlap_Sz0(nspin, dim_Sz0, psi)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: psi(dim_Sz0)
!    real (c_double) :: local_overlap_Sz0
!
!    integer (c_int) :: i, j, k, r, l, config(nspin), states(dim_Sz0)
!    real (c_double) :: LO
!
!    call zero_mag_states(nspin, dim_Sz0, states)
!    LO = 0
!    do i = 1, dim_Sz0
!      l = states(i)
!      call decode(l, nspin, config)
!      do k = 1, nspin/2
!        if(config(2*k)/=config(2*k-1)) then
!          j = i + (config(2*k-1) - config(2*k)) * (2**(2*k-1) - 2**(2*k-2))
!          r = binsearch(j,states)
!          LO = LO + abs(psi(i))**2 - psi(i) * dconjg(psi(r))
!        endif
!      enddo
!    enddo
!    LO = LO/(nspin/2)
!    local_overlap_Sz0 = LO
!
!  end function
!
!  function sigmaz_corr_c(nspin, dim, q, p, psi)
!
!    integer (c_int), intent(in) :: nspin, dim, q, p
!    complex (c_double_complex), intent(in) :: psi(dim)
!    real (c_double) :: sigmaz_corr_c
!
!    integer (c_int) :: i, k, config(nspin)
!    real (c_double) :: corr, avgq, avgp
!
!    avgq = 0
!    avgp = 0
!    corr = 0
!    do i = 1, dim
!      
!      call decode(i-1,nspin,config)
!
!      corr = corr + abs(psi(i))**2 * (1 - 2 * config(q)) * (1 - 2 * config(p))
!      avgq = avgq + abs(psi(i))**2 * (1 - 2 * config(q))
!      avgp = avgp + abs(psi(i))**2 * (1 - 2 * config(p))
!        
!    enddo
!
!    sigmaz_corr_c = corr - avgq * avgp
!
!
!  end function
!
!  function sigmaz_corr_c_Sz0(nspin, dim_Sz0, q, p, psi_Sz0)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0, q, p
!    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
!    real (c_double) :: sigmaz_corr_c_Sz0
!
!    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)
!    real (c_double) :: corr, avgq, avgp
!
!    avgq = 0
!    avgp = 0
!    corr = 0
!    call zero_mag_states(nspin, dim_Sz0, states)
!    do i = 1, dim_Sz0
!      
!      l = states(i)
!      call decode(l,nspin,config)
!
!      corr = corr + abs(psi_Sz0(i))**2 * (1 - 2 * config(q)) * (1 - 2 * config(p))
!      avgq = avgq + abs(psi_Sz0(i))**2 * (1 - 2 * config(q))
!      avgp = avgp + abs(psi_Sz0(i))**2 * (1 - 2 * config(p))
!        
!    enddo
!
!    sigmaz_corr_c_Sz0 = corr - avgq * avgp
!
!
!  end function
!
!  subroutine local_zmag_Sz0(nspin, dim_Sz0, psi_Sz0, sigmaz)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
!    real (c_double), intent(out) :: sigmaz(nspin)
!
!    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)
!
!    call zero_mag_states(nspin, dim_Sz0, states)
!    sigmaz = 0
!    do i = 1, dim_Sz0
!      l = states(i)
!      call decode(l, nspin, config)
!      do k = 1, nspin
!        sigmaz(k) = sigmaz(k) + abs(psi_Sz0(i))**2 * (1-2*config(k))
!      enddo
!    enddo
!
!  end subroutine
!
!
!
!
!  function exact_energy(nspin, V_int, h_z, i)
!
!    integer (c_int), intent(in) :: nspin, i
!    real (c_double), intent(in) :: V_int(nspin), h_z(nspin)
!    real (c_double) :: exact_energy
!
!    integer (c_int) :: k, config(nspin)
!    real (c_double) :: E
!
!    call decode(i, nspin, config)
!    config = 1 - 2*config
!    E = 0
!    do k = 1, nspin - 1
!      E = E + V_int(k) * config(k) * config(k+1) + h_z(k) * config(k)
!    enddo
!    k = nspin
!    E = E + h_z(k) * config(k)
!
!    exact_energy = E
!
!  end function
!
!  subroutine exact_energies_Sz0(nspin, dim_Sz0, V_int, h_z, E)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: V_int(nspin), h_z(nspin)
!    real (c_double), intent(out) :: E(dim_Sz0)
!
!    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)
!
!    call zero_mag_states(nspin, dim_Sz0, states)
!    do i = 1, dim_Sz0
!
!      l = states(i)
!      call decode(l, nspin, config)
!      E(i) = exact_energy(nspin, V_int, h_z, l)
!
!    enddo
!
!  end subroutine
!
!  function exact_quasi_energy(nspin, V_int, h_z, i)
!
!    integer (c_int), intent(in) :: nspin, i
!    real (c_double), intent(in) :: V_int(nspin-1), h_z(nspin)
!    real (c_double) :: exact_quasi_energy
!
!    integer (c_int) :: k, m, config(nspin)
!    real (c_double) :: QE
!
!    call decode(i, nspin, config)
!    m = 0
!    do k = 1, nspin/2
!      m = m + 2**(2*k-2) * config(2*k) + 2**(2*k-1) * config(2*k-1)
!    enddo
!    QE = (exact_energy(nspin, V_int, h_z, i) + exact_energy(nspin, V_int, h_z, m)) / 2
!    exact_quasi_energy = QE
!
!  end function
!
!
!  subroutine exact_quasi_energies_Sz0(nspin, dim_Sz0, V_int, h_z, QE) !E, Es, QE, QE_alt)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: V_int(nspin-1), h_z(nspin)
!    real (c_double), intent(out) :: QE(dim_Sz0)
!
!    integer (c_int) :: i, k, l, m, config(nspin), states(dim_Sz0)
!    !real (c_double) :: E(dim_Sz0), Es(dim_Sz0), QE_alt(dim_Sz0)
!
!    QE = 0
!    !QE_alt = 0
!    !E = 0
!    !Es = 0
!    call zero_mag_states(nspin, dim_Sz0, states)
!    do i = 1, dim_Sz0
!
!      l = states(i)
!      call decode(l, nspin, config)
!      m = 0
!      do k = 1, nspin/2
!        m = m + 2**(2*k-2) * config(2*k) + 2**(2*k-1) * config(2*k-1)
!      enddo
!      QE(i) = (exact_energy(nspin, V_int, h_z, l) + exact_energy(nspin, V_int, h_z, m)) / 2
!
!      !print "(*(I0))", config(:)
!      !config = 1 - 2*config
!      !do k = 1, nspin/2 - 1
!      !  print *, "k = ", k
!      !  E(i) = E(i) + V_int(2*k-1) * config(2*k-1) * config(2*k) + V_int(2*k) * config(2*k) * config(2*k+1) + &
!      !   &  h_z(2*k-1) * config(2*k-1) + h_z(2*k) * config(2*k) 
!      !  Es(i) = Es(i) + V_int(2*k-1) * config(2*k) * config(2*k-1) + V_int(2*k) * config(2*k-1) * config(2*k+2) + &
!      !   &  h_z(2*k-1) * config(2*k) + h_z(2*k) * config(2*k-1) 
!
!      !  print *, "Normal contribution V"
!      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", V_int(2*k-1) * config(2*k-1) * config(2*k), &
!      !    & V_int(2*k) * config(2*k) * config(2*k+1), &
!      !    & 2*k-1, config(2*k-1), 2*k, config(2*k), 2*k+1, config(2*k+1)
!      !  print *, "Swap contribution V"
!      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", V_int(2*k-1) * config(2*k) * config(2*k-1), &
!      !    & V_int(2*k) * config(2*k-1) * config(2*k+2), &
!      !    & 2*k, config(2*k), 2*k-1, config(2*k-1), 2*k+2, config(2*k+2)
!      !  print *, "Normal contribution h"
!      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", h_z(2*k-1) * config(2*k-1), h_z(2*k) * config(2*k), &
!      !    & 2*k-1, config(2*k-1), 2*k, config(2*k)
!      !  print *, "Swap contribution h"
!      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", h_z(2*k-1) * config(2*k), h_z(2*k) * config(2*k-1), &
!      !    & 2*k, config(2*k), 2*k-1, config(2*k-1)
!      !  print *, ""
!
!      !enddo
!      !k = nspin/2
!      !  print *, "k = ", k
!      !E(i) = E(i) + V_int(2*k-1) * config(2*k-1) * config(2*k) + &
!      !  & + h_z(2*k-1) * config(2*k-1) + h_z(2*k) * config(2*k)
!      !Es(i) = Es(i) + V_int(2*k-1) * config(2*k) * config(2*k-1) + &
!      !  & + h_z(2*k-1) * config(2*k) + h_z(2*k) * config(2*k-1)
!      !  print *, "Normal contribution V"
!      !  print "( F20.15,4X, *(2I0,2X)  )", V_int(2*k-1) * config(2*k-1) * config(2*k), &
!      !    & 2*k-1, config(2*k-1), 2*k, config(2*k)
!      !  print *, "Swap contribution V"
!      !  print "( F20.15,4X, *(2I0,2X)  )", V_int(2*k-1) * config(2*k) * config(2*k-1), &
!      !    & 2*k, config(2*k), 2*k-1, config(2*k-1)
!      !  print *, "Normal contribution h"
!      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", h_z(2*k-1) * config(2*k-1), h_z(2*k) * config(2*k), &
!      !    & 2*k-1, config(2*k-1), 2*k, config(2*k)
!      !  print *, "Swap contribution h"
!      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", h_z(2*k-1) * config(2*k), h_z(2*k) * config(2*k-1), &
!      !    & 2*k, config(2*k), 2*k-1, config(2*k-1)
!      !  print *, ""
!
!      !print "( F15.10,4X,*(I0)  )", E(i), (1 - config(:))/2
!      !print *, E(i), Es(i), exact_energy(nspin, V_int, h_z, l), exact_energy(nspin, V_int, h_z, m)
!
!      !QE_alt(i) = (E(i) + Es(i)) / 2
!
!      !do k = 1, nspin/2 - 1
!      !  QE(i) = QE(i) + ( ( V_int(2*k-1) * 2 * config(2*k-1) * config(2*k) + &
!      !    & V_int(2*k) * ( config(2*k-1) * config(2*k+2) + config(2*k) * config(2*k+1) )) + &
!      !    & ( h_z(2*k-1) + h_z(2*k) ) * (config(2*k) + config(2*k-1)) )/2
!      !enddo
!      !k = nspin/2
!      !QE(i) = QE(i) + ( ( V_int(2*k-1) * 2 * config(2*k-1) * config(2*k) ) + &
!      !  & ( h_z(2*k-1) + h_z(2*k) ) * (config(2*k) + config(2*k-1)) )/2
!
!      !print "( A,4X, A,4X, A  )", "(E + Esigma)/2", "QE", "config"
!      !print "( F15.10,4X, F15.10,4X, *(I0)  )", QE_alt(i), QE(i), (1 - config(:))/2 
!    enddo
!
!  end subroutine
!
!
!  subroutine time_avg(option, steps, start, avg, sigma, t_avg, t_sigma)
!    character, intent(in) :: option*1
!    integer (c_int), intent(in) :: steps, start
!    real (c_double), intent(in) :: avg(steps), sigma(steps)
!    real (c_double), intent(out) :: t_avg, t_sigma
!    integer (c_int) :: i, j, k
!    
!    t_avg = 0
!    t_sigma = 0
!    if (option == 'F') then
!      do j = start, steps
!        t_avg = t_avg + avg(j)
!        t_sigma = t_sigma + sigma(j)**2
!      enddo
!    else if (option == 'T') then
!      do j = start, steps
!        t_avg = t_avg + 2*(mod(j,2)-0.5) * avg(j)
!        t_sigma = t_sigma + sigma(j)**2
!      enddo
!    endif
!    t_avg = t_avg/real(steps-start+1,c_double)
!    t_sigma = sqrt(t_sigma/real(steps-start+1,c_double))
!
!  end subroutine time_avg
!
!
!  integer (c_int) function binsearch(val, array)
!  
!  
!    implicit none
!    integer (c_int), intent(in) :: val, array(:)
!    integer (c_int) :: mid, start, finish, range
!    
!    binsearch = -1
!    start = 1
!    finish = size(array)
!    
!    range = finish - start
!    mid = (finish + start) / 2
!    
!    do while (array(mid) /= val .and. range > 0)
!    
!      if (val > array(mid)) then
!        start = mid + 1
!      else
!        finish = mid - 1
!      end if
!      
!      range = finish - start
!      mid = (start + finish) / 2
!    
!    end do
!    
!    if (array(mid) == val) then
!      binsearch = mid
!    end if
!  
!  end function binsearch




end module genmat
