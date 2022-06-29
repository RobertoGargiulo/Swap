module genmat

  use iso_c_binding
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)



contains

  subroutine decode(decimal, nbits, bitstring)

    integer (c_int), intent(in) :: decimal, nbits
    integer (c_int), intent(out) :: bitstring(nbits)
    integer (c_int) :: i, j, k, m

    bitstring = 0
    do i = 1, nbits
      if (btest(decimal, i - 1)) bitstring(i) = 1
    enddo
  end subroutine decode

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

  subroutine buildSx(nspin, dim, Sx)

    integer (c_int) :: nspin, dim
    real (c_double):: Sx(dim,dim)

    integer (c_int) :: config(nspin)
    integer (c_int) :: i, j, k, m

    Sx = 0
    !print*, "Start Sx"
    
    do i = 0, dim - 1
      
      call decode(i, nspin, config)

      do k = 1, nspin
        j = i + (1-2*config(k)) * 2**(k-1)
        Sx(j+1,i+1) = Sx(j+1,i+1) + 1
      enddo

    enddo

  end subroutine buildSx

  subroutine buildHx(nspin,dim,h_x,Hx)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: h_x(nspin)
    real (c_double), intent(out) :: Hx(dim,dim)

    integer (c_int) :: config(nspin)
    integer (c_int) :: i, j, k, m


    Hx = 0
    do i = 0, dim - 1

      call decode(i, nspin, config)

      do k = 1, nspin

        j = i + (1-2*config(k)) * 2**(k-1)
        Hx(j+1,i+1) = Hx(j+1,i+1) +  h_x(k)

      enddo
    enddo

  end subroutine buildHx

  subroutine buildSxp(nspin,dim,p,Sxp)

    integer (c_int), intent(in) :: nspin, dim, p
    real (c_double), intent(out) :: Sxp(dim,dim)

    integer (c_int) :: config(nspin)
    integer (c_int) :: i, j, k, m

    Sxp = 0
    do i = 0, dim - 1

      call decode(i, nspin, config)

      j = i + (1-2*config(p)) * 2**(p-1)
      Sxp(j+1,i+1) = 1

    enddo

  end subroutine buildSxp

  subroutine buildSyp(nspin,dim,p,Syp)

    integer (c_int), intent(in) :: nspin, dim, p
    complex (c_double_complex), intent(out) :: Syp(dim,dim)

    integer (c_int) :: config(nspin)
    integer (c_int) :: i, j, k, m

    Syp = 0
    do i = 0, dim - 1

      call decode(i, nspin, config)

      j = i + (1-2*config(p)) * 2**(p-1)
      Syp(j+1,i+1) = C_UNIT * (1._c_double - 2._c_double * config(p))

    enddo

  end subroutine buildSyp

  subroutine buildSy(nspin, dim, Sy)

    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(out) :: Sy(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    Sy = 0
    do i = 0, dim - 1

      call decode(i, nspin, config)
 
      do k = 1, nspin
 
         j = i + (1-2*config(k)) * 2**(k-1)

         Sy(j+1,i+1) = Sy(j+1,i+1) + C_UNIT * (1._c_double - 2._c_double * config(k))

      enddo
    enddo

  end subroutine buildSy


  subroutine buildSz(nspin, dim, Sz)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(out) :: Sz(dim,dim)
    integer (c_int) :: i, j, k, m

    integer :: config(nspin)

    Sz = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin
        Sz(i+1,i+1) = Sz(i+1,i+1) + (1 - 2 * config(k))
      enddo

    enddo
  end subroutine buildSz

  subroutine buildHJ(nspin, dim, Jint, HJ)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1)
    real (c_double), intent(out) :: HJ(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    HJ = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin - 1
        HJ(i+1,i+1) = HJ(i+1,i+1) + Jint(k) * (1 - 2 * config(k)) * &
         & (1 - 2 * config(k+1))
      enddo
    enddo

    !print*,"End buildHJ"
  end subroutine buildHJ


  subroutine buildHNayak( nspin, dim, Jint, hx, hz, H )

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1), hx(nspin), hz(nspin)
    real (c_double), intent(out) :: H(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    H = 0
    do i = 0, dim - 1

      call decode(i, nspin, config)

      do k = 1, nspin - 1 !Open Boundary Conditions

        H(i+1,i+1) = H(i+1,i+1) + Jint(k) * (1 - 2 * config(k)) * &
        & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))

        j = i + (1-2*config(k)) * 2**(k-1)
        H(j+1,i+1) = H(j+1,i+1) + hx(k)

      enddo
      H(i+1,i+1) = H(i+1,i+1) + hz(nspin) * (1 - 2 * config(nspin))
      j = i + (1-2*config(nspin)) * 2**(nspin-1)
      H(j+1,i+1) = H(j+1,i+1) + hx(nspin)
    enddo
  end subroutine buildHNayak


  subroutine buildHSwap(nspin, dim, HSwap)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(out) :: HSwap(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
    
    HSwap = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin/2

        j = i + (1-2*config(2*k-1))*2**(2*k-2) + (1-2*config(2*k))*2**(2*k-1)

        HSwap(i+1,i+1) = HSwap(i+1,i+1) + (1 - 2 * config(2*k-1)) * &
          & (1 - 2 * config(2*k))-1

        HSwap(j+1,i+1) = HSwap(j+1,i+1) + 2 * (config(2*k) - config(2*k-1))**2

      enddo
    enddo

  end subroutine buildHSwap


  subroutine buildHMBL(nspin, dim, Jint, Vint, hz, H)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m
    
    H = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin-1

        j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)

        H(i+1,i+1) = H(i+1,i+1) + Vint(k) * (1 - 2 * config(k)) * &
          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
        H(j+1,i+1) = H(j+1,i+1) + Jint(k) * 2 * (config(k) - config(k+1))**2
      enddo
      k = nspin
      H(i+1,i+1) = H(i+1,i+1) + hz(k) * (1 - 2 * config(k))
    enddo

  end subroutine buildHMBL



  subroutine buildSPARSE_HMBL(nspin, dim, Jint, Vint, hz, H, ROWS, COLS)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H((nspin+1)*dim/2)
    integer (c_int), intent(out) :: ROWS((nspin+1)*dim/2), COLS((nspin+1)*dim/2)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m
    
    H = 0
    ROWS = 0
    COLS = 0
    m = 0
    do i = 1, dim

    call decode(i-1,nspin,config)

      m = m+1
      do k = 1, nspin-1
        H(m) = H(m) + Vint(k) * (1 - 2 * config(k)) * &
          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
      enddo
      k = nspin
      H(m) = H(m) + hz(k) * (1 - 2 * config(k))

      ROWS(m) = i
      COLS(m) = i

      do k = 1, nspin-1

        if (config(k)/=config(k+1)) then
          m = m+1
          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
          H(m) = H(m) + Jint(k) * 2 * (config(k) - config(k+1))**2
          COLS(m) = i
          ROWS(m) = j
        endif
      enddo

    enddo
    !print *, m
    !print *, "H_MBL_SPARSE ="
    !do m = 1, (nspin+1)*dim/2
    !  if(abs(H(m))<0e-10)  print *, H(m), ROWS(m), COLS(m), m
    !enddo

  end subroutine buildSPARSE_HMBL


  subroutine buildHESPARSE_HMBL(nspin, dim, Jint, Vint, hz, H, ROWS, COLS)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H((nspin+3)*dim/4)
    integer (c_int), intent(out) :: ROWS((nspin+3)*dim/4), COLS((nspin+3)*dim/4)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m
    
    H = 0
    ROWS = 0
    COLS = 0
    m = 0
    do i = 1, dim

    call decode(i-1,nspin,config)

      m = m+1
      do k = 1, nspin-1
        H(m) = H(m) + Vint(k) * (1 - 2 * config(k)) * &
          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
      enddo
      k = nspin
      H(m) = H(m) + hz(k) * (1 - 2 * config(k))

      ROWS(m) = i
      COLS(m) = i

      do k = 1, nspin-1

        if (config(k)/=config(k+1)) then
          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
          if(j>i) then
            m = m+1
            H(m) = H(m) + Jint(k) * 2 * (config(k) - config(k+1))**2
            COLS(m) = i
            ROWS(m) = j
          endif
        endif
      enddo

    enddo
    print *, m
    print *, "H_MBL_SPARSE ="
    do m = 1, (nspin+1)*dim/2
      !if(abs(H(m))<0e-10)  print *, H(m), ROWS(m), COLS(m), m
      print *, H(m), ROWS(m), COLS(m), m
    enddo

  end subroutine buildHESPARSE_HMBL










  subroutine zero_mag_states(nspin, dim, states)
    integer (c_int), intent(in) :: nspin, dim
    integer (c_int), intent(out) :: states(dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    k = 0
    do i = 0, dim-1

      call decode(i, nspin, config)

      if (sum(config)==nspin/2) then
        k = k+1
        states(k) = i
        !print *, states(p), p
        !print '(1X,I0)', config(:)
      endif
    enddo

  end subroutine zero_mag_states



  subroutine magntz(i, nspin, mag)

    integer(c_int), intent(in) :: i, nspin
    integer(c_int) :: config(nspin)
    real(c_double) :: mag
    integer (c_int) :: j, k, m

    call decode(i, nspin, config)

    mag = 0
    do k = 1, nspin
      mag = mag + (1 - 2 * config(k))
    enddo
    
  end subroutine magntz


  real function mag_z(nspin, dim, state)

    integer(c_int), intent(in) :: nspin, dim
    complex(c_double_complex), intent(in) :: state(dim)
    real(c_double) :: mag, magaux
    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    mag = 0
    do i = 1, dim
      call decode(i-1,nspin,config)
      magaux = 0
      do k = 1, nspin
        magaux = magaux + (1._c_double - 2._c_double * config(k))
        !print *, i, k, config(:)
      enddo
      magaux = magaux * abs(state(i))**2
      mag = mag + magaux
      !print *, i, mag, config(:)
    enddo
    mag = mag/nspin
    mag_z = mag

  end function mag_z

  real function mag_z_p(nspin, dim, state, p)

    integer(c_int), intent(in) :: nspin, dim, p
    complex(c_double_complex), intent(in) :: state(dim)
    real(c_double) :: mag
    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    mag = 0
    do i = 1, dim
      call decode(i-1,nspin,config)
      mag = mag + abs(state(i))**2 * (1 - 2*config(p))
    enddo
    mag_z_p = mag

  end function mag_z_p

  real function mag_stag_z(nspin, dim, state)

    integer(c_int), intent(in) :: nspin, dim
    complex(c_double_complex), intent(in) :: state(dim)
    real(c_double) :: mag
    real(c_double) :: magaux
    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    mag = 0
    do i = 1, dim
      call decode(i-1,nspin,config)
      magaux = 0
      do k = 1, nspin
        magaux = magaux + (-1)**k * (1._c_double - 2._c_double * config(k))
        !print *, i, k, config(:)
      enddo
      magaux = magaux * abs(state(i))**2
      mag = mag + magaux
      !print *, i, mag, config(:)
    enddo
    mag = mag/nspin
    mag_stag_z = mag

  end function mag_stag_z
  
  real function imbalance(nspin, dim, state)

    integer(c_int), intent(in) :: nspin, dim
    complex(c_double_complex), intent(in) :: state(dim)
    real(c_double) :: imb, imbaux, mag, magaux
    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    mag = 0
    imb = 0
    do i = 1, dim
      call decode(i-1,nspin,config)
      imbaux = 0
      magaux = 0
      do k = 1, nspin
        imbaux = imbaux + (-1)**k * (1 - 2 * config(k))
        magaux = magaux + (1 - 2*config(k))
        !print *, i, k, config(:)
      enddo
      magaux = magaux * abs(state(i))**2
      imbaux = imbaux * abs(state(i))**2
      mag = mag + magaux
      imb = imb + imbaux
      !print *, i, mag, config(:)
    enddo
    mag = mag/nspin
    imb = imb/nspin
    imb = imb/(1+mag)
    imbalance = imb


  end function imbalance

  subroutine buildProdState(nspin, dim, alpha, beta, state)

    !Builds a product state of the form: (alpha|up> + beta|down) \otimes (alpha|up> + beta|down) \otimes ...
    
    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: alpha, beta
    complex (c_double_complex), intent(out) :: state(dim)
    integer :: config(nspin), un_vec(nspin), n_down
    complex (c_double_complex) :: alpha_n, beta_n
    real (c_double) :: norm
    integer (c_int) :: i, j, k, m

    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
    alpha_n = alpha/norm
    beta_n = beta/norm

    state = 0
    un_vec = 1
    do i = 1, dim
      call decode(i-1,nspin,config)
      n_down = dot_product(un_vec,config)
      state(i) = alpha**(nspin-n_down) * beta**n_down

    enddo

  end subroutine buildProdState

  subroutine buildStaggState(nspin,dim,alpha,beta,state)

    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: alpha, beta
    complex (c_double_complex), intent(out) :: state(dim)
    
    integer :: config(nspin), i, k
    complex (c_double_complex) :: alpha_n, beta_n
    real (c_double) :: norm

    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"

    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
    alpha_n = alpha/norm
    beta_n = beta/norm
    state = 1
    do i = 1, dim

      call decode(i-1,nspin,config)

      do k = 1, nspin/2
        state(i) = state(i) * ( (1-config(2*k-1)) * alpha_n + config(2*k-1) * beta_n ) * &
          & ( config(2*k) * alpha_n + (1-config(2*k)) * beta_n )
      enddo
    enddo


  end subroutine buildStaggState

  subroutine time_avg(steps, start, avg, sigma, t_avg, t_sigma)
    integer (c_int), intent(in) :: steps, start
    real (c_double), intent(in) :: avg(steps), sigma(steps)
    real (c_double), intent(out) :: t_avg, t_sigma
    integer (c_int) :: i, j, k
    
    t_avg = 0
    t_sigma = 0
    do j = start, steps
      t_avg = t_avg + avg(j)
      t_sigma = t_sigma + sigma(j)**2
    enddo
    t_avg = t_avg/(steps-start)
    t_sigma = sqrt(t_sigma/real(steps-start))

  end subroutine time_avg








end module genmat
