module genmat

  use exponentiate
  use iso_c_binding
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)



contains

  subroutine decode(decimal, nbits, bitstring)

    integer (c_int), intent(in) :: decimal, nbits
    integer (c_int), intent(out) :: bitstring(nbits)
    integer(c_int) :: i

    bitstring = 0
    do i = 1, nbits
      if (btest(decimal, i - 1)) bitstring(i) = 1
    enddo
  end subroutine decode

  subroutine buildSx(nspin, dim, Sx)

    integer (c_int) :: nspin, dim
    real (c_double):: Sx(dim,dim)

    integer (c_int) :: config(nspin)
    integer (c_int) :: i, j, k

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
    integer(c_int) :: i, j, k

    !print*,"Start buildHx"

    Hx = 0
    !print*,"Start buildHx2"

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
    integer(c_int) :: i, j, k

    if (dim > 2**nspin) stop "Error: Value of dimension exceeds 2^nspin"
    if (p > nspin) stop "Error: Value of Spin Number 'p' exceeds totale spin number"


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
    integer(c_int) :: i, j, k

    if (dim > 2**nspin) stop "Error: Value of dimension exceeds 2^nspin"
    if (p > nspin) stop "Error: Value of Spin Number 'p' exceeds totale spin number"

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
    integer(c_int) :: i, j, k

    Sy = C_ZERO

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

    integer :: config(nspin)
    integer(c_int) :: i, j, k

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
    integer(c_int) :: i, j, k

    if (nspin < 2) stop "Error: Can't construct interaction without at least 2 spspins"
    !print*,"Start buildHJ"

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
    integer(c_int) :: i, j, k

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

  subroutine buildUFNayak(nspin, dim, Jint, h_x, h_z, T0, T1, UF )

    integer (c_int), intent(in) :: nspin, dim
    real(c_double), intent(in)  :: Jint(nspin-1), h_x(nspin), h_z(nspin), T0, T1
    complex(c_double_complex), intent(out)  :: UF(dim,dim)

    real (c_double) :: H(dim,dim), E(dim), W(dim,dim)
    complex (c_double_complex) :: Ux(dim,dim), U_nayak(dim,dim)

    !Sx si può leggere da file invece che essere generato
    !filestring = 'matrices/spin/Sx_nspin...'
    !open(newunit = u_int, file=filestring)
    !read (u_int, *) Sx
    !close(u_int)

    call buildSx(nspin, dim, H)
    call diagSYM( 'V', dim, H, E, W )
    call expSYM( dim, C_UNIT*T1, E, W, Ux )

    call buildHNayak( nspin, dim, Jint, h_x, h_z, H  )
    call diagSYM( 'V', dim, H, E, W )
    call expSYM( dim, -C_UNIT*T0, E, W, U_nayak )

    UF = matmul(U_nayak, Ux)

  end subroutine buildUFNayak



  subroutine buildHSwap(nspin, dim, HSwap)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(out) :: HSwap(dim,dim)

    integer :: config(nspin)
    integer(c_int) :: i, j, k

    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
    !print*,"Start buildHJ"
    
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

  subroutine buildUFSwap(nspin, dim, Jint, h_x, h_z, T0, T1, UF )

    integer (c_int), intent(in) :: nspin, dim
    real(c_double), intent(in)  :: Jint(nspin-1), h_x(nspin), h_z(nspin), T0, T1
    complex(c_double_complex), intent(out)  :: UF(dim,dim)

    real (c_double) :: H(dim,dim), E(dim), W(dim,dim)
    complex (c_double_complex) :: USwap(dim,dim), U_nayak(dim,dim)

    !Sx si può leggere da file invece che essere generato
    !filestring = 'matrices/spin/Sx_nspin...'
    !open(newunit = u_int, file=filestring)
    !read (u_int, *) Sx
    !close(u_int)

    call buildHSwap(nspin, dim, H)
    call diagSYM( 'V', dim, H, E, W )
    call expSYM( dim, -C_UNIT*T1, E, W, USwap )

    call buildHNayak( nspin, dim, Jint, h_x, h_z, H )
    call diagSYM( 'V', dim, H, E, W )
    call expSYM( dim, -C_UNIT*T0, E, W, U_nayak )

    UF = matmul(U_nayak, USwap)

  end subroutine buildUFSwap


  subroutine buildHMBL(nspin, dim, Jint, Vint, hx, hz, H)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hx(nspin), hz(nspin)
    real (c_double), intent(out) :: H(dim,dim)

    integer :: config(nspin)
    integer(c_int) :: i, j, k, m

    
    H = 0

    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin-1

        m = i + (1-2*config(k)) * 2**(k-1)
        j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)

        H(i+1,i+1) = H(i+1,i+1) + Jint(k) * (1 - 2 * config(k)) * &
          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
        H(j+1,i+1) = H(j+1,i+1) + Vint(k) * 2 * (config(k) - config(k+1))**2
        H(m+1,i+1) = H(m+1,i+1) + hx(k)
      enddo
      k = nspin
      m = i + (1-2*config(k)) * 2**(k-1)
      H(i+1,i+1) = H(i+1,i+1) + hz(k) * (1 - 2 * config(k))
      H(m+1,i+1) = H(m+1,i+1) + hx(k)
    enddo

  end subroutine buildHMBL



  subroutine magntz(i, nspin, mag)

    integer(c_int), intent(in) :: i, nspin
    integer(c_int) :: config(nspin)
    real(c_double) :: mag
    integer :: k

    call decode(i, nspin, config)

    mag = 0
    do k = 1, nspin
      mag = mag + (1 - 2 * config(k))
    enddo
    
  end subroutine magntz


  real function mag_z(nspin, dim, state)

    integer(c_int), intent(in) :: nspin, dim
    complex(c_double_complex), intent(in) :: state(dim)
    real(c_double) :: mag
    real(c_double) :: magaux
    integer :: i, j, k, config(nspin)

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
    integer :: i, config(nspin)

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
    integer :: i, j, k, config(nspin)

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
    integer :: i, j, k, config(nspin)

    mag = 0
    imb = 0
    do i = 0, dim-1
      call decode(i,nspin,config)
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


  end function mag_stag_z

  subroutine buildProdState(nspin, dim, alpha, beta, state)

    !Builds a product state of the form: (alpha|up> + beta|down) \otimes (alpha|up> + beta|down) \otimes ...
    
    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: alpha, beta
    complex (c_double_complex), intent(out) :: state(dim)
    integer :: i, config(nspin), un_vec(nspin), n_down
    complex (c_double_complex) :: alpha_n, beta_n
    real (c_double) :: norm

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
    
    integer :: i, k, config(nspin)
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




end module genmat
