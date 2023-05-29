module matrices
  
  use functions
  !use ifport
  use iso_c_binding
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  integer (c_int), private, parameter :: dimSpinHalf = 2
  real (c_double), private, parameter :: tol = 1.0e-10


contains


  subroutine buildSx(nspin, dim, Sx)

    integer (c_int) :: nspin, dim
    real (c_double):: Sx(dim,dim)

    integer (c_int) :: config(nspin)
    integer (c_int) :: i, j, k, m

    Sx = 0
    
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

  subroutine buildHJ(nspin, dim, Jzz, HJ)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jzz(nspin-1)
    real (c_double), intent(out) :: HJ(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    HJ = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin - 1
        HJ(i+1,i+1) = HJ(i+1,i+1) + Jzz(k) * (1 - 2 * config(k)) * &
         & (1 - 2 * config(k+1))
      enddo
    enddo

  end subroutine buildHJ


  subroutine buildHNayak( nspin, dim, hx, Jzz, hz, H )

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jzz(nspin-1), hx(nspin), hz(nspin)
    real (c_double), intent(out) :: H(dim,dim)

    integer :: config(nspin), spin(nspin)
    integer (c_int) :: i, j, k, m

    H = 0
    do i = 0, dim - 1

      call decode(i, nspin, config)
      spin = 1 - 2 * config

      do k = 1, nspin - 1 !Open Boundary Conditions

        H(i+1,i+1) = H(i+1,i+1) + Jzz(k) * spin(k) * spin(k+1) + &
        & hz(k) * spin(k)

        j = i + (1-2*config(k)) * 2**(k-1)
        H(j+1,i+1) = H(j+1,i+1) + hx(k)

      enddo
      H(i+1,i+1) = H(i+1,i+1) + hz(nspin) * spin(nspin)
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


  subroutine buildHMBL(nspin, dim, Jxy, Vzz, hz, H)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m
    
    H = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin-1

        j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)

        H(i+1,i+1) = H(i+1,i+1) + Vzz(k) * (1 - 2 * config(k)) * &
          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
        H(j+1,i+1) = H(j+1,i+1) + Jxy(k) * 2 * (config(k) - config(k+1))**2
      enddo
      k = nspin
      H(i+1,i+1) = H(i+1,i+1) + hz(k) * (1 - 2 * config(k))
    enddo

  end subroutine buildHMBL



  subroutine buildSPARSE_HMBL(nspin, dim, nz_dim, Jxy, Vzz, hz, H, ROWS, COLS)

    integer (c_int), intent(in) :: nspin, dim, nz_dim
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H(nz_dim)
    integer (c_int), intent(out) :: ROWS(nz_dim), COLS(nz_dim)

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
        H(m) = H(m) + Vzz(k) * (1 - 2 * config(k)) * &
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
          H(m) = H(m) + Jxy(k) * 2 * (config(k) - config(k+1))**2
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


  subroutine buildHESPARSE_HMBL(nspin, dim, Jxy, Vzz, hz, H, ROWS, COLS)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1), hz(nspin)
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
        H(m) = H(m) + Vzz(k) * (1 - 2 * config(k)) * &
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
            H(m) = H(m) + Jxy(k) * 2 * (config(k) - config(k+1))**2
            COLS(m) = i
            ROWS(m) = j
          endif
        endif
      enddo

    enddo
    !print *, m
    !print *, "H_MBL_SPARSE ="
    !do m = 1, (nspin+1)*dim/2
      !if(abs(H(m))<0e-10)  print *, H(m), ROWS(m), COLS(m), m
      !print *, H(m), ROWS(m), COLS(m), m
    !enddo

  end subroutine buildHESPARSE_HMBL


  subroutine buildSz0_SPARSE_HMBL(nspin, dim_Sz0, nz_dim, Jxy, Vzz, hz, H, ROWS, COLS)

    !

    integer (c_int), intent(in) :: nspin, dim_Sz0, nz_dim
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H(nz_dim)
    integer (c_int), intent(out) :: ROWS(nz_dim), COLS(nz_dim)
    

    integer :: config(nspin), states(dim_Sz0)
    integer (c_int) :: i, j, k, m, n, l, r

    H = 0
    ROWS = 0
    COLS = 0
    m = 0
    call zero_mag_states(nspin, dim_Sz0, states)


    do l = 1, dim_Sz0

      i = states(l)
      call decode(i,nspin,config)

      m = m+1
      do k = 1, nspin-1
        H(m) = H(m) + Vzz(k) * (1 - 2 * config(k)) * &
          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))
      enddo
      k = nspin
      H(m) = H(m) + hz(k) * (1 - 2 * config(k))

      ROWS(m) = l
      COLS(m) = l

      do k = 1, nspin-1

        if (config(k)/=config(k+1)) then
          m = m+1
          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
          H(m) = H(m) + Jxy(k) * 2 * (config(k) - config(k+1))**2
          n = binsearch(j,states)
          COLS(m) = l
          ROWS(m) = n
        endif
      enddo

    enddo
    n = m

  end subroutine buildSz0_SPARSE_HMBL


  subroutine buildSz0_HMBL(nspin, dim_Sz0, Jxy, Vzz, hz, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)

    integer :: config(nspin), states(dim_Sz0)
    integer (c_int) :: i, j, k, m, n, l, r, rflag
    
    H = 0
    call zero_mag_states(nspin, dim_Sz0, states)
    do l = 1, dim_Sz0

      i = states(l)
      call decode(i,nspin,config)

      do k = 1, nspin-1

        H(l,l) = H(l,l) + Vzz(k) * (1 - 2 * config(k)) * &
          & (1 - 2 * config(k+1)) + hz(k) * (1 - 2 * config(k))

        if (config(k)/=config(k+1)) then
          j = i + (1-2*config(k))*2**(k-1) + (1-2*config(k+1))*2**(k)
          r = binsearch(j,states)

          H(r,l) = H(r,l) + Jxy(k) * 2 * (config(k) - config(k+1))**2
        endif
      enddo
      k = nspin
      H(l,l) = H(l,l) + hz(k) * (1 - 2 * config(k))
    enddo

  end subroutine buildSz0_HMBL

  subroutine buildSz0_HMBL_NNN(nspin, dim_Sz0, Jxy, Vzz, Vzz2, hz, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1), Vzz2(nspin-2), hz(nspin)
    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)

    integer :: config(nspin), states(dim_Sz0)
    integer (c_int) :: i, j, k, m, n, l, r, rflag
    
    H = 0
    call zero_mag_states(nspin, dim_Sz0, states)
    do l = 1, dim_Sz0

      i = states(l)
      call decode(i,nspin,config)
      config = 1 - 2*config

      do k = 1, nspin-2

        H(l,l) = H(l,l) + Vzz(k) * config(k) * &
          & config(k+1) + Vzz2(k) * config(k) * config(k+2) + hz(k) * config(k)

        if (config(k)/=config(k+1)) then
          j = i + config(k)*2**(k-1) + config(k+1)*2**(k)
          r = binsearch(j,states)

          H(r,l) = H(r,l) + Jxy(k) * (config(k) - config(k+1))**2/2
        endif
      enddo

      k = nspin - 1
      H(l,l) = H(l,l) + Vzz(k) * config(k) * &
        & config(k+1) + hz(k) * config(k)

      if (config(k)/=config(k+1)) then
        j = i + config(k)*2**(k-1) + config(k+1)*2**(k)
        r = binsearch(j,states)

        H(r,l) = H(r,l) + Jxy(k) * (config(k) - config(k+1))**2/2
      endif

      k = nspin
      H(l,l) = H(l,l) + hz(k) * config(k)
    enddo

  end subroutine buildSz0_HMBL_NNN

  subroutine buildSz0_HMBL_LR(nspin, dim_Sz0, Jxy, Vzz, hz, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1,nspin), hz(nspin)
    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)

    integer :: config(nspin), states(dim_Sz0), spin(nspin), inverse(dimSpinHalf**nspin)
    integer (c_int) :: i, j, l, r, k, q 

    H = 0
    call basis_Sz0_inv(nspin, dim_Sz0, states, inverse)
    do l = 1, dim_Sz0

      i = states(l)
      call decode(i,nspin,config)
      spin = 1 - 2*config

      do k = 1, nspin-1

        H(l,l) = H(l,l) + hz(k) * spin(k)

        do q = k+1, nspin
          H(l,l) = H(l,l) + Vzz(k,q) * spin(k) * spin(q)
        enddo

        if (spin(k)/=spin(k+1)) then
          j = i + (1 - 2*config(k))*2**(k-1) + (1 - 2*config(k+1))*2**(k)
          r = inverse(j+1) 

          H(r,l) = H(r,l) + Jxy(k) * (spin(k) - spin(k+1))**2/2
        endif
      enddo
      k = nspin
      H(l,l) = H(l,l) + hz(k) * spin(k)
    enddo

  end subroutine


  subroutine buildSz0_HSwap(nspin, dim_Sz0, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)

    integer :: config(nspin), states(dim_Sz0)
    integer (c_int) :: i, j, k, m, n, l, r, rflag
 
    H = 0
    call zero_mag_states(nspin, dim_Sz0, states)
    do l = 1, dim_Sz0

      i = states(l)
      call decode(i,nspin,config)

      do k = 1, nspin/2

        H(l,l) = H(l,l) + (1 - 2 * config(2*k-1)) * &
          & (1 - 2 * config(2*k))-1

        if (config(2*k-1)/=config(2*k)) then
          j = i + (1-2*config(2*k-1))*2**(2*k-2) + (1-2*config(2*k))*2**(2*k-1)
          r = binsearch(j,states)
          

          H(r,l) = H(r,l) + 2 * (config(2*k-1) - config(2*k))**2
        endif
      enddo
    enddo

  end subroutine buildSz0_HSwap

  subroutine buildHFlip(nspin, dim, H)

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(out) :: H(dim,dim)

    integer :: config(nspin)
    integer (c_int) :: i, j, k, m

    H = 0
    do i = 0, dim - 1

      call decode(i,nspin,config)

      do k = 1, nspin

        j = i + (1-2*config(k))*2**(k-1)
        H(j+1,i+1) = H(j+1,i+1) + 1

      enddo
    enddo

  end subroutine

  subroutine buildHKhemani( nspin, dim, Jzz, hx, hy, hz, H )

    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: Jzz(nspin-1), hx(nspin), hy(nspin), hz(nspin)
    complex (c_double_complex), intent(out) :: H(dim,dim)

    integer :: config(nspin), spin(nspin)
    integer (c_int) :: i, j, k, m

    H = 0
    do i = 0, dim - 1

      call decode(i, nspin, config)
      spin = 1 - 2 * config

      do k = 1, nspin - 1 !Open Boundary Conditions

        H(i+1,i+1) = H(i+1,i+1) + Jzz(k) * spin(k) * spin(k+1) + &
        & hz(k) * spin(k)

        j = i + (1-2*config(k)) * 2**(k-1)
        H(j+1,i+1) = H(j+1,i+1) + hx(k) + hy(k) * (C_UNIT * spin(k) )

      enddo
      k = nspin
      H(i+1,i+1) = H(i+1,i+1) + hz(k) * spin(k)
      j = i + (1-2*config(k)) * 2**(k-1)
      H(j+1,i+1) = H(j+1,i+1) + hx(k) + hy(k) * (-C_UNIT * spin(k) )
    enddo
  end subroutine

  subroutine buildSz_HSwap(nspin, dim_Sz, Sz, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz, Sz
    real (c_double), intent(out) :: H(dim_Sz,dim_Sz)

    integer :: config(nspin), states(dim_Sz), inverse(dimSpinHalf**nspin)
    integer (c_int) :: i, j, k, m, n, l, r, rflag
 
    H = 0
    call basis_Sz_inv(nspin, dim_Sz, Sz, states, inverse)
    do l = 1, dim_Sz

      i = states(l)
      call decode(i,nspin,config)

      do k = 1, nspin/2

        H(l,l) = H(l,l) + (1 - 2 * config(2*k-1)) * &
          & (1 - 2 * config(2*k))-1

        if (config(2*k-1)/=config(2*k)) then
          j = i + (1-2*config(2*k-1))*2**(2*k-2) + (1-2*config(2*k))*2**(2*k-1)
          r = inverse(j+1)
          

          H(r,l) = H(r,l) + 2 * (config(2*k-1) - config(2*k))**2
        endif
      enddo
    enddo

  end subroutine buildSz_HSwap

  subroutine buildSz_HMBL_LR(nspin, dim_Sz, Sz, Jxy, Vzz, hz, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz, Sz
    real (c_double), intent(in) :: Jxy(nspin-1), Vzz(nspin-1,nspin), hz(nspin)
    real (c_double), intent(out) :: H(dim_Sz,dim_Sz)

    integer :: config(nspin), states(dim_Sz), spin(nspin), inverse(dimSpinHalf**nspin)
    integer (c_int) :: i, j, l, r, k, q 

    H = 0
    call basis_Sz_inv(nspin, dim_Sz, Sz, states, inverse)
    do l = 1, dim_Sz

      i = states(l)
      call decode(i,nspin,config)
      spin = 1 - 2*config

      do k = 1, nspin-1

        H(l,l) = H(l,l) + hz(k) * spin(k)

        do q = k+1, nspin
          H(l,l) = H(l,l) + Vzz(k,q) * spin(k) * spin(q)
        enddo

        if (spin(k)/=spin(k+1)) then
          j = i + (1 - 2*config(k))*2**(k-1) + (1 - 2*config(k+1))*2**(k)
          r = inverse(j+1) 

          H(r,l) = H(r,l) + Jxy(k) * (spin(k) - spin(k+1))**2/2
        endif
      enddo
      k = nspin
      H(l,l) = H(l,l) + hz(k) * spin(k)
    enddo

  end subroutine

  subroutine print_hamiltonian_Sz(nspin, dim_Sz, Sz, H)

    integer (c_int), intent(in) :: nspin, dim_Sz, Sz
    real (c_double_complex), intent(in) :: H(dim_Sz, dim_Sz)

    integer :: i, j, r, l, config(nspin), config2(nspin)
    integer :: idxSz(dim_Sz)

    call basis_Sz(nspin, dim_Sz, Sz, idxSz)

    print "(4X,1(A6,X),2(A3,3X),A4)", "H(r,l)", "i", "r/l", "conf"
    do l = 1, dim_Sz
      i = idxSz(l)
      call decode(i,nspin,config)

      if(abs(H(l,l)) > tol) then
        print "(4X,1(F6.2,X),2(I3,3X),*(I0))", H(l,l), l, i, config(:)
      endif

      do r = 1, dim_Sz

        if(r==l) then 
          cycle
        else if (abs(H(r,l)) > tol) then
          j = idxSz(r)
          call decode(j,nspin, config2)
          print "(4X,1(F6.2,X),2(I3,3X),*(I0))", H(r,l), l, i, config(:)
          print "(4X,1(A6,X),2(I3,3X),*(I0))", "", r, j, config2(:)
        endif

      enddo

    enddo
    print *, ""


  end subroutine

  subroutine print_unitary_Sz(nspin, dim_Sz, Sz, U)

    integer (c_int), intent(in) :: nspin, dim_Sz, Sz
    complex (c_double_complex), intent(in) :: U(dim_Sz, dim_Sz)

    integer :: i, j, r, l, config(nspin), config2(nspin)
    integer :: idxSz(dim_Sz)

    call basis_Sz(nspin, dim_Sz, Sz, idxSz)

    print "(4X,1(A18,X),2(A3,3X),A4)", "U(r,l)", "i", "r/l", "conf"
    do l = 1, dim_Sz
      i = idxSz(l)
      call decode(i,nspin,config)

      if(abs(U(l,l)) > tol) then
        print "(4X,1((f8.2f8.2x'i'),X),2(I3,3X),*(I0))", U(l,l), l, i, config(:)
      endif

      do r = 1, dim_Sz

        if(r==l) then 
          cycle
        else if (abs(U(r,l)) > tol) then
          j = idxSz(r)
          call decode(j,nspin, config2)
          print "(4X,1(f8.2f8.2x'i',X),2(I3,3X),*(I0))", U(r,l), l, i, config(:)
          print "(4X,1(A18,X),2(I3,3X),*(I0))", "", r, j, config2(:)
        endif

      enddo

    enddo
    print *, ""


  end subroutine

end module matrices
