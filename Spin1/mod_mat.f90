!This module contains all procedures regarding the generation of matrices
!Spefically, it contains the procedures for spin hamiltonains used in time evolution:
! in the Full Hilbert space, Sz=0 subspace, generic Sz subspace
! both dense and sparse versions


module matrices

  !use ifport
  use iso_c_binding
  use printing
  use functions
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
  subroutine buildSz0_HMBL(nspin, dim_Sz0, Jint, Vint, hz, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Jint(nspin-1), Vint(nspin-1), hz(nspin)
    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)

    integer :: config(nspin), states(dim_Sz0), config2(nspin)
    integer :: inverse(dimSpin1**nspin)
    integer (c_int) :: i, j, k, m, n, l, r
    integer (c_int) :: idx1, idx2, idxp1, idxp2, s1, s2
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
    call zero_mag_states_inv(nspin, dim_Sz0, states, inverse)
    !print "(4X,1(A6,X),3(A3,3X),A4)", "H(r,l)", "k", "i", "r/l", "conf"
    do l = 1, dim_Sz0

      i = states(l)
      call decode(i,nspin,config)

      do k = 1, nspin-1

        s1 = 1 - config(k)
        s2 = 1 - config(k+1)
        H(l,l) = H(l,l) + Vint(k) * s1*s2 + hz(k) * s1
        !print "(4X,1(F6.2,X),3(I3,3X),*(I0))", H(l,l), k, i, l, config(:)

        do idxp1 = 0, dimSpin1-1
          do idxp2 = 0, dimSpin1-1

            idx1 = config(k)
            idx2 = config(k+1)
            if( (idxp1 + idxp2 == idx1 + idx2) .AND. (idxp1 /= idx1) .AND. (idxp2 /= idx2) ) then
              j = i + (idxp1 - idx1) * dimSpin1**(k-1) + (idxp2 - idx2) * dimSpin1**k
              r = inverse(j+1) 
              H(r,l) = H(r,l) + Jint(k) * OPRxy(idxp1+1, idxp2+1, idx1+1, idx2+1)
              call decode(j,nspin, config2)
              !print "(4X,1(F6.2,X),3(I3,3X),*(I0))", H(r,l), k, i, l, config(:)
              !print "(4X,1(A6,X),3(I3,3X),*(I0))", "", k, j, r, config2(:)
            endif

          enddo
          !idx1 = config(k)
          !j = i + (idxp1 - idx1) * dimSpin1**(k-1)
          !H(j+1,i+1) = H(j+1,i+1) + hz(k) * OPRz(idxp1+1, idx1+1)

        enddo
      enddo
      k = nspin
      s1 = 1 - config(k)
      H(l,l) = H(l,l) + hz(k) * s1
      !print "(4X,1(F5.1,X),3(I3,3X),*(I0))", H(l,l), k, i, l, config(:)
    enddo

  end subroutine buildSz0_HMBL
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
  subroutine buildSz0_HSwap(nspin, dim_Sz0, H)

    !dim_Sz0 is the dimensione of the Sz=0 subsapce

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(out) :: H(dim_Sz0,dim_Sz0)

    integer :: config(nspin), states(dim_Sz0), config2(nspin), inverse(dimSpin1**nspin)
    integer (c_int) :: i, j, k, m, n, l, r
    integer (c_int) :: idx1, idx2, idxp1, idxp2, s1, s2, idxs1, idxs2
    real (c_double) :: OPR1(dimSpin1,dimSpin1,dimSpin1,dimSpin1)
    real (c_double) :: OPR2(dimSpin1,dimSpin1,dimSpin1,dimSpin1)
    real (c_double) :: OPR(dimSpin1,dimSpin1,dimSpin1,dimSpin1)

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
             !& KDelta(idxp1,idx1) * KDelta(idxp2,idx2)


          enddo
        enddo
      enddo
    enddo
    
    H = 0
    call zero_mag_states_inv(nspin, dim_Sz0, states, inverse)
    !print "(4X,1(A6,X),3(A3,3X),A4)", "H(r,l)", "k", "i", "r/l", "conf"
    do l = 1, dim_Sz0

      i = states(l)
      call decode(i,nspin,config)

      do k = 1, nspin/2

        s1 = 1 - config(2*k-1)
        s2 = 1 - config(2*k)
        do idxp1 = 0, dimSpin1-1
          do idxp2 = 0, dimSpin1-1

            idx1 = config(2*k-1)
            idx2 = config(2*k)
            if( (idxp1 + idxp2 == idx1 + idx2) ) then
              j = i + (idxp1 - idx1) * dimSpin1**(2*k-2) + (idxp2 - idx2) * dimSpin1**(2*k-1)
              r = inverse(j+1) 
              H(r,l) = H(r,l) + OPR(idxp1+1, idxp2+1, idx1+1, idx2+1)
              call decode(j,nspin, config2)
              !print "(4X,1(F6.2,X),3(I3,3X),*(I0))", H(r,l), k, i, l, config(:)
              !print "(4X,1(A6,X),3(I3,3X),*(I0))", "", k, j, r, config2(:)
            endif

          enddo
        enddo
      enddo
    enddo

  end subroutine buildSz0_HSwap


  subroutine print_hamiltonian_Sz0(nspin, dim_Sz0, H)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double_complex), intent(in) :: H(dim_Sz0, dim_Sz0)

    integer :: i, j, r, l, config(nspin), config2(nspin)
    integer :: idxSz0(dim_Sz0)

    call zero_mag_states(nspin, dim_Sz0, idxSz0)

    print "(4X,1(A6,X),2(A3,3X),A4)", "H(r,l)", "i", "r/l", "conf"
    do l = 1, dim_Sz0
      i = idxSz0(l)
      call decode(i,nspin,config)

      if(abs(H(l,l)) > tol) then
        print "(4X,1(F6.2,X),2(I3,3X),*(I0))", H(l,l), l, i, config(:)
      endif

      do r = 1, dim_Sz0

        if(r==l) then 
          cycle
        else if (abs(H(r,l)) > tol) then
          j = idxSz0(r)
          call decode(j,nspin, config2)
          print "(4X,1(F6.2,X),2(I3,3X),*(I0))", H(r,l), l, i, config(:)
          print "(4X,1(A6,X),2(I3,3X),*(I0))", "", r, j, config2(:)
        endif

      enddo

    enddo
    print *, ""


  end subroutine


  subroutine print_unitary_Sz0(nspin, dim_Sz0, U)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: U(dim_Sz0, dim_Sz0)

    integer :: i, j, r, l, config(nspin), config2(nspin)
    integer :: idxSz0(dim_Sz0)

    call zero_mag_states(nspin, dim_Sz0, idxSz0)

    print "(4X,1(A18,X),2(A3,3X),A4)", "U(r,l)", "i", "r/l", "conf"
    do l = 1, dim_Sz0
      i = idxSz0(l)
      call decode(i,nspin,config)

      if(abs(U(l,l)) > tol) then
        print "(4X,1((f8.2f8.2x'i'),X),2(I3,3X),*(I0))", U(l,l), l, i, config(:)
      endif

      do r = 1, dim_Sz0

        if(r==l) then 
          cycle
        else if (abs(U(r,l)) > tol) then
          j = idxSz0(r)
          call decode(j,nspin, config2)
          print "(4X,1(f8.2f8.2x'i',X),2(I3,3X),*(I0))", U(r,l), l, i, config(:)
          print "(4X,1(A18,X),2(I3,3X),*(I0))", "", r, j, config2(:)
        endif

      enddo

    enddo
    print *, ""


  end subroutine

end module matrices
