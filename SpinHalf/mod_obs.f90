!This module contains various procedures for the computation of simple observables:
! Local Magnetization < sigma_k^z >
! Imbalance I = 1/L * sum_k (-1)^k < sigma_k^z >
! Local Imbalance LI = 1/(L/2) * sum_k < (sigma_(2k)^z - sigma_(2k-1)^z)^2 >
! Correlation function < sigma_k^z * sigma_q^z > - < sigma_k^z > * < sigma_q^z >
! Time average of observable with its error:
!             < O > = sum_i O_i / N,  sigma(O) = < O^2 > - (< O >)^2
!   Options: start from non-zero time step, average over (-1)^i O_i to verify TTSB
! Exact energies of integrable model E({s_k}) = sum_k h_k s_k + V_k s_k s_{k+1}
!     and exact quasi-energies QE({s_k}) = E({s_k}) + E({s_k<->swap)}

! There are variants for both full Hilbert space and Sz subspaces (specifically Sz=0)

module observables
  
  use iso_c_binding, dp => c_double, ip => c_int, cp => c_double_complex
  use printing
  use functions, search => binsearch_closest_from_above
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  integer (c_int), private, parameter :: dimSpin1 = 3
  real (c_double), parameter, private :: pi = 4.d0 * datan(1.d0)

contains

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
      enddo
      magaux = magaux * abs(state(i))**2
      mag = mag + magaux
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
      enddo
      magaux = magaux * abs(state(i))**2
      mag = mag + magaux
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
      enddo
      magaux = magaux * abs(state(i))**2
      imbaux = imbaux * abs(state(i))**2
      mag = mag + magaux
      imb = imb + imbaux
    enddo
    mag = mag/nspin
    imb = imb/nspin
    imb = imb/(1+mag)
    imbalance = imb


  end function imbalance

  real function imbalance_Sz0(nspin, dim_Sz0, state)

    integer(c_int), intent(in) :: nspin, dim_Sz0
    complex(c_double_complex), intent(in) :: state(dim_Sz0)
    real(c_double) :: imb, imbaux
    integer :: config(nspin)
    integer (c_int) :: i, j, k, m, l, indx(dim_Sz0)

    imb = 0
    call zero_mag_states(nspin, dim_Sz0, indx)

    do i = 1, dim_Sz0

      l = indx(i)
      call decode(l,nspin,config)
      imbaux = 0
      do k = 1, nspin
        imbaux = imbaux + (-1)**k * (1 - 2 * config(k))
      enddo
      imbaux = imbaux * abs(state(i))**2
      imb = imb + imbaux
    enddo
    imb = imb/nspin
    imbalance_Sz0 = imb

  end function imbalance_Sz0

  real function imbalance_sq_Sz0(nspin, dim_Sz0, state)

    integer(c_int), intent(in) :: nspin, dim_Sz0
    complex(c_double_complex), intent(in) :: state(dim_Sz0)
    real(c_double) :: imb, imbaux
    integer :: config(nspin)
    integer (c_int) :: i, k1, k2, l, indx(dim_Sz0)

    imb = 0
    call zero_mag_states(nspin, dim_Sz0, indx)

    do i = 1, dim_Sz0

      l = indx(i)
      call decode(l,nspin,config)
      imbaux = 0
      do k1 = 1, nspin
        do k2 = 1, nspin
          imbaux = imbaux + (-1)**k1 * (1 - 2 * config(k1)) * (-1)**k2 * (1 - 2 * config(k2))
        enddo
      enddo
      imbaux = imbaux * abs(state(i))**2
      imb = imb + imbaux
    enddo
    imb = imb/nspin**2
    imbalance_sq_Sz0 = imb

  end function imbalance_sq_Sz0


  function local_imbalance(nspin, dim, psi) result(LI)
    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: psi(dim)

    integer (c_int) :: i, k, config(nspin)
    real (c_double) :: LI, LI_part

    LI = 0
    do i = 1, dim
    call decode(i-1, nspin, config)
      LI_part = 0
      do k = 1, nspin/2
        LI_part = LI_part + (config(2*k) - config(2*k-1))**2
      enddo
      LI = LI + abs(psi(i))**2 * 2 * LI_part / nspin
    enddo

  end function

  function local_imbalance_Sz0(nspin, dim_Sz0, psi_Sz0) result(LI)
    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)

    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)
    real (c_double) :: LI, LI_part

    call zero_mag_states(nspin, dim_Sz0, states)
    LI = 0
    do i = 1, dim_Sz0
      l = states(i)
      call decode(l, nspin, config)
      LI_part = 0
      do k = 1, nspin/2
        LI_part = LI_part + (config(2*k) - config(2*k-1))**2
      enddo
      LI = LI + abs(psi_Sz0(i))**2 * 2 * LI_part / nspin
    enddo

  end function

  function local_overlap(nspin, dim, psi) result(LO)
    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: psi(dim)

    integer (c_int) :: i, j, k, config(nspin)
    real (c_double) :: LO

    LO = 0
    do i = 1, dim
      call decode(i-1, nspin, config)
      do k = 1, nspin/2
        if (config(2*k)/=config(2*k-1)) then
          j = i + (config(2*k-1) - config(2*k)) * (2**(2*k-1) - 2**(2*k-2))
          LO = LO + abs(psi(i))**2 - psi(i) * dconjg(psi(j))
        endif
      enddo
    enddo
    LO = LO/(nspin/2)

  end function

  function local_overlap_Sz0(nspin, dim_Sz0, psi) result(LO)
    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: psi(dim_Sz0)

    integer (c_int) :: i, j, k, r, l, config(nspin), states(dim_Sz0)
    real (c_double) :: LO

    call zero_mag_states(nspin, dim_Sz0, states)
    LO = 0
    do i = 1, dim_Sz0
      l = states(i)
      call decode(l, nspin, config)
      do k = 1, nspin/2
        if(config(2*k)/=config(2*k-1)) then
          j = i + (config(2*k-1) - config(2*k)) * (2**(2*k-1) - 2**(2*k-2))
          r = binsearch(j,states)
          LO = LO + abs(psi(i))**2 - psi(i) * dconjg(psi(r))
        endif
      enddo
    enddo
    LO = LO/(nspin/2)

  end function

  function sigmaz_corr_c(nspin, dim, q, p, psi)

    integer (c_int), intent(in) :: nspin, dim, q, p
    complex (c_double_complex), intent(in) :: psi(dim)
    real (c_double) :: sigmaz_corr_c

    integer (c_int) :: i, k, config(nspin)
    real (c_double) :: corr, avgq, avgp

    avgq = 0
    avgp = 0
    corr = 0
    do i = 1, dim
      
      call decode(i-1,nspin,config)

      corr = corr + abs(psi(i))**2 * (1 - 2 * config(q)) * (1 - 2 * config(p))
      avgq = avgq + abs(psi(i))**2 * (1 - 2 * config(q))
      avgp = avgp + abs(psi(i))**2 * (1 - 2 * config(p))
        
    enddo

    sigmaz_corr_c = corr - avgq * avgp


  end function

  function sigmaz_corr_c_Sz0(nspin, dim_Sz0, q, p, psi_Sz0)

    integer (c_int), intent(in) :: nspin, dim_Sz0, q, p
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: sigmaz_corr_c_Sz0

    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)
    real (c_double) :: corr, avgq, avgp

    avgq = 0
    avgp = 0
    corr = 0
    call zero_mag_states(nspin, dim_Sz0, states)
    do i = 1, dim_Sz0
      
      l = states(i)
      call decode(l,nspin,config)

      corr = corr + abs(psi_Sz0(i))**2 * (1 - 2 * config(q)) * (1 - 2 * config(p))
      avgq = avgq + abs(psi_Sz0(i))**2 * (1 - 2 * config(q))
      avgp = avgp + abs(psi_Sz0(i))**2 * (1 - 2 * config(p))
        
    enddo

    sigmaz_corr_c_Sz0 = corr - avgq * avgp

  end function

  function sigmaz_corr_matrix_Sz0(nspin, dim_Sz0, psi_Sz0) result(CORR)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: CORR(nspin, nspin)

    integer :: q, p

    do q = 1, nspin
      do p = 1, nspin
        CORR(q,p) = sigmaz_corr_c_Sz0(nspin, dim_Sz0, q, p, psi_Sz0)
      enddo
    enddo
    
  end function

  function sigmaz_tot_corr_Sz0(nspin, dim_Sz0, psi_Sz0) result(tot_corr)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: tot_corr

    real (c_double) :: CORR(nspin, nspin)
    integer :: q, p

    CORR = sigmaz_corr_matrix_Sz0(nspin, dim_Sz0, psi_Sz0)
    
    tot_corr = 0
    do q = 1, nspin
      do p = 1, nspin
        tot_corr = tot_corr + abs(CORR(q,p))
      enddo
    enddo

  end function


  subroutine local_zmag_Sz0(nspin, dim_Sz0, psi_Sz0, sigmaz)
    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double), intent(out) :: sigmaz(nspin)

    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)

    call zero_mag_states(nspin, dim_Sz0, states)
    sigmaz = 0
    do i = 1, dim_Sz0
      l = states(i)
      call decode(l, nspin, config)
      do k = 1, nspin
        sigmaz(k) = sigmaz(k) + abs(psi_Sz0(i))**2 * (1-2*config(k))
      enddo
    enddo

  end subroutine

  subroutine gap_ratio(dim, energies, r_avg, r_sq)

    integer, intent(in) :: dim
    real (c_double), intent(in) :: energies(dim)
    real (c_double), intent(out) :: r_avg, r_sq

    integer :: n
    real (c_double) :: gap1, gap2

    r_avg = 0
    r_sq = 0
    do n = 1, dim - 2

      gap1 = energies(n+1) - energies(n)
      gap2 = energies(n+2) - energies(n+1)

      r_avg = r_avg + min(gap1,gap2)/max(gap1,gap2)
      r_sq = r_sq + (min(gap1,gap2)/max(gap1,gap2))**2

    enddo

    r_avg = r_avg / (dim-2)
    r_sq = r_sq / (dim-2)

  end subroutine gap_ratio

  function pi_pair(dim, QE) result(beta)

    integer (c_int), intent(in) :: dim
    real (c_double), intent(in) :: QE(dim)
    integer (c_int) :: alpha, beta(dim)
    real (c_double) :: val

    do alpha = 1, dim
      val = mod(QE(alpha) + 2*pi, 2*pi) - pi
      if (val > QE(dim)) then 
        val = val - 2*pi !If (QE(alpha)+pi)_1 > QE(beta) for all beta, then look from below by shifting by -2pi
        !print *, "Problem: alpha = ", alpha, "; a(alpha) = ", QE(alpha)
        !print *, ""
      endif

      beta(alpha) = search(val, QE)
      !print "(2(A12),4(A22,4X))", "i", "alpha", "a(i)", "val", "(a(alpha) + pi)_1", "a(alpha)"
      !do i = 1, size(QE)
      !  print *, i, alpha, QE(i), val, mod(QE(alpha) + 2*pi, 2*pi) - pi, QE(alpha)
      !enddo

      !print "(2(A12),4(A22,4X))", "beta", "alpha", "a(beta)", "val", "(a(alpha) + pi)_1", "a(alpha)"
      !print *, beta, alpha, QE(beta), val, mod(QE(alpha) + 2*pi, 2*pi) - pi, QE(alpha)

      !print *, "Check: "
      !print "(1(A12),1(A22,4X))", "i", "a(i) - val"!, "|a(i) - a(alpha)|-pi"
      !do i = 1, size(QE)
      !  print *, i, QE(i) - val!, abs(QE(i) - QE(alpha)) - pi
      !enddo
      !print *, ""
      !print *, ""
      !print *, ""

    enddo

  end function

  subroutine log_gap_difference(dim, energies, log_pair_avg, log_near_avg, log_avg, log_sq)

    !Computes the averages of the logarithm of gaps of neighbouring eigenvalues 
    !and paired eigenvalues (separated by half the spectrum)
    ! log_near_avg = < log( E_{n+1} - E_n ) >
    ! log_pair_avg = < log( E_{n+dim/2} - E_n ) >
    ! log_avg = log_pair_avg - log_near_avg
    ! log_sq = < (log Delta^alpha - log Delta_0^alpha)^2 >


    integer (c_int), intent(in) :: dim
    real (c_double), intent(in) :: energies(dim)
    real (c_double), intent(out) :: log_avg, log_pair_avg, log_near_avg, log_sq
    real (c_double) :: pair(dim), near(dim)

    integer (c_int) :: alpha, beta, pi_paired(dim)

    pi_paired = pi_pair(dim, energies)

    do alpha = 1, dim 

      !beta = merge(alpha+1,alpha+1-dim,alpha<=dim-1)
      !near(alpha) = energies(beta) - energies(alpha)
      
      if (alpha<=dim-1) then
        beta = alpha + 1
        near(alpha) = energies(beta) - energies(alpha)
      else if (alpha>dim-1) then
        beta = alpha + 1 - dim
        near(alpha) = energies(beta) + 2*pi - energies(alpha)
      endif

      !beta = merge(alpha+dim/2, alpha-dim/2, alpha+dim/2<=dim)
      beta = pi_paired(alpha)
      pair(alpha) = abs(abs(energies(beta) - energies(alpha)) - pi)

    enddo

    log_pair_avg = sum(log(pair)) / dim
    log_near_avg = sum(log(near)) / dim
    log_avg = log_pair_avg - log_near_avg
    log_sq = sum((log(pair) - log(near))**2) / dim

  end subroutine log_gap_difference

  function IPR(psi)
    complex (c_double_complex), intent(in) :: psi(:)
    real (c_double) :: IPR
    integer (c_int) :: dim, i

    dim = size(psi)
    !print *, dim
    IPR = 0
    do i = 1, dim
      IPR = IPR + abs(psi(i))**4
      !print*, IPR, abs(psi(i))
    enddo

  end function IPR

  subroutine finite_imbalance_states(nspin, dim, IMB, states)
    integer (c_int), intent(in) :: nspin, dim
    real (c_double), intent(in) :: IMB
    integer (c_int), intent(out) :: states(dim)

    integer (c_int) :: i, k, l, config(nspin), sum_even, sum_odd
    real (c_double) :: loc_imb

    l = 0
    print *, "       i     l  config"
    do i = 0, dim-1
      call decode(i, nspin, config)
      sum_odd = 0
      sum_even = 0
      !print *, i, sum_even, sum_odd
      do k = 1, nspin/2
        sum_odd = sum_odd + config(2*k-1)
        sum_even = sum_even + config(2*k)
        !print *, k, sum_even, sum_odd
      enddo
      loc_imb = real(sum_odd - sum_even, c_double) / real(nspin - sum_odd - sum_even, c_double)
      !print *, loc_imb
      if(abs(IMB) < 1.0e-10 .AND. sum_odd == sum_even) then
        l = l+1
        states(l) = i
        print "(2(8X,I0),4X,*(I0))", i, l, config(:)
      else if(abs(IMB) > 1.0e-10 .AND. abs(loc_imb - IMB) < 1.0e-10) then
        l = l+1
        states(l) = i
        print "(2(8X,I0),4X,*(I0))", i, l, config(:)
      endif
    enddo

  end subroutine finite_imbalance_states

  subroutine large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, states)
    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: IMB, LI
    integer (c_int), intent(out) :: states(dim_Sz0)

    integer :: i, l, m, idxSz0(dim_Sz0)
    real (c_double) :: IMBc, LIc

    call zero_mag_states(nspin, dim_Sz0, idxSz0)
    m = 0
    states = 0
    !print "( 4X,A4, 4X,A4, 4X,A6, 4X,A6 )", "m", "l", "IMB", "LI"
    do i = 1, dim_Sz0

      l = idxSz0(i)
      IMBc = imbalance_basis(nspin, l)
      LIc = local_imbalance_basis(nspin, l)
      if (IMB >= 0) then
        if (IMBc >= IMB .AND. LIc >= LI) then
          m = m+1
          states(m) = l
          !print "(4X,I4, 4X,I4, 4X,F6.3, 4X,F6.3, 4X,F6.3)", m, l, IMBc, LIc
        endif
      else if (IMB < 0) then
        if (IMBc <= IMB .AND. LIc >= LI) then
          m = m+1
          states(m) = l
          !print "(4X,I4, 4X,I4, 4X,F6.3, 4X,F6.3)", m, l, IMBc, LIc
        endif
      endif
    enddo

    !print *, "states = ", states(1:m)

  end subroutine

  subroutine finite_imbalance_states_Sz0(nspin, dim_Sz0, IMB, states)
    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: IMB
    integer (c_int), intent(out) :: states(dim_Sz0)

    integer (c_int) :: i, k, m, l, config(nspin), sum_even, sum_odd, idxSz0(dim_Sz0)
    real (c_double) :: loc_imb

    m = 0
    call zero_mag_states(nspin, dim_Sz0, idxSz0)
    print *, "       i     l  config"
    do i = 1, dim_Sz0
      l = idxSz0(i)
      call decode(l, nspin, config)
      sum_odd = 0
      sum_even = 0
      !print *, i, sum_even, sum_odd
      do k = 1, nspin/2
        sum_odd = sum_odd + config(2*k-1)
        sum_even = sum_even + config(2*k)
        !print *, k, sum_even, sum_odd
      enddo
      loc_imb = real(sum_odd - sum_even, c_double) / real(nspin - sum_odd - sum_even, c_double)
      !print *, loc_imb
      if(abs(IMB) < 1.0e-10 .AND. sum_odd == sum_even) then
        m = m+1
        states(m) = l
        print "(2(8X,I0),4X,*(I0))", l, m, config(:)
      else if(abs(IMB) > 1.0e-10 .AND. abs(loc_imb - IMB) < 1.0e-10) then
        m = m+1
        states(m) = l
        print "(2(8X,I0),4X,*(I0))", l, m, config(:)
      endif
    enddo

  end subroutine


  function imbalance_basis(nspin, i)

    integer(c_int), intent(in) :: nspin, i
    real (c_double) :: imb, imbalance_basis
    integer :: config(nspin), k

    call decode(i,nspin,config)
    imb = 0
    do k = 1, nspin
        imb = imb + (-1)**k * (1 - 2 * config(k))
        !print *, k, imb
    enddo
    imb = imb/(2*nspin - 2*sum(config))
    imbalance_basis = imb

  end function imbalance_basis

  function imbalance_sq_basis(nspin, i)

    integer(c_int), intent(in) :: nspin, i
    real (c_double) :: imb, imbalance_sq_basis
    integer :: config(nspin), k

    call decode(i,nspin,config)
    imb = 0
    do k = 1, nspin
        imb = imb + (-1)**k * (1 - 2 * config(k))
        !print *, k, imb
    enddo
    imb = imb/(2*nspin - 2*sum(config))
    imbalance_sq_basis = imb**2

  end function imbalance_sq_basis

  function local_imbalance_basis(nspin, i)
    integer (c_int), intent(in) :: nspin, i
    real (c_double) :: local_imbalance_basis

    integer (c_int) :: k, config(nspin)
    real (c_double) :: LI

    LI = 0
    call decode(i, nspin, config)
    do k = 1, nspin/2

      LI = LI + (config(2*k) - config(2*k-1))**2
    enddo
    LI = 2 * LI / nspin
    local_imbalance_basis = LI

  end function

  subroutine buildLarge_IMB_LI_State_basis_Sz0(nspin, dim_Sz0, i, IMB, LI, psi)
    integer (c_int), intent(in) :: nspin, dim_Sz0, i
    real (c_double), intent(in) :: IMB, LI
    complex (c_double_complex), intent(out) :: psi(dim_Sz0)

    integer (c_int) :: idxSz0(dim_Sz0), nz, j, l, config(nspin)
    integer (c_int), allocatable :: idx(:)
    real (c_double) :: rand

    call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, idxSz0)

    do j = 1, dim_Sz0
      if(idxSz0(j) == 0) then
        nz = j-1
        exit
      endif
    enddo
    print *, "nz = ", nz, "i = ", i
    if (i>nz) stop "Error: Not enough number of states."

    l = idxSz0(i)
    call zero_mag_states(nspin, dim_Sz0, idxSz0)
    j = binsearch(l, idxSz0)
    psi = 0
    psi(j) = 1

    call decode(l, nspin, config)
    print *, "build Large IMB, LI basis state"
    print *, IMB, LI
    print *, imbalance_basis(nspin,l), local_imbalance_basis(nspin,l), l
    1 format ( 4X,F6.3, 4X,F6.3, 4X,I4, 4X,I0 )

  end subroutine


  function exact_energy(nspin, Vzz, hz, i)

    integer (c_int), intent(in) :: nspin, i
    real (c_double), intent(in) :: Vzz(nspin), hz(nspin)
    real (c_double) :: exact_energy

    integer (c_int) :: k, config(nspin)
    real (c_double) :: E

    call decode(i, nspin, config)
    config = 1 - 2*config
    E = 0
    do k = 1, nspin - 1
      E = E + Vzz(k) * config(k) * config(k+1) + hz(k) * config(k)
    enddo
    k = nspin
    E = E + hz(k) * config(k)

    exact_energy = E

  end function

  subroutine exact_energies_Sz0(nspin, dim_Sz0, Vzz, hz, E)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Vzz(nspin), hz(nspin)
    real (c_double), intent(out) :: E(dim_Sz0)

    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)

    call zero_mag_states(nspin, dim_Sz0, states)
    do i = 1, dim_Sz0

      l = states(i)
      call decode(l, nspin, config)
      E(i) = exact_energy(nspin, Vzz, hz, l)

    enddo

  end subroutine

  function exact_quasi_energy(nspin, Vzz, hz, i) result(QE)

    integer (c_int), intent(in) :: nspin, i
    real (c_double), intent(in) :: Vzz(nspin-1), hz(nspin)

    integer (c_int) :: k, m, config(nspin)
    real (c_double) :: QE

    call decode(i, nspin, config)
    m = 0
    do k = 1, nspin/2
      m = m + 2**(2*k-2) * config(2*k) + 2**(2*k-1) * config(2*k-1)
    enddo
    QE = (exact_energy(nspin, Vzz, hz, i) + exact_energy(nspin, Vzz, hz, m)) / 2

  end function


  subroutine exact_quasi_energies_Sz0(nspin, dim_Sz0, Vzz, hz, QE) !E, Es, QE, QE_alt)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Vzz(nspin-1), hz(nspin)
    real (c_double), intent(out) :: QE(dim_Sz0)

    integer (c_int) :: i, k, l, m, config(nspin), states(dim_Sz0)
    !real (c_double) :: E(dim_Sz0), Es(dim_Sz0), QE_alt(dim_Sz0)

    QE = 0
    !QE_alt = 0
    !E = 0
    !Es = 0
    call zero_mag_states(nspin, dim_Sz0, states)
    do i = 1, dim_Sz0

      l = states(i)
      call decode(l, nspin, config)
      m = 0
      do k = 1, nspin/2
        m = m + 2**(2*k-2) * config(2*k) + 2**(2*k-1) * config(2*k-1)
      enddo
      QE(i) = (exact_energy(nspin, Vzz, hz, l) + exact_energy(nspin, Vzz, hz, m)) / 2
      QE(i) = real(C_UNIT*log(exp(-C_UNIT*QE(i))), kind(1.0_c_double))

      write (*,"(*(I0))",advance='no'), config(:)
      write (*,*) QE(i)
      !config = 1 - 2*config
      !do k = 1, nspin/2 - 1
      !  print *, "k = ", k
      !  E(i) = E(i) + Vzz(2*k-1) * config(2*k-1) * config(2*k) + Vzz(2*k) * config(2*k) * config(2*k+1) + &
      !   &  hz(2*k-1) * config(2*k-1) + hz(2*k) * config(2*k) 
      !  Es(i) = Es(i) + Vzz(2*k-1) * config(2*k) * config(2*k-1) + Vzz(2*k) * config(2*k-1) * config(2*k+2) + &
      !   &  hz(2*k-1) * config(2*k) + hz(2*k) * config(2*k-1) 

      !  print *, "Normal contribution V"
      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", Vzz(2*k-1) * config(2*k-1) * config(2*k), &
      !    & Vzz(2*k) * config(2*k) * config(2*k+1), &
      !    & 2*k-1, config(2*k-1), 2*k, config(2*k), 2*k+1, config(2*k+1)
      !  print *, "Swap contribution V"
      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", Vzz(2*k-1) * config(2*k) * config(2*k-1), &
      !    & Vzz(2*k) * config(2*k-1) * config(2*k+2), &
      !    & 2*k, config(2*k), 2*k-1, config(2*k-1), 2*k+2, config(2*k+2)
      !  print *, "Normal contribution h"
      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", hz(2*k-1) * config(2*k-1), hz(2*k) * config(2*k), &
      !    & 2*k-1, config(2*k-1), 2*k, config(2*k)
      !  print *, "Swap contribution h"
      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", hz(2*k-1) * config(2*k), hz(2*k) * config(2*k-1), &
      !    & 2*k, config(2*k), 2*k-1, config(2*k-1)
      !  print *, ""

      !enddo
      !k = nspin/2
      !  print *, "k = ", k
      !E(i) = E(i) + Vzz(2*k-1) * config(2*k-1) * config(2*k) + &
      !  & + hz(2*k-1) * config(2*k-1) + hz(2*k) * config(2*k)
      !Es(i) = Es(i) + Vzz(2*k-1) * config(2*k) * config(2*k-1) + &
      !  & + hz(2*k-1) * config(2*k) + hz(2*k) * config(2*k-1)
      !  print *, "Normal contribution V"
      !  print "( F20.15,4X, *(2I0,2X)  )", Vzz(2*k-1) * config(2*k-1) * config(2*k), &
      !    & 2*k-1, config(2*k-1), 2*k, config(2*k)
      !  print *, "Swap contribution V"
      !  print "( F20.15,4X, *(2I0,2X)  )", Vzz(2*k-1) * config(2*k) * config(2*k-1), &
      !    & 2*k, config(2*k), 2*k-1, config(2*k-1)
      !  print *, "Normal contribution h"
      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", hz(2*k-1) * config(2*k-1), hz(2*k) * config(2*k), &
      !    & 2*k-1, config(2*k-1), 2*k, config(2*k)
      !  print *, "Swap contribution h"
      !  print "( F20.15,4X, F20.15,4X, *(2I0,2X)  )", hz(2*k-1) * config(2*k), hz(2*k) * config(2*k-1), &
      !    & 2*k, config(2*k), 2*k-1, config(2*k-1)
      !  print *, ""

      !print "( F15.10,4X,*(I0)  )", E(i), (1 - config(:))/2
      !print *, E(i), Es(i), exact_energy(nspin, Vzz, hz, l), exact_energy(nspin, Vzz, hz, m)

      !QE_alt(i) = (E(i) + Es(i)) / 2

      !do k = 1, nspin/2 - 1
      !  QE(i) = QE(i) + ( ( Vzz(2*k-1) * 2 * config(2*k-1) * config(2*k) + &
      !    & Vzz(2*k) * ( config(2*k-1) * config(2*k+2) + config(2*k) * config(2*k+1) )) + &
      !    & ( hz(2*k-1) + hz(2*k) ) * (config(2*k) + config(2*k-1)) )/2
      !enddo
      !k = nspin/2
      !QE(i) = QE(i) + ( ( Vzz(2*k-1) * 2 * config(2*k-1) * config(2*k) ) + &
      !  & ( hz(2*k-1) + hz(2*k) ) * (config(2*k) + config(2*k-1)) )/2

      !print "( A,4X, A,4X, A  )", "(E + Esigma)/2", "QE", "config"
      !print "( F15.10,4X, F15.10,4X, *(I0)  )", QE_alt(i), QE(i), (1 - config(:))/2 
    enddo

  end subroutine

  function exact_energy_LR(nspin, Vzz, hz, i) result(E)

    integer (c_int), intent(in) :: nspin, i
    real (c_double), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)

    integer (c_int) :: k, q, config(nspin), spin(nspin)
    real (c_double) :: E

    call decode(i, nspin, config)
    spin = 1 - 2*config
    E = 0
    do k = 1, nspin - 1
      E = E + hz(k) * spin(k)
      do q = k+1, nspin
        E = E + Vzz(k,q) * spin(k) * spin(q)
      enddo
    enddo
    k = nspin
    E = E + hz(k) * spin(k)

  end function

  function exact_energies_Sz0_LR(nspin, dim_Sz0, Vzz, hz) result(E)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)

    real (c_double) :: E(dim_Sz0)
    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0)

    call zero_mag_states(nspin, dim_Sz0, states)
    do i = 1, dim_Sz0

      l = states(i)
      call decode(l, nspin, config)
      E(i) = exact_energy_LR(nspin, Vzz, hz, l)

    enddo

  end function

  function exact_quasi_energy_LR(nspin, Vzz, hz, i) result(QE)

    integer (c_int), intent(in) :: nspin, i
    real (c_double), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)

    real (c_double) :: QE
    integer (c_int) :: k, j, config(nspin), config_swap(nspin)

    call decode(i, nspin, config)
    j = 0
    do k = 1, nspin/2
      j = j + config(2*k) * 2**(2*k-2)  + config(2*k-1) * 2**(2*k-1)
    enddo
    !call decode(j, nspin, config_swap)
    !print *, "Check swap: "
    !print "(A,2X,*(I0))", "config: ", config(:)
    !print "(A,2X,*(I0))", "config_swap: ", config_swap(:)

    QE = (exact_energy_LR(nspin, Vzz, hz, i) + exact_energy_LR(nspin, Vzz, hz, j)) / 2
    !print "(*(A26))", "E(i)", "E(j)", "QE = (E(i)+E(j))/2"
    !print *, exact_energy_LR(nspin, Vzz, hz, i), exact_energy_LR(nspin, Vzz, hz, j), QE

  end function


  function exact_quasi_energies_Sz0_LR(nspin, dim_Sz0, Vzz, hz) result(QE) !E, Es, QE, QE_alt)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    real (c_double), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)
    real (c_double) :: QE(dim_Sz0)

    integer (c_int) :: i, l, config(nspin), states(dim_Sz0)
    !real (c_double) :: mu

    QE = 0
    call zero_mag_states(nspin, dim_Sz0, states)
    !print *, "(Start) Exact Quasienergies"
    do l = 1, dim_Sz0

      !print *, "Quasienergy at l = ", l
      i = states(l)
      call decode(i, nspin, config)
      !mu = exact_quasi_energy_LR(nspin, Vzz, hz, i)
      QE(l) = exact_quasi_energy_LR(nspin, Vzz, hz, i)
      QE(l) = real(C_UNIT*log(exp(-C_UNIT*QE(l))), kind(dp))

      !print *, "Energies:    E_exact,   mu = E_exact + (E_pair) "
      !print *, exact_energy_LR(nspin, Vzz, hz, i), mu
      !print *, "First BZ (mu): "
      !print *, real(C_UNIT*log(exp(-C_UNIT*mu)), kind(dp)), mod(mu + pi,2*pi) - pi

      !print "(*(A16))", "i", "l", "config", "QE"
      !write (*,"(2(4X,I12))",advance='no') i, l
      !write (*,"(4X,*(I0))",advance='no'), config(:)
      !write (*,*) QE(l)

      !print *, ""
      !print *, ""

    enddo

    end function

  subroutine time_avg(option, steps, start, avg, sigma, t_avg, t_sigma)
    character, intent(in) :: option*1
    integer (c_int), intent(in) :: steps, start
    real (c_double), intent(in) :: avg(steps), sigma(steps)
    real (c_double), intent(out) :: t_avg, t_sigma
    integer (c_int) :: i, j, k
    
    t_avg = 0
    t_sigma = 0
    if (option == 'F') then
      do j = start, steps
        t_avg = t_avg + avg(j) !Time average of avg(t)
        t_sigma = t_sigma + sigma(j)**2
      enddo
    else if (option == 'T') then
      do j = start, steps
        t_avg = t_avg + 2*(mod(j,2)-0.5) * avg(j) !Time average of (-1)^t avg(t)
        t_sigma = t_sigma + sigma(j)**2
      enddo
    endif
    t_avg = t_avg/real(steps-start+1,c_double)
    t_sigma = sqrt(t_sigma/real(steps-start+1,c_double))

  end subroutine time_avg

end module observables
