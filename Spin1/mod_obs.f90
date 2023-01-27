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
  
  use iso_c_binding
  use printing
  use functions
  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double) 
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  integer (c_int), private, parameter :: dimSpin1 = 3

contains

  function sigmaz_Sz0(nspin, dim_Sz0, psi_Sz0)
    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: sigmaz_Sz0(nspin)

    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0), s
    real (c_double) :: sigmaz(nspin)

    call zero_mag_states(nspin, dim_Sz0, states)
    sigmaz = 0
    do l = 1, dim_Sz0
      i = states(l)
      call decode(i, nspin, config)
      do k = 1, nspin
        s = 1 - config(k)
        sigmaz(k) = sigmaz(k) + abs(psi_Sz0(l))**2 * s
      enddo
    enddo
    sigmaz_Sz0 = sigmaz

  end function

  function imbalance_Sz0(nspin, dim_Sz0, state)

    integer(c_int), intent(in) :: nspin, dim_Sz0
    complex(c_double_complex), intent(in) :: state(dim_Sz0)
    real (c_double) :: imbalance_Sz0
    real(c_double) :: imb, imbaux
    integer (c_int) :: i, k, l, indx(dim_Sz0), config(nspin), s

    imb = 0
    call zero_mag_states(nspin, dim_Sz0, indx)

    do l = 1, dim_Sz0

      i = indx(l)
      call decode(i,nspin,config)
      imbaux = 0
      do k = 1, nspin
        s = 1 - config(k)
        imbaux = imbaux + (-1)**k * s
      enddo
      imbaux = imbaux * abs(state(l))**2
      imb = imb + imbaux
    enddo
    imb = imb/nspin
    imbalance_Sz0 = imb

  end function imbalance_Sz0

  function local_imbalance_Sz0(nspin, dim_Sz0, psi_Sz0)
    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: local_imbalance_Sz0

    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0), s1, s2
    real (c_double) :: LI, LI_part

    call zero_mag_states(nspin, dim_Sz0, states)
    LI = 0
    do l = 1, dim_Sz0
      i = states(l)
      call decode(i, nspin, config)
      LI_part = 0
      do k = 1, nspin/2
        s1 = 1 - config(2*k-1)
        s2 = 1 - config(2*k)
        LI_part = LI_part + (s2 - s1)**2
      enddo
      LI = LI + abs(psi_Sz0(l))**2 * LI_part / (nspin/2)
    enddo
    local_imbalance_Sz0 = LI

  end function

  function sigmaz_corr_c_Sz0(nspin, dim_Sz0, q, p, psi_Sz0)

    integer (c_int), intent(in) :: nspin, dim_Sz0, q, p
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: sigmaz_corr_c_Sz0

    integer (c_int) :: i, k, l, config(nspin), states(dim_Sz0), sp, sq
    real (c_double) :: corr, avgq, avgp

    avgq = 0
    avgp = 0
    corr = 0
    call zero_mag_states(nspin, dim_Sz0, states)
    do l = 1, dim_Sz0
      
      i = states(l)
      call decode(i,nspin,config)

      sq = 1 - config(q)
      sp = 1 - config(p)
      corr = corr + abs(psi_Sz0(l))**2 * sq * sp
      avgq = avgq + abs(psi_Sz0(l))**2 * sq
      avgp = avgp + abs(psi_Sz0(l))**2 * sp
        
    enddo

    sigmaz_corr_c_Sz0 = corr - avgq * avgp


  end function




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
!
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
        t_avg = t_avg + avg(j)
        t_sigma = t_sigma + sigma(j)**2
      enddo
    else if (option == 'T') then
      do j = start, steps
        t_avg = t_avg + 2*(mod(j,2)-0.5) * avg(j)
        t_sigma = t_sigma + sigma(j)**2
      enddo
    endif
    t_avg = t_avg/real(steps-start+1,c_double)
    t_sigma = sqrt(t_sigma/real(steps-start+1,c_double))

  end subroutine time_avg


end module
