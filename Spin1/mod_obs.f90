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
  
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  use printing
  use functions, search => binsearch_closest_in_circle
  implicit none

  complex (dcp), private, parameter :: C_ZERO = dcmplx(0._dp, 0._dp)
  complex (dcp), private, parameter :: C_ONE = dcmplx(1._dp, 0._dp) 
  complex (dcp), private, parameter :: C_UNIT = dcmplx(0._dp, 1._dp)

  integer (ip), private, parameter :: dimSpin1 = 3
  real (dp), parameter, private :: pi = 4._dp * atan(1._dp)

contains

  function sigmaz_Sz(nspin, dim_Sz, Sz, psi_Sz) result(sigmaz)
    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    complex (dcp), intent(in) :: psi_Sz(dim_Sz)
    real (dp) :: sigmaz(nspin)

    integer (ip) :: i, k, l, config(nspin), states(dim_Sz), s

    call basis_Sz(nspin, dim_Sz, Sz, states)
    sigmaz = 0
    do l = 1, dim_Sz
      i = states(l)
      call decode(i, nspin, config)
      do k = 1, nspin
        s = 1 - config(k)
        sigmaz(k) = sigmaz(k) + abs(psi_Sz(l))**2 * s
      enddo
    enddo

  end function


  function sigmaz_Sz0(nspin, dim_Sz0, psi_Sz0)
    integer (ip), intent(in) :: nspin, dim_Sz0
    complex (dcp), intent(in) :: psi_Sz0(dim_Sz0)
    real (dp) :: sigmaz_Sz0(nspin)

    integer (ip) :: i, k, l, config(nspin), states(dim_Sz0), s
    real (dp) :: sigmaz(nspin)

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

    integer(ip), intent(in) :: nspin, dim_Sz0
    complex(dcp), intent(in) :: state(dim_Sz0)
    real (dp) :: imbalance_Sz0
    real(dp) :: imb, imbaux
    integer (ip) :: i, k, l, indx(dim_Sz0), config(nspin), s

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
    integer (ip), intent(in) :: nspin, dim_Sz0
    complex (dcp), intent(in) :: psi_Sz0(dim_Sz0)
    real (dp) :: local_imbalance_Sz0

    integer (ip) :: i, k, l, config(nspin), states(dim_Sz0), s1, s2
    real (dp) :: LI, LI_part

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

    integer (ip), intent(in) :: nspin, dim_Sz0, q, p
    complex (dcp), intent(in) :: psi_Sz0(dim_Sz0)
    real (dp) :: sigmaz_corr_c_Sz0

    integer (ip) :: i, k, l, config(nspin), states(dim_Sz0), sp, sq
    real (dp) :: corr, avgq, avgp

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
!    integer(ip), intent(in) :: i, nspin
!    integer(ip) :: config(nspin)
!    real(dp) :: mag
!    integer (ip) :: j, k, m
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
!    integer(ip), intent(in) :: nspin, dim
!    complex(dcp), intent(in) :: state(dim)
!    real(dp) :: mag, magaux
!    integer :: config(nspin)
!    integer (ip) :: i, j, k, m
!
!    mag = 0
!    do i = 1, dim
!      call decode(i-1,nspin,config)
!      magaux = 0
!      do k = 1, nspin
!        magaux = magaux + (1._dp - 2._dp * config(k))
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
!    integer(ip), intent(in) :: nspin, dim, p
!    complex(dcp), intent(in) :: state(dim)
!    real(dp) :: mag
!    integer :: config(nspin)
!    integer (ip) :: i, j, k, m
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
!    integer(ip), intent(in) :: nspin, dim
!    complex(dcp), intent(in) :: state(dim)
!    real(dp) :: mag
!    real(dp) :: magaux
!    integer :: config(nspin)
!    integer (ip) :: i, j, k, m
!
!    mag = 0
!    do i = 1, dim
!      call decode(i-1,nspin,config)
!      magaux = 0
!      do k = 1, nspin
!        magaux = magaux + (-1)**k * (1._dp - 2._dp * config(k))
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
!    integer(ip), intent(in) :: nspin, dim
!    complex(dcp), intent(in) :: state(dim)
!    real(dp) :: imb, imbaux, mag, magaux
!    integer :: config(nspin)
!    integer (ip) :: i, j, k, m
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
!    integer(ip), intent(in) :: nspin, dim_Sz0
!    complex(dcp), intent(in) :: state(dim_Sz0)
!    real(dp) :: imb, imbaux
!    integer :: config(nspin)
!    integer (ip) :: i, j, k, m, l, indx(dim_Sz0)
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
!    integer(ip), intent(in) :: nspin, dim_Sz0
!    complex(dcp), intent(in) :: state(dim_Sz0)
!    real(dp) :: imb, imbaux
!    integer :: config(nspin)
!    integer (ip) :: i, k1, k2, l, indx(dim_Sz0)
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
!    integer(ip), intent(in) :: nspin, i
!    real (dp) :: imb, imbalance_basis
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
!    integer(ip), intent(in) :: nspin, i
!    real (dp) :: imb, imbalance_sq_basis
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
!    integer (ip), intent(in) :: nspin, i
!    real (dp) :: local_imbalance_basis
!
!    integer (ip) :: k, config(nspin)
!    real (dp) :: LI
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
!    integer (ip), intent(in) :: nspin, dim
!    complex (dcp), intent(in) :: psi(dim)
!    real (dp) :: local_imbalance
!
!    integer (ip) :: i, k, config(nspin)
!    real (dp) :: LI, LI_part
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
!    integer (ip), intent(in) :: nspin, dim
!    complex (dcp), intent(in) :: psi(dim)
!    real (dp) :: local_overlap
!
!    integer (ip) :: i, j, k, config(nspin)
!    real (dp) :: LO
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
!    integer (ip), intent(in) :: nspin, dim_Sz0
!    complex (dcp), intent(in) :: psi(dim_Sz0)
!    real (dp) :: local_overlap_Sz0
!
!    integer (ip) :: i, j, k, r, l, config(nspin), states(dim_Sz0)
!    real (dp) :: LO
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
!    integer (ip), intent(in) :: nspin, dim, q, p
!    complex (dcp), intent(in) :: psi(dim)
!    real (dp) :: sigmaz_corr_c
!
!    integer (ip) :: i, k, config(nspin)
!    real (dp) :: corr, avgq, avgp
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
  function exact_energy(nspin, Vz, hz, sigmaz) result(E)

    integer (ip), intent(in) :: nspin
    real (dp), intent(in) :: Vz(nspin-1), hz(nspin)
    integer (ip), intent(in) :: sigmaz(nspin)
    real (dp) :: E

    integer (ip) :: k

    E = 0
    do k = 1, nspin - 1
      E = E + Vz(k) * sigmaz(k) * sigmaz(k+1) + hz(k) * sigmaz(k)
    enddo
    k = nspin
    E = E + hz(k) * sigmaz(k)

  end function


  function exact_energies_Sz(nspin, dim_Sz, Sz, Vz, hz) result(E)

    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    real (dp), intent(in) :: Vz(nspin-1), hz(nspin)
    real (dp) :: E(dim_Sz)

    integer (ip) :: i, k, l, config(nspin), idxSz(dim_Sz), sigmaz(nspin)

    call basis_Sz(nspin, dim_Sz, Sz, idxSz)
    do l = 1, dim_Sz

      i = idxSz(l)
      call decode(i, nspin, config)
      sigmaz = 1 - config
      E(l) = exact_energy(nspin, Vz, hz, sigmaz)

    enddo

  end function 

  function exact_quasi_energy(nspin, Vz, hz, sigmaz) result(QE)

    integer (ip), intent(in) :: nspin
    real (dp), intent(in) :: Vz(nspin-1), hz(nspin)
    integer (ip), intent(in) :: sigmaz(nspin)
    real (dp) :: QE

    integer (ip) :: k, sigmaz_swap(nspin)

    do k = 1, nspin/2
      sigmaz_swap(2*k-1) = sigmaz(2*k)
      sigmaz_swap(2*k) = sigmaz(2*k-1)
    enddo
    QE = ( exact_energy(nspin, Vz, hz, sigmaz) + exact_energy(nspin, Vz, hz, sigmaz_swap) ) / 2
    QE = real( C_UNIT * log( exp(-C_UNIT*QE ) ), dp )

  end function

  function exact_quasi_energies_Sz(nspin, dim_Sz, Sz, Vz, hz) result(QE)

    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    real (dp), intent(in) :: Vz(nspin-1), hz(nspin)
    real (dp) :: QE(dim_Sz)

    integer (ip) :: i, k, l, config(nspin), idxSz(dim_Sz), sigmaz(nspin)

    call basis_Sz(nspin, dim_Sz, Sz, idxSz)
    do l = 1, dim_Sz

      i = idxSz(l)
      call decode(i, nspin, config)
      sigmaz = 1 - config
      QE(l) = exact_quasi_energy(nspin, Vz, hz, sigmaz)

    enddo

  end function 
!
!
!  subroutine exact_quasi_energies_Sz0(nspin, dim_Sz0, V_int, h_z, QE) !E, Es, QE, QE_alt)
!
!    integer (ip), intent(in) :: nspin, dim_Sz0
!    real (dp), intent(in) :: V_int(nspin-1), h_z(nspin)
!    real (dp), intent(out) :: QE(dim_Sz0)
!
!    integer (ip) :: i, k, l, m, config(nspin), states(dim_Sz0)
!    !real (dp) :: E(dim_Sz0), Es(dim_Sz0), QE_alt(dim_Sz0)
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
    integer (ip), intent(in) :: steps, start
    real (dp), intent(in) :: avg(steps), sigma(steps)
    real (dp), intent(out) :: t_avg, t_sigma
    integer (ip) :: i, j, k
    
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
    t_avg = t_avg/real(steps-start+1,dp)
    t_sigma = sqrt(t_sigma/real(steps-start+1,dp))

  end subroutine time_avg


  !------- Spectral Properties -------!

  subroutine gap_ratio(energies, r_avg, r_sq)

    real (dp), intent(in) :: energies(:)
    real (dp), intent(out) :: r_avg, r_sq

    integer :: n, dim
    real (dp) :: gap1, gap2

    dim = size(energies)
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

  function IPR(psi)
    complex (dcp), intent(in) :: psi(:)
    real (dp) :: IPR
    integer (ip) :: dim, i

    dim = size(psi)
    !print *, dim
    IPR = 0
    do i = 1, dim
      IPR = IPR + abs(psi(i))**4
      !print*, IPR, abs(psi(i))
    enddo

  end function IPR

  function pi_pair(QE) result(beta)

    real (dp), intent(in) :: QE(:)

    integer (ip) :: alpha, dim
    integer (ip), allocatable :: beta(:)
    real (dp) :: val

    dim = size(QE)
    allocate(beta(dim))

    do alpha = 1, dim
      val = mod(QE(alpha) + 2*pi, 2*pi) - pi

      !print *, "alpha = ", alpha, "QE(alpha) = ", QE(alpha), "val =", val
      beta(alpha) = search(val, QE) !<--------- If you want the "minimal pi-distance from above"

    enddo

  end function

  subroutine log_gap_difference(QE, log_pair_avg, log_pair_sq, log_near_avg, log_near_sq, log_avg, log_sq)

    !Computes the averages of the logarithm of gaps of neighbouring eigenvalues 
    !and paired eigenvalues (where pi-angle-distance has been minimized)
    ! log_near_avg = < log( E_{alpha+1} - E_alpha ) >
    ! log_pair_avg = < log( E_beta - E_alpha ) >
    ! log_avg = log_pair_avg - log_near_avg
    ! log_sq = < (log Delta^alpha - log Delta_0^alpha)^2 >
    ! If there is spectral pairing: Delta^alpha / Delta_0^alpha      -> 0  as L-> +infty
    !                               log(Delta^alpha / Delta_0^alpha) -> -infty


    real (dp), intent(in) :: QE(:)
    real (dp), intent(out) :: log_pair_avg, log_pair_sq, log_near_avg, log_near_sq, log_avg, log_sq

    real (dp), allocatable :: pair(:), near(:)
    integer (ip), allocatable :: pi_paired(:)

    integer (ip) :: dim, alpha, beta !, beta1, beta2
    real (dp) :: val

    dim = size(QE)
    allocate(pair(dim), near(dim), pi_paired(dim))

    pi_paired = pi_pair(QE)

    !print "(2(A12),3(A26))", "alpha", "beta", "E(alpha)", "(E(alpha)+pi)_1", "E(beta)"
    !do alpha = 1, dim
    !  beta = pi_paired(alpha)
    !  print *, alpha, beta, QE(alpha), mod(QE(alpha) + 2*pi, 2*pi) - pi, QE(beta)
    !enddo

    !beta1 = 0
    !print "(2(A12),3(A26))", "alpha", "beta", "dQE(beta-1)", "dQE(beta)", "dQE(beta+1)"
    !do alpha = 1, dim
    !  beta = pi_paired(alpha)
    !  val = mod(QE(alpha) + 2*pi, 2*pi) - pi
    !  !print *, alpha, beta, abs(QE(beta1)-val), abs(QE(beta) - val), abs(QE(beta2)-val)
    !  print *, alpha, beta, minloc(abs(QE-val)), minloc( abs(exp(C_UNIT*(QE-val)) - 1), dim )
    !  beta1 = beta1 + (beta - minloc( abs(exp(C_UNIT*(QE-val)) - 1), dim ) )
    !enddo
    !print *, "Total difference between binsearch and minloc|e^(i(QE-val))-1|:", beta1


    !print "(*(A26))", "Delta_pi", "Delta_0", "log(Delta_pi)", "log(Delta_0)"
    do alpha = 1, dim 

      beta = merge(alpha+1,1,alpha.ne.dim)
      near(alpha) = QE(beta) + merge(0._dp,2*pi,alpha.ne.dim) - QE(alpha)

      beta = pi_paired(alpha)
      pair(alpha) = abs(abs(QE(beta) - QE(alpha)) - pi)

      !print *, pair(alpha), near(alpha), log(pair(alpha)), log(near(alpha))

    enddo
    pair = max(pair, epsilon(pair))
    near = max(near, epsilon(near))

    log_pair_avg = sum(log(pair)) / dim
    log_pair_sq = sum((log(pair))**2) / dim

    log_near_avg = sum(log(near)) / dim
    log_near_sq = sum((log(near))**2) / dim

    log_avg = log_pair_avg - log_near_avg
    log_sq = sum((log(pair) - log(near))**2) / dim

  end subroutine

  subroutine log_gap_difference_half_spectrum_shift(QE, log_pair_avg, log_pair_sq, log_near_avg, log_near_sq, log_avg, log_sq)

    !Computes the averages of the logarithm of gaps of neighbouring eigenvalues 
    !and paired eigenvalues (separated by half the spectrum)
    ! log_near_avg = < log( E_{alpha+1} - E_alpha ) >
    ! log_pair_avg = < log( E_{alpha+N/2} - E_alpha ) >
    ! log_avg = log_pair_avg - log_near_avg
    ! log_sq = < (log Delta^alpha - log Delta_0^alpha)^2 >
    ! If there is spectral pairing: Delta^alpha / Delta_0^alpha      -> 0  as L-> +infty
    !                               log(Delta^alpha / Delta_0^alpha) -> -infty


    real (dp), intent(in) :: QE(:)
    real (dp), intent(out) :: log_pair_avg, log_pair_sq, log_near_avg, log_near_sq, log_avg, log_sq

    real (dp), allocatable :: pair(:), near(:)
    integer (ip), allocatable :: pi_paired(:)

    integer (ip) :: dim, alpha, beta !, beta1, beta2
    real (dp) :: val

    dim = size(QE)
    allocate(pair(dim), near(dim), pi_paired(dim))


    !print "(*(A26))", "(Delta_pi)_half", "Delta_0", "log(Delta_pi)_half", "log(Delta_0)"
    do alpha = 1, dim 

      beta = merge(alpha+1,1,alpha.ne.dim)
      near(alpha) = QE(beta) + merge(0._dp,2*pi,alpha.ne.dim) - QE(alpha)

      beta = merge(alpha+dim/2, alpha-dim/2, alpha<=dim/2)
      pair(alpha) = abs(abs(QE(beta) - QE(alpha)) - pi)
      !print *, pair(alpha), near(alpha), log(pair(alpha)), log(near(alpha))

    enddo
    pair = max(pair, epsilon(pair))
    near = max(near, epsilon(near))

    log_pair_avg = sum(log(pair)) / dim
    log_pair_sq = sum((log(pair))**2) / dim

    log_near_avg = sum(log(near)) / dim
    log_near_sq = sum((log(near))**2) / dim

    log_avg = log_pair_avg - log_near_avg
    log_sq = sum((log(pair) - log(near))**2) / dim

  end subroutine

  function exact_energy_LR(nspin, Vzz, hz, i) result(E)

    integer (ip), intent(in) :: nspin, i
    real (dp), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)

    integer (ip) :: k, q, config(nspin), spin(nspin)
    real (dp) :: E

    call decode(i, nspin, config)
    spin = 1 - config
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

  function exact_quasi_energy_LR(nspin, Vzz, hz, i) result(QE)

    integer (ip), intent(in) :: nspin, i
    real (dp), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)

    real (dp) :: QE
    integer (ip) :: k, j, config(nspin)!, config_swap(nspin)

    call decode(i, nspin, config)
    j = 0
    do k = 1, nspin/2
      j = j + config(2*k) * dimSpin1**(2*k-2)  + config(2*k-1) * dimSpin1**(2*k-1)
    enddo
    !call decode(j, nspin, config_swap)
    !print "(A15, *(I0))", "config = ", config(:)
    !print "(A15, *(I0))", "config_swap = ", config_swap(:)

    QE = (exact_energy_LR(nspin, Vzz, hz, i) + exact_energy_LR(nspin, Vzz, hz, j)) / 2 + &
      & merge(pi,0.0_dp,j>i)
    QE = real(C_UNIT*log(exp(-C_UNIT*QE)), kind=dp)

  end function

  function exact_quasi_energies_LR_Sz(nspin, dim_Sz, Sz, Vzz, hz) result(QE)

    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    real (dp), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)
    real (dp) :: QE(dim_Sz)

    integer (ip) :: i, k, l, config(nspin), idxSz(dim_Sz)
    !real (dp) :: mu

    QE = 0
    call basis_Sz(nspin, dim_Sz, Sz, idxSz)
    !print *, "(Start) Exact Quasienergies"
    do l = 1, dim_Sz

      !print *, "Quasienergy at l = ", l
      i = idxSz(l)
      call decode(i, nspin, config)
      !mu = exact_quasi_energy_LR(nspin, Vzz, hz, i)
      QE(l) = exact_quasi_energy_LR(nspin, Vzz, hz, i)

      !print *, "Energies:    E_exact,   mu = E_exact + (E_pair) "
      !print *, exact_energy_LR(nspin, Vzz, hz, i), mu
      !print *, "First BZ (mu): "
      !print *, real(C_UNIT*log(exp(-C_UNIT*mu)), kind=dp), mod(mu + pi,2*pi) - pi

      !print "(*(A16))", "i", "l", "config", "QE"
      !write (*,"(2(4X,I12))",advance='no') i, l
      !write (*,"(4X,*(I0))",advance='no'), config(:)
      !write (*,*) QE(l)

      !print *, ""
      !print *, ""

    enddo

  end function

  function exact_energies_LR_Sz(nspin, dim_Sz, Sz, Vzz, hz) result(E)

    integer (ip), intent(in) :: nspin, dim_Sz, Sz
    real (dp), intent(in) :: Vzz(nspin-1,nspin), hz(nspin)
    real (dp) :: E(dim_Sz)

    integer (ip) :: i, k, l, config(nspin), idxSz(dim_Sz)

    E = 0
    call basis_Sz(nspin, dim_Sz, Sz, idxSz)
    do l = 1, dim_Sz

      i = idxSz(l)
      E(l) = exact_energy_LR(nspin, Vzz, hz, i)

    enddo

  end function

end module
