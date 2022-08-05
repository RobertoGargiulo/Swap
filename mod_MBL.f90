module MBL_subrtns

  use iso_c_binding
  use genmat
  use exponentiate
  implicit none

contains

  subroutine gap_ratio_R(dim, energies, r_avg, r_sigma)

    integer, intent(in) :: dim
    real (c_double), intent(in) :: energies(dim)

    integer :: n
    real (c_double) :: r_avg, r_sigma, gap1, gap2

    r_avg = 0
    r_sigma = 0
    do n = 1, dim - 2

      gap1 = energies(n+1) - energies(n)
      gap2 = energies(n+2) - energies(n+1)

      r_avg = r_avg + min(gap1,gap2)/max(gap1,gap2)
      r_sigma = r_sigma + (min(gap1,gap2)/max(gap1,gap2))**2

    enddo

    r_avg = r_avg/real(dim-2)
    r_sigma = sqrt( real(dim-2)/real(dim-3) * ( r_sigma/real(dim-2) - r_avg**2  ) ) / sqrt(real(dim-2))

  end subroutine gap_ratio_R

  !real function avg_gap_ratio_C(nspin, quasi_energies)

    !

  !end function avg_gap_ratio_C


  subroutine left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)

    !Computes the reduced density matrix for a pure state of the first (from the left) 'nspinA' spins in the chain
    ! rho_A = Tr_B(rho) = Tr_B( |psi> <psi| )
    !Traces out the remaining (nspin - nspin_A) spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, nspin_A, dim, dim_A
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_A(dim_A,dim_A)

    integer :: jA, kA, iB, dim_B

    dim_B = dim/dim_A
    rho_A = 0

    do jA = 1, dim_A
      do kA = 1, dim_A
        do iB = 1, dim_B

          rho_A(jA,kA) = rho_A(jA,kA) + psi(jA + 2**(nspin_A) * iB) * dconjg( psi(kA + 2**(nspin_A) * iB ) )

        enddo
      enddo
    enddo

  end subroutine left_reduced_DM

  subroutine right_reduced_DM(nspin, nspin_B, dim, dim_B, psi, rho_B)

    !Computes the reduced density matrix for a pure state of the first (from the left) 'nspinA' spins in the chain
    ! rho_A = Tr_B(rho) = Tr_B( |psi> <psi| )
    !Traces out the remaining (nspin - nspin_A) spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, nspin_B, dim, dim_B
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_B(dim_B,dim_B)

    integer :: jB, kB, iA, dim_A

    dim_A = dim/dim_B
    rho_B = 0

    do jB = 1, dim_B
      do kB = 1, dim_B
        do iA = 1, dim_A

          rho_B(jB,kB) = rho_B(jB,kB) + psi(iA + 2**(nspin-nspin_B) * jB) * dconjg( psi(iA + 2**(nspin-nspin_B) * kB ) )

        enddo
      enddo
    enddo

  end subroutine right_reduced_DM




  real function entanglement(dim, rho)

    !Computes the Entanglement (von Neumann) Entropy for a given density matrix rho of dimension dim
    ! EE = S = -Tr(rho * ln(rho))
    ! For a reduced density matrix S_A = -Tr_A(rho_A ln(rho_A))

    integer, intent(in) :: dim
    real (c_double_complex), intent(in) :: rho(dim,dim)

    real (c_double) :: EE, prob(dim), W(dim,dim)
    integer :: i

    call diagSYM( 'V', dim, rho, prob, W)

    EE = 0
    do i = 1, dim
      EE = EE - prob(i) * log(prob(i))
    enddo
  
    entanglement = EE

  end function entanglement




  
end module MBL_subrtns


module MBL

  use MBL_subrtns

  interface gap_ratio
    module procedure gap_ratio_R!, level_statistic_C
  end interface

end module MBL



