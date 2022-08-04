module MBL_subrtns

  use iso_c_binding
  use genmat
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


  subroutine reduced_DM(nspin, nspin_A, dim_A, psi, rho_A)

    !Computes the reduced density matrix for a pure state for the part of length nspin_A
    ! rho_A = Tr_B(rho) = Tr_B( |psi> <psi| )
    !Traces out the remaining (nspin - nspin_A) spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, dim_A
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_A(dim_A,dim_A)

    integer :: i, j, k, iA, jA, iB, jB, kB, dim_B

    dim_B = dim/dim_A
    rho_A = 0

    do jA = 1, dim_A
      do kA = 1, dim_A
        do iB = 1, dim_B

          rho_A(jA,kA) = rho_A(jA,kA) + psi(jA + 2**(nspin_A) * iB) * dconjg( psi(kA + 2**(nspin_A) * iB )

        enddo
      enddo
    enddo

  end subroutine reduced_DM





  subroutine entanglement(dim, rho, EE)

    !Computes the Entanglement Entropy for a given density matrix rho of dimension dim
    ! EE = S = -Tr(rho * ln(rho))

    integer, intent(in) :: dim
    complex (c_double_complex), intent(in) :: rho(dim,dim)
    real (c_double), intent(out) :: EE(dim)


  end subroutine entanglement




  
end module MBL_subrtns


module MBL

  use MBL_subrtns

  interface gap_ratio
    module procedure gap_ratio_R!, level_statistic_C
  end interface

end module MBL



