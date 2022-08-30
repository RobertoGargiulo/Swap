module MBL_subrtns

  use iso_c_binding
  use genmat
  use exponentiate
  use printing
  implicit none

  real (c_double), parameter, private :: pi = 4.d0 * datan(1.d0)

contains

  subroutine gap_ratio_R(dim, energies, r_avg, r_sq)

    integer, intent(in) :: dim
    real (c_double), intent(in) :: energies(dim)

    integer :: n
    real (c_double) :: r_avg, r_sq, gap1, gap2

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

  end subroutine gap_ratio_R

  !real function avg_gap_ratio_C(nspin, quasi_energies)

    !

  !end function avg_gap_ratio_C

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

    integer (c_int) :: i, j

    do i = 1, dim 
      
      if (i<=dim-1) then
        j = i + 1
        near(i) = energies(j) - energies(i)
      else if (i>dim-1) then
        j = i + 1 -dim
        near(i) = energies(j) + 2*pi - energies(i)
      endif
  
      if (i<=dim/2) then
        j = i + dim/2
        pair(i) = abs( energies(j) - energies(i) - pi )
      else if (i>dim/2) then
        j = i - dim/2
        pair(i) = abs( energies(j) + 2*pi - energies(i) - pi  )
      endif

    enddo

    log_pair_avg = sum(log(pair)) / dim
    log_near_avg = sum(log(near)) / dim
    log_avg = log_pair_avg - log_near_avg
    log_sq = sum((log(pair) - log(near))**2) / dim

  end subroutine log_gap_difference



  subroutine left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)

    !Computes the reduced density matrix for a pure state of the first (from the left) 'nspinA' spins in the chain
    ! rho_A = Tr_B(rho) = Tr_B( |psi> <psi| )
    !Traces out the remaining (nspin - nspin_A) spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, nspin_A, dim, dim_A
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_A(dim_A,dim_A)

    integer :: jA, kA, iB, dim_B, i, j

    i = nspin
    dim_B = dim/dim_A
    rho_A = 0
    do jA = 0, dim_A-1
      do kA = 0, dim_A-1
        do iB = 0, dim_B-1

          i = jA + 2**(nspin_A) * iB + 1
          j = kA + 2**(nspin_A) * iB + 1

          rho_A(jA+1,kA+1) = rho_A(jA+1,kA+1) + psi(i) * dconjg( psi(j) )

        enddo
      enddo
    enddo

  end subroutine left_reduced_DM

  subroutine right_reduced_DM(nspin, nspin_B, dim, dim_B, psi, rho_B)

    !Computes the reduced density matrix for a pure state of the last (from the left) 'nspinA' spins in the chain
    ! rho_B = Tr_A(rho) = Tr_A( |psi> <psi| )
    !Traces out the remaining (nspin - nspin_B) spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, nspin_B, dim, dim_B
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_B(dim_B,dim_B)

    integer :: jB, kB, iA, dim_A, i, j

    dim_A = dim/dim_B
    rho_B = 0

    do jB = 0, dim_B-1
      do kB = 0, dim_B-1
        do iA = 0, dim_A-1

          i = iA + 2**(nspin-nspin_B) * jB + 1
          j = iA + 2**(nspin-nspin_B) * kB + 1

          rho_B(jB+1,kB+1) = rho_B(jB+1,kB+1) + psi(i) * dconjg( psi(j) )

        enddo
      enddo
    enddo

  end subroutine right_reduced_DM

  subroutine edges_reduced_DM(nspin, nspin_A, nspin_B, dim, dim_A, dim_B, psi, rho_AB)

    !Computes the reduced density matrix for a pure state of the first (from the left) 'nspinA' spins and last 'nspinB' spins in the chain
    ! rho_AB = Tr_C(rho) where C is the remaining region
    !Traces out the remaining (nspin - nspin_A - nspin_B) spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, nspin_A, nspin_B, dim, dim_A, dim_B
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_AB(dim_A*dim_B,dim_A*dim_B)

    integer :: jAB, kAB, jA, jB, kA, kB, iC, dim_C, i, j

    dim_C = dim/(dim_A*dim_B)
    rho_AB = 0

    do jA = 0, dim_A-1
      do jB = 0, dim_B-1
        do kA = 0, dim_A-1
          do kB = 0, dim_B-1
            do iC = 0, dim_C-1

              jAB = jA + 2**(nspin_A)*jB + 1
              kAB = kA + 2**(nspin_A)*kB + 1

              i = jA + 2**(nspin_A)*iC + 2**(nspin-nspin_B) * jB + 1
              j = kA + 2**(nspin_A)*iC + 2**(nspin-nspin_B) * kB + 1 

              rho_AB(jAB,kAB) = rho_AB(jAB,kAB) + psi(i) * dconjg( psi(j) )
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine edges_reduced_DM



  function entanglement(dim, rho)

    !Computes the Entanglement (von Neumann) Entropy for a given density matrix rho of dimension dim, using Exact Diagonalization
    ! EE = S(rho) = -Tr(rho * ln(rho))
    ! For a reduced density matrix S_A = -Tr_A(rho_A ln(rho_A))

    integer, intent(in) :: dim
    complex (c_double_complex), intent(in) :: rho(dim,dim)
    real (c_double) :: entanglement

    real (c_double) :: EE, prob(dim)
    complex (c_double_complex), allocatable :: W(:,:)
    integer :: i

    allocate(W(dim,dim))
    call diagHE( 'N', dim, rho, prob, W)
    deallocate(W)

    EE = 0
    do i = 1, dim
      if(prob(i) > 0 ) then
        EE = EE - prob(i) * log(prob(i))
        !if (prob(i) > 1.0e-3) print *, prob(i), i
      endif
    enddo
  
    entanglement = EE

  end function entanglement


  function mutual_information(nspin, nspin_A, dim, dim_A, psi)

    !Computes the long-range mutual information between the two edges.
    !The edges are regions (spin chains) A,B both of length nspin_A and local dimension dim_A (here case of Full Hilbert Space)
    !psi is the (pure) state of the spin chain
    ! MI = S_A + S_B - S_AB
    ! MI is the maximum information of A we can find from B, and viceversa

    integer, intent(in) :: nspin_A, nspin, dim_A, dim
    complex (c_double_complex), intent(in) :: psi(dim)
    real (c_double) :: mutual_information

    real (c_double) :: MI
    complex (c_double_complex) :: rho_A(dim_A,dim_A), rho_B(dim_A,dim_A), rho_AB(dim_A*dim_A,dim_A*dim_A)

    call left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)
    call right_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_B)
    call edges_reduced_DM(nspin, nspin_A, nspin_A, dim, dim_A, dim_A, psi, rho_AB)

    MI = entanglement(dim_A,rho_A) + entanglement(dim_A,rho_B) - entanglement(dim_A*dim_A,rho_AB)
    !print *, nspin_A, entanglement(dim_A,rho_A), entanglement(dim_A,rho_B), entanglement(dim_A*dim_A,rho_AB), MI
    
    mutual_information = MI

  end function mutual_information

  function mutual_information_Sz0(nspin, nspin_A, dim, dim_Sz0, dim_A, psi_Sz0)

    !Computes the long-range mutual information between the two edges.
    !The edges are regions (spin chains) A,B both of length nspin_A and local dimension dim_A (here case of subspace Sz=0)
    !psi_Sz0 is the (pure) state of the spin chain
    ! MI = S_A + S_B - S_AB
    ! MI is the maximum information of A we can find from B, and viceversa

    integer, intent(in) :: nspin_A, nspin, dim_A, dim, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: mutual_information_Sz0

    real (c_double) :: MI, MBEE
    complex (c_double_complex) :: rho_A(dim_A,dim_A), rho_B(dim_A,dim_A), rho_AB(dim_A*dim_A,dim_A*dim_A)!, rho(dim_Sz0,dim_Sz0)
    complex (c_double_complex) :: psi(dim)
    integer :: i, l, states(dim_Sz0), config(nspin)

    real (c_double) :: var1

    call buildState_Sz0_to_FullHS(nspin, dim, dim_Sz0, psi_Sz0, psi)

    call left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)
    call right_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_B)
    call edges_reduced_DM(nspin, nspin_A, nspin_A, dim, dim_A, dim_A, psi, rho_AB)
    
    !print *, "rho_A = "
    !call printmat(dim_A, rho_A, 'C')
    !print *, "rho_B = "
    !call printmat(dim_A, rho_B, 'C')
    !print *, "rho_AB = "
    !call printmat(dim_A*dim_A, rho_AB, 'C')

    MI = entanglement(dim_A,rho_A) + entanglement(dim_A,rho_B) - entanglement(dim_A*dim_A,rho_AB)

    !print *, "S_A         S_B          S_{AB}"
    !print *, nspin_A, entanglement(dim_A,rho_A), entanglement(dim_A,rho_B), entanglement(dim_A*dim_A,rho_AB), MI
    !print *, MI, imbalance_sq_Sz0(nspin, dim_Sz0, psi_Sz0), IPR(psi)
    
    mutual_information_Sz0 = MI

  end function mutual_information_Sz0


  function max_bipartite_entropy_Sz0(nspin, dim, dim_Sz0, psi_Sz0)

    !Computes the long-range mutual information between the two edges.
    !The edges are regions (spin chains) A,B both of length nspin_A and local dimension dim_A (here case of subspace Sz=0)
    !psi_Sz0 is the (pure) state of the spin chain
    ! MI = S_A + S_B - S_AB
    ! MI is the maximum information of A we can find from B, and viceversa

    integer, intent(in) :: nspin, dim, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: max_bipartite_entropy_Sz0

    real (c_double) :: MBEE, pMBEE, BEE(nspin-1), BEEvar
    integer (c_int) :: nspin_A, dim_A
    complex (c_double_complex), allocatable :: rho_A(:,:) !, rho(dim_Sz0,dim_Sz0)
    complex (c_double_complex) :: psi(dim)


    call buildState_Sz0_to_FullHS(nspin, dim, dim_Sz0, psi_Sz0, psi)
    MBEE = 0
    pMBEE = 0
    !print *, "    nspin_A        dim_A     BEE"
    do nspin_A = 1, nspin - 1
      dim_A = 2**nspin_A
      allocate(rho_A(dim_A,dim_A))
      if (nspin_A <= nspin/2) then
        !print *, "rho_A allocated"
        call left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)
      else 
        call right_reduced_DM(nspin, nspin-nspin_A, dim, dim/dim_A, psi, rho_A)
      endif
      BEEvar = entanglement(dim_A,rho_A)
      BEE(nspin_A) = BEEvar
      if (BEEvar > pMBEE) MBEE = BEEvar
      pMBEE = BEEvar
      !print *, "    nspin_A        dim_A     BEE"
      !print *, nspin_A, dim_A, BEE(nspin_A)
      !call printmat(dim_A, rho_A, 'A')
      deallocate(rho_A)
    enddo

    max_bipartite_entropy_Sz0 = maxval(BEE)
    !print *, max_bipartite_entropy_Sz0, MBEE


  end function max_bipartite_entropy_Sz0


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

  integer function int_2dto1d(int_2d)

    implicit none
    integer, intent(in) :: int_2d(2)
    integer :: int_1d
    

    if (int_2d(1) < int_2d(2)) then
      int_1d = int_2d(2)**2 + int_2d(1)
    else if (int_2d(1) >= int_2d(2)) then
      int_1d = int_2d(1)**2 + int_2d(1) + int_2d(2)
    endif

    int_2dto1d = int_1d

  end function int_2dto1d


  function int_1dto2d(int_1d)

    implicit none
    integer, intent(in) :: int_1d
    integer :: int_2d(2)
    integer :: int_1dto2d(2)

    if (int_1d - floor(sqrt(real(int_1d)))**2 < floor(sqrt(real(int_1d))) ) then
      int_2d(1) = int_1d - floor(sqrt(real(int_1d)))**2
      int_2d(2) = floor(sqrt(real(int_1d)))
    else
      int_2d(1) = floor(sqrt(real(int_1d)))
      int_2d(2) = int_1d - floor(sqrt(real(int_1d)))**2 - floor(sqrt(real(int_1d)))
    endif

    int_1dto2d = int_2d

  end function int_1dto2d

  
end module MBL_subrtns


module MBL

  use MBL_subrtns

  interface gap_ratio
    module procedure gap_ratio_R!, level_statistic_C
  end interface

end module MBL



