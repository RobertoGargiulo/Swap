module entanglement

  use iso_c_binding
  use functions, only: decode, rebuild_state => buildState_Sz0_to_FullHS
  use exponentiate, only: diagHE
  use printing, only: printmat
  implicit none

  real (c_double), parameter, private :: pi = 4.d0 * datan(1.d0)

contains




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

  subroutine odd_reduced_DM(nspin, dim, dim_A, psi, rho_A)

    !Computes the reduced density matrix for a pure state of the odd spins in the chain
    ! rho_A = Tr_B(rho) = Tr_B( |psi> <psi| )
    !Traces out the remaining nspin/2 spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, dim, dim_A
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_A(dim_A,dim_A)

    integer :: jA, kA, iB, dim_B, i, j, configA1(nspin/2), configA2(nspin/2), configB(nspin/2), k

    i = nspin
    dim_B = dim/dim_A
    rho_A = 0
    do jA = 0, dim_A-1
      do kA = 0, dim_A-1
        do iB = 0, dim_B-1

          call decode(jA, nspin/2, configA1)
          call decode(kA, nspin/2, configA2)
          call decode(iB, nspin/2, configB)
          i = 0
          j = 0
          do k = 1, nspin/2
            i = i + configA1(k) * 2**(2*k-2) + configB(k) * 2**(2*k-1)
            j = j + configA2(k) * 2**(2*k-2) + configB(k) * 2**(2*k-1)
          enddo

          rho_A(jA+1,kA+1) = rho_A(jA+1,kA+1) + psi(i+1) * dconjg( psi(j+1) )

        enddo
      enddo
    enddo

  end subroutine

  subroutine even_reduced_DM(nspin, dim, dim_B, psi, rho_B)

    !Computes the reduced density matrix for a pure state of the even spins in the chain
    ! rho_A = Tr_B(rho) = Tr_B( |psi> <psi| )
    !Traces out the remaining nspin/2 spins
    !Works in the full Hilbert Space (with the correspondence j<->{b_k})

    integer, intent(in) :: nspin, dim, dim_B
    complex (c_double_complex), intent(in) :: psi(dim)
    complex (c_double_complex), intent(out) :: rho_B(dim_B,dim_B)

    integer :: jB, kB, iA, dim_A, i, j, configA(nspin/2), configB1(nspin/2), configB2(nspin/2), k

    i = nspin
    dim_A = dim/dim_B
    rho_B = 0
    do jB = 0, dim_B-1
      do kB = 0, dim_B-1
        do iA = 0, dim_A-1

          call decode(jB, nspin/2, configB1)
          call decode(kB, nspin/2, configB2)
          call decode(iA, nspin/2, configA)
          i = 0
          j = 0
          do k = 1, nspin/2
            i = i + configA(k) * 2**(2*k-2) + configB1(k) * 2**(2*k-1)
            j = j + configA(k) * 2**(2*k-2) + configB2(k) * 2**(2*k-1)
          enddo

          rho_B(jB+1,kB+1) = rho_B(jB+1,kB+1) + psi(i+1) * dconjg( psi(j+1) )

        enddo
      enddo
    enddo

  end subroutine




  function vN_entropy(dim, rho)

    !Computes the Entanglement (von Neumann) Entropy for a given density matrix rho of dimension dim, using Exact Diagonalization
    ! EE = S(rho) = -Tr(rho * ln(rho))
    ! For a reduced density matrix S_A = -Tr_A(rho_A ln(rho_A))

    integer, intent(in) :: dim
    complex (c_double_complex), intent(in) :: rho(dim,dim)
    real (c_double) :: vN_entropy

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
  
    vN_entropy = EE

  end function vN_entropy


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

    MI = vN_entropy(dim_A,rho_A) + vN_entropy(dim_A,rho_B) - vN_entropy(dim_A*dim_A,rho_AB)
    !print *, nspin_A, vN_entropy(dim_A,rho_A), vN_entropy(dim_A,rho_B), vN_entropy(dim_A*dim_A,rho_AB), MI
    
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

    call rebuild_state(nspin, dim, dim_Sz0, psi_Sz0, psi)

    call left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)
    call right_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_B)
    call edges_reduced_DM(nspin, nspin_A, nspin_A, dim, dim_A, dim_A, psi, rho_AB)
    
    !print *, "rho_A = "
    !call printmat(dim_A, rho_A, 'C')
    !print *, "rho_B = "
    !call printmat(dim_A, rho_B, 'C')
    !print *, "rho_AB = "
    !call printmat(dim_A*dim_A, rho_AB, 'C')

    MI = vN_entropy(dim_A,rho_A) + vN_entropy(dim_A,rho_B) - vN_entropy(dim_A*dim_A,rho_AB)

    !print *, "S_A         S_B          S_{AB}"
    !print *, nspin_A, vN_entropy(dim_A,rho_A), vN_entropy(dim_A,rho_B), vN_entropy(dim_A*dim_A,rho_AB), MI
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


    call rebuild_state(nspin, dim, dim_Sz0, psi_Sz0, psi)
    MBEE = 0
    pMBEE = 0
    !print *, "    nspin_A        dim_A     BEE"
    do nspin_A = 1, nspin - 1
      dim_A = 2**nspin_A
      if (nspin_A <= nspin/2) then
        allocate(rho_A(dim_A,dim_A))
        call left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)
        BEEvar = vN_entropy(dim_A,rho_A)
      else 
        allocate(rho_A(dim/dim_A,dim/dim_A))
        call right_reduced_DM(nspin, nspin-nspin_A, dim, dim/dim_A, psi, rho_A)
        BEEvar = vN_entropy(dim/dim_A,rho_A)
      endif
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

  function comb_entropy(nspin, dim, psi)

    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: psi(dim)
    real (c_double) :: comb_entropy

    integer (c_int) :: dim_A
    complex (c_double_complex), allocatable :: rho_A(:,:)
    real (c_double) :: CE!, CE2

    dim_A = 2**(nspin/2)
    allocate(rho_A(dim_A,dim_A))

    call odd_reduced_DM(nspin, dim, dim_A, psi, rho_A)
    CE = vN_entropy(dim_A, rho_A)

    !call even_reduced_DM(nspin, dim, dim_A, psi, rho_A)
    !CE2 = vN_entropy(dim_A, rho_A)

    deallocate(rho_A)

    comb_entropy = CE
    

  end function
    
  function comb_entropy_Sz0(nspin, dim, dim_Sz0, psi_Sz0)

    integer (c_int), intent(in) :: nspin, dim, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: comb_entropy_Sz0

    integer (c_int) :: dim_A
    complex (c_double_complex), allocatable :: rho_A(:,:), psi(:)
    real (c_double) :: CE!1, CE2

    allocate(psi(dim))
    call rebuild_state(nspin, dim, dim_Sz0, psi_Sz0, psi)

    dim_A = 2**(nspin/2)
    allocate(rho_A(dim_A,dim_A))

    call odd_reduced_DM(nspin, dim, dim_A, psi, rho_A)
    CE = vN_entropy(dim_A, rho_A)

    !call even_reduced_DM(nspin, dim, dim_A, psi, rho_A)
    !CE2 = vN_entropy(dim_A, rho_A)

    deallocate(rho_A)

    !print *, CE1, CE2
    comb_entropy_Sz0 = CE
    

  end function

  function avg_block_entropy_Sz0(nspin, dim, dim_Sz0, psi_Sz0)

    !Computes the average block entropy over the possible bipartitions A,B
    !A has left edge coinciding the left edge of the chain, and its right edge can vary over odd sites
    !rho_A -> left_reduced_DM
    !L_A = 1, 3, 5, ..., L-1

    integer, intent(in) :: nspin, dim, dim_Sz0
    complex (c_double_complex), intent(in) :: psi_Sz0(dim_Sz0)
    real (c_double) :: avg_block_entropy_Sz0

    real (c_double) :: avgBE
    integer (c_int) :: nspin_A, dim_A
    complex (c_double_complex), allocatable :: rho_A(:,:) !, rho(dim_Sz0,dim_Sz0)
    complex (c_double_complex) :: psi(dim)


    call rebuild_state(nspin, dim, dim_Sz0, psi_Sz0, psi)
    avgBE = 0
    !print *, "    nspin_A        dim_A     BEE"
    do nspin_A = 1, nspin - 1, 2
      dim_A = 2**nspin_A
      if (nspin_A <= nspin/2) then
        allocate(rho_A(dim_A,dim_A))
        call left_reduced_DM(nspin, nspin_A, dim, dim_A, psi, rho_A)
        avgBE = avgBE + vN_entropy(dim_A,rho_A)
      else 
        allocate(rho_A(dim/dim_A,dim/dim_A))
        call right_reduced_DM(nspin, nspin-nspin_A, dim, dim/dim_A, psi, rho_A)
        avgBE = avgBE + vN_entropy(dim/dim_A,rho_A)
      endif
      deallocate(rho_A)
    enddo
    avgBE = avgBE/(nspin/2)

    avg_block_entropy_Sz0 = avgBE

  end function

end module entanglement



