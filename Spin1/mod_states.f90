!This module contains various procedures used for building specific initial states:
! One-spin product state (|phi>)^(\otimes L)
! Neel state (up, down, up, down, ...)
! Class of LI=1, -1<I<1 state ( a|up,down> + b|down,up> + c|0,0> ) \otimes ...
! Class of 0<LI<1, I=0 state ( a|up,down> + b|down,up> + c|0,0> ) \otimes (<-> swap) \otimes ...

!There are variants for both the Hilbert space and the Sz subspaces (specifically Sz=0)


module states


  use iso_c_binding
  use printing
  use functions
  implicit none

  integer (c_int), private, parameter :: dimSpin1 = 3
  real (c_double), private, parameter :: tol = 1.0e-6

contains

  subroutine buildProdState(nspin, dim, alpha, state)

    !Builds a product state of the form: |phi> \otimes |phi> \otimes ...
    ! where: |phi> = alpha_1|+1> + alpha_2|0> + alpha_3|-1>
 
    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: alpha(dimSpin1)
    complex (c_double_complex), intent(out) :: state(dim)
    integer :: config(nspin), num(dimSpin1)
    complex (c_double_complex) :: alpha_n(dimSpin1)
    real (c_double) :: norm
    integer (c_int) :: i, k, idx

    alpha_n = alpha/dot_product(alpha,alpha)
    !print *, alpha_n

    !num() contains the number of all spins with a certain orientation
    ! (up, zero, down) <-> (0,1,2) in a given configuration (encoded in an integer i)
    ! It counts the number of 0, 1, 2 in 'config'
    state = 1
    do i = 1, dim
      call decode(i-1,nspin,config)
      !print *, "config = ", config(:)
      num = 0
      do k = 1, nspin
        num(config(k)+1) = num(config(k)+1) + 1
      enddo
      !print *, "num = ", num(:)

      do idx = 1, dimSpin1
        state(i) = state(i) * alpha_n(idx) ** num(idx)
        !print *, "alpha_n(", idx ,") =", alpha_n(idx) ** num(idx)
        !print *, "state(i) = ", state(i)
      enddo
      !print *, "state = ", state(i)

    enddo

  end subroutine buildProdState

!  subroutine buildNayakState(nspin,dim,alpha,beta,state)
!
!    integer (c_int), intent(in) :: nspin, dim
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim)
!    
!    integer :: config(nspin), i, k
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    do i = 1, dim
!
!      call decode(i-1,nspin,config)
!
!      do k = 1, nspin/2
!        state(i) = state(i) * ( (1-config(2*k-1)) * alpha_n + config(2*k-1) * beta_n ) * &
!          & ( config(2*k) * alpha_n + (1-config(2*k)) * beta_n )
!      enddo
!    enddo
!
!
!  end subroutine buildNayakState
!
!  subroutine buildNayakState_Sz0(nspin,dim_Sz0,alpha,beta,state)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim_Sz0)
!    
!    integer :: config(nspin), i, k, l, indx(dim_Sz0)
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    call zero_mag_states(nspin, dim_Sz0, indx)
!    do i = 1, dim_Sz0
!      l = indx(i)
!      call decode(l,nspin,config)
!
!      do k = 1, nspin/2
!        state(i) = state(i) * ( (1-config(2*k-1)) * alpha_n + config(2*k-1) * beta_n ) * &
!          & ( config(2*k) * alpha_n + (1-config(2*k)) * beta_n )
!      enddo
!    enddo
!
!
!  end subroutine buildNayakState_Sz0
!
!
!  subroutine buildI0LI1ProdState_Sz0(nspin, dim_Sz0, alpha, beta, state)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim_Sz0)
!    
!    integer (c_int) :: config(nspin), i, k, l, indx(dim_Sz0)
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    call zero_mag_states(nspin, dim_Sz0, indx)
!    do i = 1, dim_Sz0
!      l = indx(i)
!      call decode(l,nspin,config)
!
!      do k = 1, nspin/4
!        state(i) = state(i) * ( (1-config(4*k-3)) * config(4*k-2) * alpha_n + &
!          &  config(4*k-3) * (1-config(4*k-2)) * beta_n ) * &
!          & ( (1 - config(4*k-1)) * config(4*k) * beta_n + &
!          &  config(4*k-1) * (1-config(4*k)) * alpha_n )
!      enddo
!      if(mod(nspin/2,2) /= 0) then
!        k = nspin/4 + 1
!        state(i) = state(i) * ( (1-config(4*k-3)) * config(4*k-2) * alpha_n + &
!          &  config(4*k-3) * (1-config(4*k-2)) * beta_n )
!      endif
!    enddo
!
!
!  end subroutine
!
  subroutine buildNeelState_Sz0(nspin, dim_Sz0, psi)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
    
    integer :: config(nspin), i, k, indx(dim_Sz0), l, inverse(dimSpin1**nspin)

    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
    call zero_mag_states_inv(nspin, dim_Sz0, indx, inverse)

    !Neel state is given by up, down, up, down, ...
    ! which (config = 1 - s) corresponds to 0, 2, 0, 2, ...
    i = 0
    do k = 1, nspin/2
      config(2*k-1) = 0
      config(2*k) = 2
      i = i + config(2*k-1) * dimSpin1**(2*k-2) + config(2*k) * dimSpin1**(2*k-1)
    enddo
    
    psi = 0
    psi(inverse(i+1)) = 1

  end subroutine buildNeelState_Sz0


  subroutine buildNeelState(nspin, dim, psi)

    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(out) :: psi(dim)
    
    integer :: config(nspin), i, k

    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even."
    if (dim /= dimSpin1 ** nspin ) stop "Error: Wrong dimension."

    !Neel state is given by up, down, up, down, ...
    ! which (config = 1 - s) corresponds to 0, 2, 0, 2, ...
    i = 0
    do k = 1, nspin/2
      config(2*k-1) = 0
      config(2*k) = 2
      i = i + config(2*k-1) * dimSpin1**(2*k-2) + config(2*k) * dimSpin1**(2*k-1)
    enddo
    
    psi = 0
    psi(i+1) = 1

  end subroutine buildNeelState

  subroutine buildUpZeroState(nspin, dim, psi)

    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(out) :: psi(dim)
    
    integer :: config(nspin), i, k

    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even."
    if (dim /= dimSpin1 ** nspin ) stop "Error: Wrong dimension."

    !Neel state is given by up, zero up, zero, ...
    ! which (config = 1 - s) corresponds to 0, 1, 0, 1, ...
    i = 0
    do k = 1, nspin/2
      config(2*k-1) = 0
      config(2*k) = 1
      i = i + config(2*k-1) * dimSpin1**(2*k-2) + config(2*k) * dimSpin1**(2*k-1)
    enddo
    
    psi = 0
    psi(i+1) = 1

  end subroutine

!
!  subroutine buildLI1ProdState_Sz0(nspin, dim_Sz0, alpha, beta, state)
!
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    complex (c_double_complex), intent(in) :: alpha, beta
!    complex (c_double_complex), intent(out) :: state(dim_Sz0)
!    
!    integer (c_int) :: config(nspin), i, k, l, indx(dim_Sz0)
!    complex (c_double_complex) :: alpha_n, beta_n
!    real (c_double) :: norm
!
!    if (mod(nspin,2)==1) stop "Error: Number of Spins must be even"
!
!    norm = sqrt(abs(alpha)**2 + abs(beta)**2)
!    alpha_n = alpha/norm
!    beta_n = beta/norm
!    state = 1
!    call zero_mag_states(nspin, dim_Sz0, indx)
!    do i = 1, dim_Sz0
!      l = indx(i)
!      call decode(l,nspin,config)
!
!      do k = 1, nspin/2
!        state(i) = state(i) * ( (1-config(2*k-1)) * config(2*k) * alpha_n + &
!          &  config(2*k-1) * (1-config(2*k)) * beta_n )
!      enddo
!
!    enddo
!
!  end subroutine
!
!
!  subroutine buildRndLarge_IMB_LI_State_Sz0(nspin, dim_Sz0, IMB, LI, psi)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: IMB, LI
!    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
!
!    integer (c_int) :: idxSz0(dim_Sz0), nz, l, i, config(nspin)
!    integer (c_int), allocatable :: idx(:)
!    real (c_double) :: rand
!
!    call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, idxSz0)
!
!    do i = 1, dim_Sz0
!      if(idxSz0(i) == 0) then
!        nz = i-1
!        exit
!      endif
!    enddo
!
!    allocate(idx(nz))
!
!    idx = idxSz0(1:nz)
!
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    call init_random_seed()
!    do i = 1, nz
!      l = idx(i)
!      call random_number(rand)
!      psi(binsearch(l,idxSz0)) = rand
!    enddo
!    psi = psi/sqrt(dot_product(psi,psi))
!
!    !print *, imbalance_Sz0(nspin, dim_Sz0, psi), local_imbalance_Sz0(nspin, dim_Sz0, psi)
!
!  end subroutine
!
!
!  subroutine buildRndLarge_IMB_LI_State_basis_Sz0(nspin, dim_Sz0, IMB, LI, psi)
!    integer (c_int), intent(in) :: nspin, dim_Sz0
!    real (c_double), intent(in) :: IMB, LI
!    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
!
!    integer (c_int) :: idxSz0(dim_Sz0), nz, l, i, config(nspin)
!    integer (c_int), allocatable :: idx(:)
!    real (c_double) :: rand
!
!    call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, idxSz0)
!
!    do i = 1, dim_Sz0
!      if(idxSz0(i) == 0) then
!        nz = i-1
!        exit
!      endif
!    enddo
!
!    allocate(idx(nz))
!
!    idx = idxSz0(1:nz)
!
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    call init_random_seed()
!    call random_number(rand)
!    l = idx(floor(nz*rand + 1))
!    i = binsearch(l,idxSz0)
!
!    psi = 0
!    psi(i) = 1
!
!    call decode(l, nspin, config)
!    !print *, IMB, LI
!    !print *, imbalance_basis(nspin,l), local_imbalance_basis(nspin,l), l
!    1 format ( 4X,F6.3, 4X,F6.3, 4X,I4, 4X,I0 )
!
!
!  end subroutine
!    
!  subroutine buildLarge_IMB_LI_State_basis_Sz0(nspin, dim_Sz0, i, IMB, LI, psi)
!    integer (c_int), intent(in) :: nspin, dim_Sz0, i
!    real (c_double), intent(in) :: IMB, LI
!    complex (c_double_complex), intent(out) :: psi(dim_Sz0)
!
!    integer (c_int) :: idxSz0(dim_Sz0), nz, j, l, config(nspin)
!    integer (c_int), allocatable :: idx(:)
!    real (c_double) :: rand
!
!    call large_IMB_LI_states_Sz0(nspin, dim_Sz0, IMB, LI, idxSz0)
!
!    do j = 1, dim_Sz0
!      if(idxSz0(j) == 0) then
!        nz = j-1
!        exit
!      endif
!    enddo
!    print *, "nz = ", nz, "i = ", i
!    if (i>nz) stop "Error: Not enough number of states."
!
!    l = idxSz0(i)
!    call zero_mag_states(nspin, dim_Sz0, idxSz0)
!    j = binsearch(l, idxSz0)
!    psi = 0
!    psi(j) = 1
!
!    call decode(l, nspin, config)
!    print *, "build Large IMB, LI basis state"
!    print *, IMB, LI
!    print *, imbalance_basis(nspin,l), local_imbalance_basis(nspin,l), l
!    1 format ( 4X,F6.3, 4X,F6.3, 4X,I4, 4X,I0 )
!
!  end subroutine
!
!

  subroutine printstate(nspin, dim, state, state_name)

    integer (c_int), intent(in) :: nspin, dim
    complex (c_double_complex), intent(in) :: state(dim)
    character(len=*), intent(in) :: state_name

    integer :: i, config(nspin)

    print *, state_name
    do i = 0, dim-1
      if (abs(state(i+1))**2 > tol) then
        call decode(i, nspin, config)
        print "( 4X,F8.4, 4X,I6, 4X,*(I0) )", abs(state(i+1))**2, i, config(:)
      endif
    enddo
    print *, ""

  end subroutine

  subroutine printstate_Sz0(nspin, dim_Sz0, state, state_name)

    integer (c_int), intent(in) :: nspin, dim_Sz0
    complex (c_double_complex), intent(in) :: state(dim_Sz0)
    character(len=*), intent(in) :: state_name

    integer :: i, l, config(nspin), idxSz0(dim_Sz0)

    print *, state_name
    call zero_mag_states(nspin, dim_Sz0, idxSz0)
    do l = 1, dim_Sz0
      if (abs(state(l))**2 > tol) then
        i = idxSz0(l)
        call decode(i, nspin, config)
        print "( 4X,F8.4, 2(4X,I6), 4X,*(I0) )", abs(state(l))**2, l, i, config(:)
      endif
    enddo
    print *, ""

  end subroutine

  subroutine printstate_Sz(nspin, dim_Sz, Sz, state, state_name)

    integer (c_int), intent(in) :: nspin, dim_Sz, Sz
    complex (c_double_complex), intent(in) :: state(dim_Sz)
    character(len=*), intent(in) :: state_name

    integer :: i, l, config(nspin), idxSz(dim_Sz)

    print *, state_name
    call basis_Sz(nspin, dim_Sz, Sz, idxSz)
    do l = 1, dim_Sz
      if (abs(state(l))**2 > tol) then
        i = idxSz(l)
        call decode(i, nspin, config)
        print "( 4X,F8.4, 2(4X,I6), 4X,*(I0) )", abs(state(l))**2, l, i, config(:)
      endif
    enddo
    print *, ""

  end subroutine


end module states
