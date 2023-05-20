program test


  use sorts, only: sort => dpquicksort
  use functions, only: search => binsearch_closest_in_circle, binom
  use observables, only: log_gap_difference, pi_pair, &
    & exact_QE => exact_quasi_energies_Sz0_LR
  use iso_c_binding, only: dp => c_double, ip => c_int, c_double_complex
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = cmplx(0._dp, 1._dp, kind=c_double_complex)
  real (dp), parameter :: pi = 4*datan(1.0_dp)
  real (dp), allocatable :: array(:), pair(:), near(:), pair1(:), pair2(:)
  real (dp) :: val, log_pair_avg, log_near_avg, log_avg, log_sq
  integer (ip) :: i, n, alpha, beta, n_dis, dim, dim_Sz0
  integer (ip), allocatable :: pi_paired(:)
  integer (ip) :: beta1, betap, betas, beta1s, beta1p, beta2
  real (dp) :: log_avg1, log_avg2

  integer (ip) :: nspin
  real (dp), allocatable :: Vzz(:,:), hz(:)

  write (*,*) "Length of Chain:"
  read (*,*) nspin
  print *, ""
  dim_Sz0 = binom(nspin, nspin/2)
  n = dim_Sz0
  

  write (*,*) "Number of Iterations:"
  read (*,*) n_dis
  print *, ""

  allocate(array(n), pi_paired(n), pair(n), near(n))
  allocate(pair1(n), pair2(n))
  allocate(Vzz(nspin-1,nspin), hz(nspin))

  log_avg = 0
  log_avg1 = 0
  log_avg2 = 0
  do i = 1, n_dis

    !call random_number(array)
    !array = 2*pi*(array - 0.5)
    !array(n/2+1) = array(n/2)
    !array(n) = mod(array(n/2-1) + 2*pi,2*pi) - pi

    call random_number(Vzz)
    call random_number(hz)

    array = exact_QE(nspin, dim_Sz0, Vzz, hz)
    call sort(array)

    !print "(1(A12),1(A22,4X))", "i", "a(i)"
    !do i = 1, size(array)
    !  print *, i, array(i)
    !enddo
    !print *, ""
    !print *, ""

    pi_paired = pi_pair(array)

    do alpha = 1, n

      val = mod( array(alpha) + 2*pi, 2*pi) - pi

      beta = merge(alpha+1,1,alpha.ne.n)
      near(alpha) = array(beta) + merge(0._dp,2*pi,alpha.ne.n) - array(alpha)

      beta = pi_paired(alpha)
      beta1 = minloc( abs(exp(C_UNIT*(array-val)) - 1), n )
      beta2 = merge(alpha+n/2, alpha-n/2, alpha<=n/2)

      betap = merge(beta-1,n, beta .ne. 1)
      betas = merge(beta+1,1, beta .ne. n)

      beta1p = merge(beta1-1,n, beta1 .ne. 1)
      beta1s = merge(beta1+1,1, beta1 .ne. n)

      pair(alpha) = abs(abs(array(beta) - array(alpha)) - pi)
      pair1(alpha) = abs(abs(array(beta1) - array(alpha)) - pi)
      pair2(alpha) = abs(abs(array(beta2) - array(alpha)) - pi)

      !print *, alpha, beta, beta1, beta2
      !print *, array(alpha), val, array(beta), array(beta1), array(beta2)

      !print *, pair(alpha), pair1(alpha), pair2(alpha)
      !print *, ""

    enddo
    pair = max(pair,epsilon(pair))
    pair1 = max(pair1,epsilon(pair1))
    pair2 = max(pair2,epsilon(pair2))

    log_avg = log_avg + sum(log(pair))/n
    log_avg1 = log_avg1 + sum(log(pair1))/n
    log_avg2 = log_avg2 + sum(log(pair2))/n

  enddo

  log_avg = log_avg/n_dis
  log_avg1 = log_avg1/n_dis
  log_avg2 = log_avg2/n_dis

  print *, "log_dis_avg = ", log_avg, log_avg1, log_avg2
  
end program
