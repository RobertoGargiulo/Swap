program test


  use sorts, only: sort => dpquicksort
  use functions, only: search => binsearch_closest_in_circle
  use observables, only: log_gap_difference, pi_pair
  use iso_c_binding, only: dp => c_double, ip => c_int, c_double_complex
  implicit none

  complex (c_double_complex), parameter :: C_UNIT = cmplx(0._dp, 1._dp)
  real (dp), parameter :: pi = 4*datan(1.0_dp)
  real (dp), allocatable :: array(:), pair(:), near(:)
  real (dp) :: val, log_pair_avg, log_near_avg, log_avg, log_sq
  integer (ip) :: i, n, alpha, beta, n_dis
  integer (ip), allocatable :: pi_paired(:)

  write (*,*) "Size of array:"
  read (*,*) n
  print *, ""

  write (*,*) "Number of Iterations:"
  read (*,*) n_dis
  print *, ""

  allocate(array(n), pi_paired(n), pair(n), near(n))

  do i = 1, n_dis

    call random_number(array)
    array = 2*pi*(array - 0.5)
    array(n/2+1) = array(n/2)
    array(n) = mod(array(n/2-1) + 2*pi,2*pi) - pi
    call sort(array)

    !print "(1(A12),1(A22,4X))", "i", "a(i)"
    !do i = 1, size(array)
    !  print *, i, array(i)
    !enddo
    !print *, ""
    !print *, ""

    !call log_gap_difference(n, array, log_pair_avg, log_near_avg, log_avg, log_sq)
    !print *, "log_avg = ", log_avg

    !pi_paired = pi_pair(n, array)

    do alpha = 1, n 

      beta = merge(alpha+1,1,alpha.ne.n)
      near(alpha) = array(beta) + merge(0._dp,2*pi,alpha.ne.n) - array(alpha)

      !beta = pi_paired(alpha)
      beta = minloc( abs(exp(C_UNIT*(array-val)) - 1), n )
      pair(alpha) = abs(abs(array(beta) - array(alpha)) - pi)

      !print *, pair(alpha), near(alpha), log(pair(alpha)), log(near(alpha))

    enddo

    log_pair_avg = sum(log(pair)) / n
    log_near_avg = sum(log(near)) / n
    log_avg = log_pair_avg - log_near_avg
    log_sq = sum((log(pair) - log(near))**2) / n

  enddo

  
end program
