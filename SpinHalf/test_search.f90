program test


  use sorts, only: sort => dpquicksort
  use functions, only: search => binsearch_closest_in_circle
  use observables, only: log_gap_difference
  use iso_c_binding, only: dp => c_double, ip => c_int
  implicit none

  real (dp), parameter :: pi = 4*datan(1.0_dp)
  real (dp), allocatable :: array(:)
  real (dp) :: val, log_pair_avg, log_near_avg, log_avg, log_sq
  integer (ip) :: i, n, alpha, beta
  integer (ip), allocatable :: pi_pair(:)

  write (*,*) "Size of array:"
  read (*,*) n
  print *, ""

  allocate(array(n), pi_pair(n))

  call random_number(array)
  array = 2*pi*(array - 0.5)
  array(n/2+1) = array(n/2)
  array(n) = mod(array(n/2-1) + 2*pi,2*pi) - pi
  call sort(array)

  print "(1(A12),1(A22,4X))", "i", "a(i)"
  do i = 1, size(array)
    print *, i, array(i)
  enddo
  print *, ""
  print *, ""

  call log_gap_difference(n, array, log_pair_avg, log_near_avg, log_avg, log_sq)
  print *, "log_avg = ", log_avg
  
  !do alpha = 1, n
  !  val = mod(array(alpha) + 2*pi, 2*pi) - pi
  !  if (val > array(n)) then 
  !    val = val - 2*pi
  !    !print *, "Problem: alpha = ", alpha, "; a(alpha) = ", array(alpha)
  !    !print *, ""
  !  endif

  !  beta = search(val, array)
  !  pi_pair(alpha) = beta
  !  !print "(2(A12),4(A22,4X))", "i", "alpha", "a(i)", "val", "(a(alpha) + pi)_1", "a(alpha)"
  !  !do i = 1, size(array)
  !  !  print *, i, alpha, array(i), val, mod(array(alpha) + 2*pi, 2*pi) - pi, array(alpha)
  !  !enddo

  !  !print "(2(A12),4(A22,4X))", "beta", "alpha", "a(beta)", "val", "(a(alpha) + pi)_1", "a(alpha)"
  !  !print *, beta, alpha, array(beta), val, mod(array(alpha) + 2*pi, 2*pi) - pi, array(alpha)

  !  !print *, "Check: "
  !  !print "(1(A12),1(A22,4X))", "i", "a(i) - val"!, "|a(i) - a(alpha)|-pi"
  !  !do i = 1, size(array)
  !  !  print *, i, array(i) - val!, abs(array(i) - array(alpha)) - pi
  !  !enddo
  !  !print *, ""
  !  !print *, ""
  !  !print *, ""

  !  beta = merge(alpha+n/2,alpha-n/2,alpha+n/2<=n)
  !  print *, alpha, pi_pair(alpha), beta, array(pi_pair(alpha))-array(alpha), &
  !   & array(pi_pair(alpha)) - val, abs(abs(array(pi_pair(alpha)) - array(alpha)) - pi), &
  !   & abs(abs(array(beta) - array(alpha)) - pi)

  !enddo

  !print *, "Pairs: "
  !print "(3(A12),1(A22,4X))", "alpha", "beta", "alpha + n/2", "a(beta) - a(alpha)"!, "|a(i) - a(alpha)|-pi"
  !do alpha = 1, n
  !  print *, alpha, pi_pair(alpha), merge(alpha+n/2,alpha-n/2,alpha+n/2<=n), array(pi_pair(alpha))-array(alpha)
  !enddo

end program
