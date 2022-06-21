program integerkinds
    use iso_c_binding
    use omp_lib
    implicit none

    integer :: i, psum, sum

!$OMP PARALLEL PRIVATE(psum) SHARED(sum)
    
    psum = 0
    sum = 0

    !$OMP DO
    do i = 1, 5
        psum = psum + i
        print *, "Thread: ", omp_get_thread_num(), "psum = ", psum, i
    enddo
    !$OMP END DO

    !$OMP CRITICAL
        sum = sum + psum
    !$OMP END CRITICAL

!$OMP END PARALLEL
    print *, "Thread :", omp_get_thread_num(), "sum = ", sum, "i = ", i

end program integerkinds

