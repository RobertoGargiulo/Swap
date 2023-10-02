program dummy

  use omp_lib
  use ifport
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  implicit none

  integer (ip)     ::  dim
  complex(dcp), dimension(:,:), allocatable :: U

  write (*,*) "dim"
  read (*,*) dim
  print*,""

  print *, "Allocating USwap."
  allocate(U(dim,dim))

  print *, "dim = ", dim, "dim^2 = ", dim**2, "size(U) = ", size(U)
  print *, "Building U..."
  U = 0
  print *, "USwap initialized (set to zero)."

  print *, "Testing matrix multiplication matmul(USwap, USwap)"
  U = matmul(U,U) 
  print *, "USwap built."

  !$OMP PARALLEL
  call init_random_seed()
  print *, "Size of Thread team: ", omp_get_num_threads()
  print *, "Current code segment is in parallel: ", omp_in_parallel()
  !$OMP do private(i, U)
  do i = 1, n_disorder
    call 

  !$OMP END DO
  !$OMP END PARALLEL

end program

