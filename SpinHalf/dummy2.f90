program dummy

  use ifport
  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  implicit none

  integer (ip)     ::  dim
  complex(dcp), dimension(:,:), allocatable :: U

  write (*,*) "dim"
  read (*,*) dim
  print*,""

  print *, "Allocating U."
  allocate(U(dim,dim))

  print *, "dim = ", dim, "dim^2 = ", dim**2, "size(U) = ", size(U)
  print *, "Building U..."
  U = 0
  print *, "U initialized (set to zero)."

  !print *, "U = ", U

  print *, "Testing matrix multiplication matmul(U, U)"
  U = matmul(U,U)
  print *, "U built."

  !print *, "U = ", U


end program
