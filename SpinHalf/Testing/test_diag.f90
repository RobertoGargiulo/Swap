program diagonalization

  use iso_c_binding, dp => c_double, ip => c_int, dcp => c_double_complex
  use exponentiate, only: diagSYM, expSYM, diagUN, diagHE, expHE
  use printing, only: printmat
  use sorts, only: sort => dpquicksort
  implicit none


  complex (dcp), parameter :: C_UNIT = dcmplx(0._dp, 1._dp)
  real (dp), parameter :: pi = 4._dp * atan(1._dp)
  complex (dcp) :: M(2,2), W(2,2), U(2,2), PH(2)
  real (dp) :: E(2), QE(2)
  integer :: i, j

  logical :: SELECT
  EXTERNAL SELECT

  M(1,1) = 1
  M(1,2) = -C_UNIT
  M(2,1) = C_UNIT
  M(2,2) = -1

  do i = 1, 2
    do j = 1, 2
      print *, i, j, M(i,j)
    enddo
  enddo

  call diagHE('V', 2, M, E, W)

  print *, "E = ", E(:)

  call expHE( 2, 1._dp, E, W, U )
  call diagUN( SELECT, 2, U, PH, W)

end program 

  

logical function SELECT(z)

  implicit none

  complex(kind=8), intent(in) :: z

  SELECT = .TRUE.
  RETURN

end



