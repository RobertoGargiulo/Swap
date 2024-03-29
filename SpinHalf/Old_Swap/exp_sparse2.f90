MODULE exp_sparse
!USE factorials
!Important: requires to be built with mataid.o expokit.o
use iso_c_binding
USE MATAID
USE EXPOKIT
implicit none
!external zgcoov
COMPLEX (c_double_complex), PRIVATE, PARAMETER :: ZERO=DCMPLX(0.0d0,0.0d0)
COMPLEX (c_double_complex), PRIVATE, PARAMETER :: C_ONE=DCMPLX(1.0d0,0.0d0)
COMPLEX (c_double_complex), PRIVATE, PARAMETER :: C_UNIT=DCMPLX(0.0d0,1.0d0)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine that calculates w=exp(t*A)*v with v and A complex
!A is nxn and is assumed to be sparse
!A has nz non-zero entries, each of one specified in Coordinate List (COO) format in (ia, ja, a) which are respectively (row, column, A(row,column))
!---Arguments------------------
!	dim         : (input) order of the principal matrix A.
!	nz_dim      : (input) number of non-zero entries of A
!	k_dim       : (input) maximum size for the Krylov basis.
!	rows	      : (input) rows of non-zero entries of A
!	cols	      : (input) columns of non-zero entries of A
!	H	          : (input) non-zero entries of A
!	init_state  : (input) given operand vector (see above)7
!	t           : (input) time at wich the solution is needed (can be < 0).
!	fin_state	  : (output) computed approximation of exp(t*H)*v.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE evolveKrylov(dim, nz_dim, k_dim, ia, ja, a, v, t, w)
	integer, intent(in):: n, nz, m
	integer, intent(in):: ia(nz), ja(nz)
	complex*16, intent(in):: a(nz), v(n)
	complex*16, intent(out) :: w(n)
	real*8, intent(in):: t
	integer, dimension(:), allocatable:: iwsp
	double precision tol, anorm
	complex*16, dimension(:), allocatable::  wsp
	integer i, itrace, iflag
	complex*16 ZERO, ONE
	parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
	!common /CMAT/ a, ia, ja, nz, n
	!	  arguments variables ...
	integer lwsp, liwsp
	lwsp = n*(m+2)+5*(m+2)**2+7
	liwsp = m+2
	allocate(iwsp(liwsp))
	allocate(wsp(lwsp))
	
	do i = 1,n
		wsp(i) = ZERO
	enddo
	do i = 1,nz
		wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
	enddo
	anorm = wsp(1)
	do i = 2,n
		if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
	enddo

	
	
	!intrinsic ABS, CMPLX, CONJG, DBLE
	
	tol = 1.0d-7
    itrace = 0
!    WRITE(*,*) 'Orcomeno'
    call ZGEXPV( n,a,ia,ja,nz, m, t,v,w, tol, anorm, wsp,lwsp, iwsp,liwsp, zgcoov, itrace, iflag )
     
end subroutine evolve


subroutine ZEXPM( n, m, t, A,LDA, B,LDB, C,LDC )
      implicit none
      double precision t
      integer n, m, LDA, LDB, LDC
      double complex A(LDA,n), B(LDA,m), C(LDC,m)

      integer ideg, iflag, iexp, ns, lwsp
!      double complex ZERO, ONE
!      integer :: ideg=6, lwsp=4*n*n+ideg+1
!      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

      double complex, allocatable :: wsp(:)
      integer, allocatable :: ipiv(:)

      ideg = 6
      lwsp = 4*n*n+ideg+1
      ALLOCATE( wsp(lwsp), ipiv(n) )

!*---  compute E = exp(tA) via Expokit's Pade algorithm
      call ZGPADM(ideg, n, t, A,LDA, wsp,lwsp, ipiv, iexp, ns, iflag)

!*---  multiply C = E*B via BLAS
      call ZGEMM('n','n', n,m,n, ONE, wsp(iexp),n, B,LDB, ZERO, C,LDC )

      DEALLOCATE( ipiv, wsp )
END subroutine ZEXPM
!
!SUBROUTINE ZESPONENZO(n, m, t, A,LDA, B,LDB, C,LDC, iord)
!
! IMPLICIT NONE
! INTEGER :: n,m
! DOUBLE PRECISION :: t
 
!
! aux = A
!
! baux = A
!
! DO ind = 1, n
!  baux = ONE + baux(ind,ind)
! ENDDO
!
! DO ind = 1, iord-1
!   call ZGEMM('n','n', n,n,n, ONE,A,n, aux,n, ZERO, aux,n )
!   baux = baux + aux / DEXP(facln(ind+1))
! ENDDO
!
! CALL ZGEMM('n','n', n,m,n, baux, wsp(iexp),n, B,LDB, ZERO, C,LDC )
!
!END SUBROUTINE ZESPINENZO

END MODULE exp_sparse
