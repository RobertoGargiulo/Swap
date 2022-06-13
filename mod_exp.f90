module exponentiate

  use iso_c_binding !Necessariamente nel modulo? Oppure basta nel programma?

  implicit none

  complex (c_double_complex), private, parameter :: C_ZERO = dcmplx(0._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_ONE = dcmplx(1._c_double, 0._c_double)
  complex (c_double_complex), private, parameter :: C_UNIT = dcmplx(0._c_double, 1._c_double)

  INTEGER, PRIVATE, PARAMETER              :: idp=KIND(1.0D0)
  REAL(KIND=idp), PRIVATE, PARAMETER       :: ZERO = 0.0D0, ONE = 1.0D0
  REAL(KIND=idp), PRIVATE, PARAMETER       :: Pi = 3.141592653589793D0
  COMPLEX(KIND=idp), PRIVATE, PARAMETER    :: II = DCMPLX(0.D0,1.D0)
  COMPLEX(KIND=idp), PRIVATE, PARAMETER    :: ZZERO = DCMPLX(0.D0,0.D0)


contains

!  subroutine allocateWork(N)
!
!    integer :: N !Dimensione spazio di Hilbert, N = dim = 2**nspin
!
!    allocate(iwork(5*N), ifail(N), work(8*N), HP(N*(N+1)/2))
!  
!  end subroutine allocateWork
!
!
!  subroutine deallocateWork(N)
!
!    integer :: N !Dimensione spazio di Hilbert, N = dim = 2**nspin
!
!    deallocate(iwork(5*N), ifail(N), work(8*N), HP(N*(N+1)/2))
!  
!  end subroutine deallocateWork


  subroutine diagSYM( JOBZ, N, H, E, W )
    
    !Computes all eigenvalues and eigenvectors of H(N,N), and puts them in E(N), W(N,N)
    implicit none

    character, intent(in) :: JOBZ*1
    integer, intent(in) :: N !Dimension of Hilbert Space
    real(c_double), intent(in) :: H(N,N)
    real(c_double), intent(out) :: E(N), W(N,N) !Siccome H non ci serve più, non possiamo usare solo H come matrice iniziale
    !e poi ci inseriamo gli autovettori?

    real(c_double) :: err
    
    integer :: iwork(5*N), ifail(N), info
    real(c_double) :: work(8*N), HP(N*(N+1)/2), VL, VU, DLAMCH

    integer :: i, j, M
    !Variabili intermedie necessarie per DSPEVX. Eliminarle tramite (de)allocateWork?
    M = N

    do j = 1, n
      do i = 1, j
        HP(i + j*(j-1)/2) = H(i,j)
        !print*, i, j, i + j*(j-1)/2, n*(n+1)/2
      enddo
    enddo
    !print *, "H in packed storage HP"

    !Check per simmetria?

    call DSPEVX(JOBZ, 'A', 'U', N, HP, VL, VU, 1, N, 2.d0*DLAMCH('S'), M, E, W, N, work, iwork, ifail, info)
    !DSPEVX(char, char, char, int, double(:), double, int, int, double, int, double(:), double(:,:), int, double(:), int(:), int(:),
    !int) 

    if (info /= 0) then
      print *, 'Diagonalizzazione Fallita. INFO = ', info
      print *, ""
      print "(1X,A3,3X,A5)", 'i','IFAIL'
      do i = 1, N
        print "(1X,I3,3X,I5)", i, ifail(i)
      enddo
    endif

  end subroutine diagSYM


  subroutine expSYM( dim, const, E, W, U )

    integer, intent(in) :: dim
    real(c_double), intent(in) :: E(dim), W(dim,dim)
    complex(c_double_complex), intent(in) :: const
    complex(c_double_complex), intent(out) :: U(dim,dim)

    complex(c_double_complex) :: Udiag(dim,dim), Uaux(dim,dim)
    integer :: i

    Udiag = 0
    do i = 1, dim
      Udiag(i,i) = exp(const*dcmplx(E(i)))
    enddo

    !U = matmul(W,matmul(Udiag,transpose(W)))

    call zgemm('N','C', dim, dim, dim, C_ONE, Udiag, dim , dcmplx(W), dim, C_ZERO, Uaux, dim)
    call zgemm('N','N', dim, dim, dim, C_ONE, dcmplx(W), dim , Uaux, dim, C_ZERO, U, dim)


  end subroutine expSYM

  subroutine diagGENERAL( JOBZ, N, U, PH, WL, WR )
    
    !Computes all eigenvalues and eigenvectors of U(N,N), and puts them in E(N), W(N,N)
    implicit none

    character, intent(in) :: JOBZ*1
    integer, intent(in) :: N !Dimension of Hilbert Space
    complex(c_double_complex), intent(in) :: U(N,N)
    complex(c_double_complex), intent(out) :: PH(N), WR(N,N), WL(N,N) !Siccome H non ci serve più, non possiamo usare solo H come matrice iniziale
    !e poi ci inseriamo gli autovettori?

    integer :: info
    complex(c_double_complex) :: work(8*N), Uaux(N,N)
    real(c_double) :: rwork(2*N)

    integer :: i
    !Variabili intermedie necessarie per DSPEVX. Eliminarle tramite (de)allocateWork?
    Uaux = U
    !Check per simmetria?
    call ZGEEV( JOBZ , JOBZ, N, Uaux, N, PH, WL, N, WR, N, work, 8*N, rwork, info  )
    !ZGEEV( char*1 JOBVL , char* JOBVR, N, dcmplx A(LDA,N), int LDA, dcmplx W(N), dcmplx VL(LDVL,N), int LDVL, dcmplx VR(LDVR,N),
    !int LDVR, dcmplx WORK(LWORK), int LWORK, dcmplx RWORK(2*N), int INFO )

    !call DSPEVX(JOBZ, 'A', 'U', N, HP, VL, VU, 1, N, 2.d0*DLAMCH('S'), M, E, W, N, work, iwork, ifail, info)
    !DSPEVX(char, char, char, int, double(:), double, int, int, double, int, double(:), double(:,:), int, double(:), int(:), int(:),
    !int) 

    if (info /= 0) then
      print *, 'Diagonalizzazione Fallita. INFO = ', info
      print *, ""
    endif

  end subroutine diagGENERAL


  subroutine diagUN( SELECT, N, U, PH, W )
    
    !Computes all eigenvalues and eigenvectors of U(N,N), and puts them in E(N), W(N,N)
    implicit none

    integer, intent(in) :: N !Dimension of Hilbert Space
    complex(c_double_complex), intent(in) :: U(N,N)
    complex(c_double_complex), intent(out) :: PH(N), W(N,N) !Siccome H non ci serve più, non possiamo usare solo H come matrice iniziale
    !e poi ci inseriamo gli autovettori?

    integer :: info, SDIM
    complex(c_double_complex) :: work(3*N), Uaux(N,N)
    real(c_double) :: rwork(N)

    integer :: i
    logical :: bwork(N), SELECT
    external :: SELECT
    Uaux = U
    call ZGEES( 'V', 'S', SELECT, N, Uaux, N, SDIM, PH, W, N, work, 3*N, rwork, bwork, info)

    !call ZGEES( char*1 JOBVS, char* SORT, logical SELECT, int N, dcmplx A(LDA,N), int LDA, int SDIM, dcmplx W(N), dcmplx VS(LDVS,N),
    !int LDVS, dcmplx WORK(LWORK), int LWORK(3*N), dreal RWORK(N), logical BWORK(N),int INFO)

    !print *, "INFO = ", info

    if (info /= 0) then
      print *, 'Diagonalizzazione Fallita. INFO = ', info
      print *, ""
    endif

  end subroutine diagUN


end module exponentiate
