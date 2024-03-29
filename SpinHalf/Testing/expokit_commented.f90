MODULE EXPOKIT
!
IMPLICIT NONE
!
CONTAINS

!----------------------------------------------------------------------|
      subroutine ZGPADM(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

      implicit none
      double precision t
      integer          ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      complex*16       H(ldh,m), wsp(lwsp)

!-----Purpose----------------------------------------------------------|
!*
!     Computes exp(t*H), the matrix exponential of a general complex 
!     matrix in full, using the irreducible rational Pade approximation
!     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
!     combined with scaling-and-squaring.
!*
!-----Arguments--------------------------------------------------------|
!*
!     ideg      : (input) the degre of the diagonal Pade to be used.
!                 a value of 6 is generally satisfactory.
!*
!     m         : (input) order of H.
!*
!     H(ldh,m)  : (input) argument matrix.
!*
!     t         : (input) time-scale (can be < 0).
!                  
!     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
!*
!     ipiv(m)   : (workspace)
!*
!*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
!                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
!                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!                 NOTE: if the routine was called with wsp(iptr), 
!                       then exp(tH) will start at wsp(iptr+iexph-1).
!*
!     ns        : (output) number of scaling-squaring used.
!*
!     iflag     : (output) exit flag.
!                       0 - no problem
!                      <0 - problem
!*
!----------------------------------------------------------------------|
!     Roger B. Sidje (rbs@maths.uq.edu.au)
!     EXPOKIT: Software Package for Computing Matrix Exponentials.
!     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
!----------------------------------------------------------------------|
!*
      integer i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
      double precision hnorm
      complex*16 cp, cq, scale, scale2, ZERO, ONE

      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX

!---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGPADM)'
!*
!---  initialise pointers ...
!*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
!*
!---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
!     and set scale = t/2^ns ...
!*
      do i = 1,m
         wsp(i) = ZERO
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
            !if (j == m) print *, real(wsp(i)), j !!!!!
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         !print*, hnorm, dble(wsp(i))  !!!!
         hnorm = MAX( hnorm,DBLE(wsp(i)) )
      enddo
      !print *, hnorm, t !!!!
      hnorm = ABS( t*hnorm )
      !print *, "hnorm = t*hnorm = ", hnorm !!!!!
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of ZGPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
      scale2 = scale*scale
!*
!---  compute Pade coefficients ...
!*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = ONE
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
!*
!---  H2 = scale2*H*H ...
!*
      call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
!*
!---  initialise p (numerator) and q (denominator) ...
!*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = ZERO
            wsp(iq + (j-1)*m + i-1) = ZERO
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
!*
!---  Apply Horner rule ...
!*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,&
                  wsp(ih2),m, ZERO,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
!*
!---  Obtain (+/-)(I + 2*(p\q)) ...
!*
      if ( iodd.ne.0 ) then
         call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m,&
                     H,ldh, ZERO,wsp(ifree),m )
         iq = ifree
      else
         call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,&
                     H,ldh, ZERO,wsp(ifree),m )
         ip = ifree
      endif
      call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
      call ZGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag.ne.0 ) stop 'Problem in ZGESV (within ZGPADM)'
      call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.ne.0 ) then
         call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
!*
!--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
!*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,&
                     ZERO,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
      END
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
      subroutine ZNCHBV( m, t, H,ldh, y, wsp )

      implicit none
      integer          m, ldh
      double precision t
      complex*16       H(ldh,m), y(m), wsp(m*(m+2))

!-----Purpose----------------------------------------------------------|
!*
!---  ZNCHBV computes y = exp(t*H)*y using the partial fraction
!     expansion of the uniform rational Chebyshev approximation
!     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg.
!     About 14-digit accuracy is expected if the matrix H is negative
!     definite. The algorithm may behave poorly otherwise. 
!*
!-----Arguments--------------------------------------------------------|
!*
!     m       : (input) order of the Hessenberg matrix H
!*
!     t       : (input) time-scaling factor (can be < 0).
!*
!     H(ldh,m): (input) upper Hessenberg matrix.
!*
!     y(m)    : (input/output) on input the operand vector,
!               on output the resulting vector exp(t*H)*y.
!*
!     wsp     : (workspace). Observe that a double precision vector of
!               length 2*m*(m+2) can be used as well when calling this
!               routine.
!*
!----------------------------------------------------------------------|
!     Roger B. Sidje (rbs@maths.uq.edu.au)
!     EXPOKIT: Software Package for Computing Matrix Exponentials.
!     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
!----------------------------------------------------------------------|
!*
      complex*16 ZERO
      integer ndeg, i, j, k, ip, ih, iy, iz
      parameter ( ndeg=7, ZERO=(0.0d0,0.0d0) )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg), tmpc

      intrinsic ABS,DBLE,CONJG,MIN
      
!---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

!---  Coefficients and poles of the partial fraction expansion...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
!*
      do ip = 1,ndeg
         theta(ndeg+ip) = CONJG( theta(ip) )
         alpha(ndeg+ip) = CONJG( alpha(ip) )
      enddo
!     
!---  Accumulation of the contribution of each pole ...
!*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,2*ndeg
         alpha(ip) = 0.5d0*alpha(ip)
!---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            wsp(iy+j-1) = wsp(iz+j-1)
            do i = 1,MIN( j+1,m )
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            do k = i,m
               wsp(ih+(j-1)*m+k-1) = ZERO
            enddo
         enddo
         do i = 1,m-1
!---        Get pivot and exchange rows ...
            if (ABS(wsp(ih+(i-1)*m+i-1)).lt.ABS(wsp(ih+(i-1)*m+i))) then
               call ZSWAP( m-i+1, wsp(ih+(i-1)*m+i-1),m, &
                          wsp(ih+(i-1)*m+i),m )
               call ZSWAP( 1, wsp(iy+i-1),1, wsp(iy+i),1 )
            endif
!---        Forward eliminiation ... 
            tmpc = wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1)
            call ZAXPY( m-i, -tmpc, wsp(ih+i*m+i-1),m, wsp(ih+i*m+i),m )
            wsp(iy+i) = wsp(iy+i) - tmpc*wsp(iy+i-1)
         enddo
!---     Backward substitution ...    
         do i = m,1,-1
            tmpc = wsp(iy+i-1)
            do j = i+1,m
               tmpc = tmpc - wsp(ih+(j-1)*m+i-1)*wsp(iy+j-1)
            enddo
            wsp(iy+i-1) = tmpc / wsp(ih+(i-1)*m+i-1)
         enddo
!---     Accumulate the partial result in y ...     
         do j = 1,m
            y(j) = y(j) + alpha(ip)*wsp(iy+j-1)
         enddo
      enddo
      END
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
      subroutine ZGEXPV( n, a, ia, ja, nz, m, t, v, w, tol, anorm,&
                        wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

      implicit none
      integer          n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp), nz
      double precision t, tol, anorm
      integer   ia(nz), ja(nz)
      complex*16       v(n), w(n), wsp(lwsp), a(nz)
      external         matvec

!-----Purpose----------------------------------------------------------|
!*
!---  ZGEXPV computes w = exp(t*A)*v
!     for a Zomplex (i.e., complex double precision) matrix A 
!*
!     It does not compute the matrix exponential in isolation but
!     instead, it computes directly the action of the exponential
!     operator on the operand vector. This way of doing so allows 
!     for addressing large sparse problems. 
!*
!     The method used is based on Krylov subspace projection
!     techniques and the matrix under consideration interacts only
!     via the external routine `matvec' performing the matrix-vector 
!     product (matrix-free method).
!*
!-----Arguments--------------------------------------------------------|
!*
!     n      : (input) order of the principal matrix A.
!*
!     m      : (input) maximum size for the Krylov basis.
!*
!     t      : (input) time at wich the solution is needed (can be < 0).
!*
!     v(n)   : (input) given operand vector.
!*
!     w(n)   : (output) computed approximation of exp(t*A)*v.
!*
!     tol    : (input/output) the requested accuracy tolerance on w. 
!              If on input tol=0.0d0 or tol is too small (tol.le.eps)
!              the internal value sqrt(eps) is used, and tol is set to
!              sqrt(eps) on output (`eps' denotes the machine epsilon).
!              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
!*
!     anorm  : (input) an approximation of some norm of A.
!*
!   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
!                                   +---------+-------+---------------+
!              (actually, ideg=6)        V        H      wsp for PADE
!                   
! iwsp(liwsp): (workspace) liwsp .ge. m+2
!*
!     matvec : external subroutine for matrix-vector multiplication.
!              synopsis: matvec( x, y )
!                        complex*16 x(*), y(*)
!              computes: y(1:n) <- A*x(1:n)
!                        where A is the principal matrix.
!*
!     itrace : (input) running mode. 0=silent, 1=print step-by-step info
!*
!     iflag  : (output) exit flag.
!              <0 - bad input arguments 
!               0 - no problem
!               1 - maximum number of steps reached without convergence
!               2 - requested tolerance was too high
!*
!-----Accounts on the computation--------------------------------------|
!     Upon exit, an interested user may retrieve accounts on the 
!     computations. They are located in the workspace arrays wsp and 
!     iwsp as indicated below: 
!*
!     location  mnemonic                 description
!     -----------------------------------------------------------------|
!     iwsp(1) = nmult, number of matrix-vector multiplications used
!     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
!     iwsp(3) = nscale, number of repeated squaring involved in Pade
!     iwsp(4) = nstep, number of integration steps used up to completion 
!     iwsp(5) = nreject, number of rejected step-sizes
!     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
!     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
!     -----------------------------------------------------------------|
!     wsp(1)  = step_min, minimum step-size used during integration
!     wsp(2)  = step_max, maximum step-size used during integration
!     wsp(3)  = x_round, maximum among all roundoff errors (lower bound) 
!     wsp(4)  = s_round, sum of roundoff errors (lower bound)
!     wsp(5)  = x_error, maximum among all local truncation errors
!     wsp(6)  = s_error, global sum of local truncation errors
!     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
!     wsp(8)  = t_now, integration domain successfully covered
!     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
!     wsp(10) = ||w||/||v||, scaled norm of the solution w.
!     -----------------------------------------------------------------|
!     The `hump' is a measure of the conditioning of the problem. The
!     matrix exponential is well-conditioned if hump = 1, whereas it is
!     poorly-conditioned if hump >> 1. However the solution can still be
!     relatively fairly accurate even when the hump is large (the hump 
!     is an upper bound), especially when the hump and the scaled norm
!     of w [this is also computed and returned in wsp(10)] are of the 
!     same order of magnitude (further details in reference below).
!*
!----------------------------------------------------------------------|
!-----The following parameters may also be adjusted herein-------------|
!*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,&
                mxreject = 0,&
                ideg     = 6,&
                delta    = 1.2d0,&
                gamma    = 0.9d0 )

!     mxstep  : maximum allowable number of integration steps.
!               The value 0 means an infinite number of steps.
! 
!     mxreject: maximum allowable number of rejections at each step. 
!               The value 0 means an infinite number of rejections.
!*
!     ideg    : the Pade approximation of type (ideg,ideg) is used as 
!               an approximation to exp(H). The value 0 switches to the
!               uniform rational Chebyshev approximation of type (14,14)
!*
!     delta   : local truncation error `safety factor'
!*
!     gamma   : stepsize `shrinking factor'
!*
!----------------------------------------------------------------------|
!     Roger B. Sidje (rbs@maths.uq.edu.au)
!     EXPOKIT: Software Package for Computing Matrix Exponentials.
!     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
!----------------------------------------------------------------------|

      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,&
             ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,&
             nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,&
                      s_error, x_error, t_now, t_new, t_step, t_old,&
                      xm, beta, break_tol, p1, p2, p3, eps, rndoff,&
                      vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hij

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
!*
!---  check restrictions on input parameters ...
!*
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGEXPV)'
!*
!---  initialisations ...
!*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      !print *, "t_out (initial) = ", t_out !!!!!
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
!*>>>  break_tol = tol
!*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call ZCOPY( n, v,1, w,1 )
      beta = DZNRM2( n, w,1 )
      vnorm = beta
      hump = beta 
!*
!---  obtain the very first stepsize ...
!*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      !print *, "anorm (t_new) = ", anorm !!!!!
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      !print *, "t_new (1) = ", t_new !!!!!
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) ! p1
      !print *, "t_new (2) = ", t_new !!!!!
!*
!---  step-by-step integration ...
!*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      !print *, "t_step (t_out>t_now) = ", t_step !!!!!
      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo
!*
!---  Arnoldi loop ...
!*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) ,a, ia, ja, nz,n)
         do i = 1,j
            hij = ZDOTC( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call ZAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DZNRM2( n, wsp(j1v),1 )
!---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            !print *, "t_step (happy breakdown) = ", t_step !!!!
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v), a, ia, ja, nz, n)
      avnorm = DZNRM2( n, wsp(j1v),1 )
!*
!---  set 1 for the 2-corrected scheme ...
!*
 300  continue
      wsp(ih+m*mh+m+1) = ONE
!*
!---  loop while ireject<mxreject until the tolerance is reached ...
!*
      ireject = 0
 401  continue
!*
!---  compute w = beta*V*exp(t_step*H)*e1 ...
!*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
!---     irreducible rational Pade approximation ...
         !print *, "Call to ZGPADM"
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,&
                     wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
!---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         enddo
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif
 402  continue
! 
!---  error estimate ...
! 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   ! beta
         p2 = ABS( wsp(iexph+m+1) ) ! beta ! avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
!*
!---  reject the step-size if the error is not acceptable ...
!   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.&
          (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma ! t_step ! (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         print *, "p1 = ", p1 !!!!!
         t_step = AINT( t_step/p1 + 0.55d0 ) ! p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZGEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
!*
!---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
!*
      mx = mbrkdwn + MAX( 0,k1-1 )
      hij = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hij,wsp(iv),n,wsp(iexph),1,ZERO,w,1 )
      beta = DZNRM2( n, w,1 )
      hump = MAX( hump, beta )
!*
!---  suggested value for the next stepsize ...
!*
      t_new = gamma ! t_step ! (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      !print *, " p1 = ", p1 !!!!
      t_new = AINT( t_new/p1 + 0.55d0 ) ! p1

      err_loc = MAX( err_loc,rndoff )
!*
!---  update the time covered ...
!*
      t_now = t_now + t_step
!*
!---  display and keep some information ...
!*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      wsp(9)  = CMPLX( hump/vnorm )
      wsp(10) = CMPLX( beta/vnorm )
      END
!----------------------------------------------------------------------|
END MODULE EXPOKIT
