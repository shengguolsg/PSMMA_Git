module BasicMM
  implicit none

contains

!!!!!!!!
  SUBROUTINE AddOneRow( LDA, K, N, A, V )
!
! ..Scalar Arguments..
    INTEGER, INTENT(IN) :: LDA, K, N
! ..
! ..Array Arguments..
    DOUBLE PRECISION, INTENT(INOUT) :: A(LDA,*), V(*)
!
! Purpose
! =======
! Update the K-th row of A by adding the vector V, 
! A(k,1:N)= A(k,1:N) + V(1:N).
! 
! =========
! Written by S.-G. Li, on May, 2nd, 2013
! ======================================
!
    INTEGER I
    
    DO I = 1, N
       A( k, I ) = A( k,I ) + V( I )
    END DO

  END SUBROUTINE AddOneRow

!!!!!!!!
  SUBROUTINE Cauchy(A,X,Y,Z,N)
! 
! Construction an Cauchy matrix A_{ij} = z_j / (x_j - y_i)
! 
! X column basis
! Y row basis
! Z nominator
!
    REAL(8)  :: A(N,*), X(*), Y(*), Z(*)
    INTEGER  :: N       ! Leading dim of A
!
    INTEGER  :: I,J
    
    DO J = 1, N
       DO I = 1, N
          A(I,J) = Z(J) / ( X(J)-Y(I) )
       END DO
    END DO
    
  END SUBROUTINE Cauchy  

!!!!!!!!
  subroutine Cauchylike(A,D,F,u,v,M,N)
! A(i,j) = u(i)v(j) / (D(i)-F(j))
    integer M,N,i,j
    double precision A(M,*), D(*),F(*),u(*),v(*)

    do j = 1,N
       do i = 1,M
          A(i,j) = u(i)*v(j) / (D(i)-F(j))
       end do
    end do
  end subroutine Cauchylike

!!!!!!!
  SUBROUTINE Cauchy_VT1(A,D,F,U,Z,DIFL,DIFR,M,N)
! 
! .. Scalar Arguments ..
    INTEGER  M, N
! .. Array Arguments ..
    DOUBLE PRECISION A(M,*),D(*),F(*),U(*),Z(*),DIFL(*),DIFR(*)
!
! Purpose
! =======
! It computes the right singular vector matrix, VT, for updating SVD problem. 
! On entry, VT( i,j ) = U(i)*V(j) / ( F(j)**2 - D(i)**2 ). This routine returns 
! the same matrix as that by Cauchy_VT. The difference is that 
!  1) D is given in this routine; Cauchy_VT needs to compute it by calling dlasd4 properly;
!  2) A is computed by using vector operations here; Cauchy_VT uses matrix-matrix multiplication,
!     it needs more memory;
!  3) Before using this routine, mdlasd31.f90 must be called first.
!
!  It is reasonable to choose M == N, though it also works when M not equal to N.
! 
! Parameters
! ==========
! A  (output)  DOUBLE PRECISION array, DIMENSION ( M, N )
!    It is the transpose of right singular vector matrix for updating SVD problem.
! 
! D  (input) DOUBLE PRECISION array, dimension ( N )
!    On entry, D are the new singular values of B in ascending order, B = [B0; Z^T] where
!    Z is an appended vector. F(i) < D(i) < F(i+1). The singular values of B0 are F.
!
! F  (input) DOUBLE PRECISION array, dimension ( M )
!    The singular values of B0 in ascending order. These are the poles of the secular equation.

! U  (input) DOUBLE PRECISION array, dimension ( M )
!    The row normalization scalar of matrix A, or the transpose of right singular vector matrix, VT.
!    
! Z  (input) DOUBLE PRECISION array, dimension ( N )
!    The appended row vector, the entries of which have been modified by mdlasd31.f90.
!
! DIFL (input) DOUBLE PRECISION array, dimension ( MAX(M,N) )
!      The difference between old singular values and new singular values
!      DIFL(J) = D(J) - F(J), positive value, the same as LAPACK. 
!
! DIFR (INPUT) DOUBLE PRECISION array, DIMENSION( MAX(M,N) )
!      The difference between old singular values and new singular values
!      DIFR(J) = D(J) - F(J+1), negative value, the same as LAPACK. DIFR( MAX(M,N) )
!      will not be referenced. 
!
! M   (input) INTEGER, row dimension of A
!
! N   (input) INTEGER, column dimension of A. 
!     N should equal to M. 
!      
! ============
! Written by S.-G. Li, in Changsha China, on Dec. 14th, 2012
! =========================================================
!
! .. Local Scalars ..
    INTEGER   I,J
    DOUBLE PRECISION  DIFLJ,ZJ,FJ
! ..
!     .. External Functions ..
    double precision   DLAMC3
    external           DLAMC3

    DO J = 1, N
         DIFLJ = DIFL( J )
         ZJ = Z(J)
         FJ = F( J )
         A(J, J ) = -U(J)*ZJ / DIFLJ / ( D( J )+FJ ) 
         DO I = 1, J - 1
            A( I,J ) = U(I)*ZJ / ( DLAMC3( FJ, -F(I+1) )-DIFR(I) ) / ( D( I )+FJ ) 
         END DO
         DO I = J + 1, M
            A( I,J ) = U(I)*ZJ / ( DLAMC3( FJ, -F(I) )-DIFL(I) ) / ( D( I )+FJ ) 
         END DO
      END DO

    END SUBROUTINE Cauchy_VT1

!!!!!!!
  subroutine CauchylikesvdU(A,D,F,U,V,DIFL,DIFR,M,N)
! A(i,j) = U(i)V(j) / (D(i)^2-F(j)^2)
! only works for M == N
! U is the normalizing scalars and V is Z
! D is the updated singular values and F is the old ones.

    integer M,N,i,j
    double precision DIFLJ,VJ,FJ
    double precision A(M,*),D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)

!     .. External Functions ..
    double precision   DLAMC3
    external           DLAMC3


    DO J = 1, N-1
         DIFLJ = DIFL( J )
         VJ = V(J)
         FJ = F( J )
         A(J, J ) = U(J)*VJ / DIFLJ / ( D( J )+FJ ) 
         DO I = 1, J - 1
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I+1),-FJ )+DIFR(I) ) / ( D( I )+FJ ) 
         END DO
         DO I = J + 1, M
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I),-FJ )+DIFL(I) ) / ( D( I )+FJ ) 
         END DO
      END DO

      DO I = 1, M
         A(I,N) = U(I)
      END DO

    end subroutine CauchylikesvdU

!!!!!!!
  subroutine CauchylikesvdU1(A,D,F,U,V,DIFL,DIFR,M,N)
! A(i,j) = U(i)V(j) / (D(i)^2-F(j)^2)
! only works for M == N
! U is the normalizing scalars and V is Z
! D is the updated singular values and F is the old ones.

    integer M,N,i,j
    integer JX(M)
    double precision DIFLJ,VJ,FJ
    double precision A(M,*),D(*),F(*),U(*),V(*),DIFL(*),DIFR(*)

!     .. External Functions ..
    double precision   DLAMC3
    external           DLAMC3

    DO J = 1, N-1
         DIFLJ = DIFL( J )
         VJ = V(J)
         FJ = F( J )
         A(J, J ) = U(J)*VJ / DIFLJ / ( D( J )+FJ ) 
         DO I = 1, J - 1
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I+1),-FJ )+DIFR(I) ) / ( D( I )+FJ ) 
         END DO
         DO I = J + 1, M
            A( I,J ) = U(I)*VJ / ( DLAMC3( F(I),-FJ )+DIFL(I) ) / ( D( I )+FJ ) 
         END DO
      END DO

      DO I = 1, M
         A(I,N) = U(I)
      END DO

!  reorder matrix A
      JX(1:M) = (/ (i, i=M,1,-1) /)
      A(1:M,1:N)= A(JX,1:N)
      A(1:M,1:M)=A(1:M,JX)

    end subroutine CauchylikesvdU1

!!!!!!!
    subroutine testOrthLeft(D,D1,ALPHA_L,ZZ,DIFL,DIFR,N)
!
! .. Scalar Arguments ..
      integer N
! .. Array Arguments ..
      double precision :: D(*),D1(*),ALPHA_L(*),ZZ(*),DIFL(*),DIFR(*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of left singular vector matrix.
! 
! .. Parameters ..
! D       (input)  double precision, the old singular values
! D1      (input)  double precision, the updated singular values
! ALPHA_L (input)  double precision, the scalars for the left singular vector matrix
! ZZ      (input)  double precision, the vector of the numerator
! DIFL    (input)  double precision
! DIFR    (input)  double precision

      integer lwork,info
      double precision err
      double precision, allocatable :: A(:,:),work(:),Dt(:)

      lwork = N*(N+1)
      allocate( A(N,N),work(lwork),Dt(N) )

!      call CauchylikesvdU(A,D1,D,ALPHA_L,ZZ,DIFL,DIFR,N,N+1) 
      call Cauchy_VT1(A,D1,D,ALPHA_L,ZZ,DIFL,DIFR,N,N) 

      ! *****************************************************
      !           Check the orthogonality of V              *
      ! *****************************************************
      Dt = 0.0D0
      call dgesvd('N','N',N,N,A,N,Dt,A,N,A,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
      write(*,*) "Othogonality of exact: ", err

      deallocate( A,work,Dt )

    end subroutine testOrthLeft

!!!!!!!
    subroutine testOrthLeft1(D,D1,ALPHA_L,ZZ,DIFL,DIFR,N)
!
! .. Scalar Arguments ..
      integer N
! .. Array Arguments ..
      double precision :: D(*),D1(*),ALPHA_L(*),ZZ(*),DIFL(*),DIFR(*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of left singular vector matrix.
! 
! .. Parameters ..
! D       (input)  double precision, the old singular values
! D1      (input)  double precision, the updated singular values
! ALPHA_L (input)  double precision, the scalars for the left singular vector matrix
! ZZ      (input)  double precision, the vector of the numerator
! DIFL    (input)  double precision
! DIFR    (input)  double precision

      integer lwork,info
      double precision err
      double precision, allocatable :: A(:,:),work(:),Dt(:)

      lwork = N*(N+1)
      allocate( A(N,N+1),work(lwork),Dt(N) )

      call CauchylikesvdU1(A,D1,D,ALPHA_L,ZZ,DIFL,DIFR,N,N+1) 

      ! *****************************************************
      !           Check the orthogonality of V              *
      ! *****************************************************
      Dt = 0.0D0
      call dgesvd('N','N',N,N+1,A,N,Dt,A,N+1,A,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      !     write(*,*) 'The svals of A are ', Dt(1:n)
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
      write(*,*) "Othogonality of exact: ", err

      deallocate( A,work,Dt )

    end subroutine testOrthLeft1

!!!!!!!
    SUBROUTINE TestOrth(V,M,N)
!
! .. Scalar Arguments ..
      integer M,N
! .. Array Arguments ..
      double precision :: V(M,*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of matrix V.
! 
! .. Local Scalars ..
      integer lwork,info,MN, ierr
      double precision err
! .. Local Arrays ..
      double precision, allocatable :: A(:,:),work(:),Dt(:)

      lwork = N*(N+1)
      IF( M < N) THEN
         MN = M
      ELSE
         MN = N
      END IF
      ALLOCATE( A(M,N),work(lwork),Dt(MN), stat=ierr )
      IF( ierr .ne. 0 ) THEN
         WRITE(*,*) 'Allocate failed in TestOrth '
         RETURN
      END IF      

! *********************************************************
!               Check the orthogonality of V              *
! *********************************************************
      Dt = 0.0D0
      A(1:M,1:N) = V(1:M,1:N)
      call dgesvd('N','N',M,N,A,M,Dt,A,M,A,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(MN)-1.0D0) )
      WRITE(*,*) "Test Othogonality of matrix V ", err

!      Write(*,*) "The svals are ", Dt(1:MN)

      DEALLOCATE( A,work,Dt )

    END SUBROUTINE TestOrth

!!!!!!!
  SUBROUTINE testlowrank( M,N,A,LDA )

        IMPLICIT NONE
!    
      INTEGER :: M,N,LDA
      DOUBLE PRECISION :: A(*)
!           
!    This routine compute the numerical rank of matrix A sequentially. 
!    
      INTEGER :: prank, LWORK, INFO
      DOUBLE PRECISION, PARAMETER :: ptol = 1.0E-13
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), SV(:), TempB(:)

      LWORK = 10*LDA
      ALLOCATE( WORK(LWORK), SV(LDA), TempB(LDA*N) )

      CALL DLACPY( 'A',M,N,A,LDA,TempB,LDA)
!      WRITE(*,*) "M,N,LDA=",M,N,LDA

      SV = 0.0D+0
      CALL DGESVD( 'N','N', M, N, TempB, LDA,SV, &
                  WORK,LDA,WORK,LDA,WORK,LWORK,INFO)
      prank = count(SV > ptol)
      WRITE(*,*) "prank=", prank,M, N

      DEALLOCATE( SV )
      DEALLOCATE( TempB )
      DEALLOCATE( WORK )

      END SUBROUTINE testlowrank

!!!!!!!
  SUBROUTINE Cauchylike_md2(A,D,W,U,V,M,N)
! 
! .. Scalar Arguments ..
    INTEGER   M, N
! .. Array Arguments ..
    DOUBLE PRECISION  A(M,*),D(*),W(*),U(*),V(*) 
!
! Purpose 
! ========
! It returns a Cauchy matrix A(i,j) = U(i)*V(j) / (D(j)**2 - W(i)**2 )
! A is an M-by-N matrix. 
! 

    INTEGER  I,J 
    DO J = 1, N
         DO I = 1, M
            A( I,J ) = U(I)*V(J) / ( (D(J)+W(I)) * (D(J)-W(I)) ) 
         END DO
      END DO

  END SUBROUTINE Cauchylike_md2

!!!!!!!
  SUBROUTINE Cauchy_VT( N, D, DSIGMA, VT, Z, INFO )
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D(*), DSIGMA(*), VT(N, *), Z(*)
!     ..
!
!  Purpose
!  =======
!
!  Cauchy_VT finds all the singular values of M =[DSIGMA; Z].  It makes the
!  appropriate calls to DLASD4. It returns the right singular matrix VT.
!  The singular values of M are stored in D. Z will be modified by Lowner Theorem. 
!  DSIGMA will also be changed slightly. 
!
!  Cauchy_VT can be used for testing the correctness of Cauchy2hssvd. 
!  The form of VT is VT( i,j ) = alpha(i)*Z(j) / ( DSIGMA(j)**2 - D(i)**2 ), where
!  alpha(i) is the normalization scalar of each row of VT. 
!
!  VT must be a square matrix. 
!
!  Arguments
!  =========
!
!  N     (input) INTEGER
!         The column dimension of the matrix M, or the dimension of matrix 
!         A = DIAG( DSIGMA )* DIAG( DSIGMA ) + Z^T * Z. 
!
!  D      (output) DOUBLE PRECISION array, dimension(N)
!         On exit the square roots of the roots of the secular equation
!         in ascending order, the singular values of M.
!
!  DSIGMA (input) DOUBLE PRECISION array, dimension(N)
!         The singular values of previous matrix in ascending order. 
!         These are the poles of the secular equation.
!
!  VT     (output) DOUBLE PRECISION array, dimension (N, N)
!         It contains the transpose of right singular matrix of M, the singular values of 
!         which are in ascending order. 
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         It contains the appended row vector and it will be updated after this routine. 
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!   Modified from mdlasd3 in LAPACK by Shengguo Li, on Dec. 13th, 2012.
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, K, ierr
      DOUBLE PRECISION   RHO, TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION, ALLOCATABLE :: U(:,:), Q(:)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLASCL, DLASD4
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      K  = N
! 
!    U    (workspace) DOUBLE PRECISION array, dimension (N,N)
!         Used to compute VT. 
!
!    Q    (workspace) DOUBLE PRECISION array, dimension (N)
!         Used to copy Z, let the updated Z have the same sign as Z. 
      ALLOCATE( U(N,N), Q(N), stat = ierr )
      IF( ierr .ne. 0 ) THEN
         WRITE(*,*) 'Allocate failed in Cauchy_VT '
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2!DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 20 I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   20 CONTINUE
!
!     Keep a copy of Z.
!
      CALL DCOPY( K, Z, 1, Q, 1 )
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Find the new singular values.
!
      DO 30 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ), &
                     VT( 1, J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
   30 CONTINUE
!
!     Compute updated Z.
!
      DO 60 I = 1, K
         Z( I ) = U( I, K )*VT( I, K )
         DO 40 J = 1, I - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J ) ) / &
                    ( DSIGMA( I )+DSIGMA( J ) ) )
   40    CONTINUE
         DO 50 J = I, K - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J+1 ) ) /  &
                    ( DSIGMA( I )+DSIGMA( J+1 ) ) )
   50    CONTINUE
         Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I ) )
   60 CONTINUE
!
!     Generate the right singular vectors V
!
      DO 90 I = 1, K
         U( 1, I ) = Z( 1 ) / U( 1, I ) / VT( 1, I )
         DO 70 J = 2, K
            U( J, I ) = Z( J ) / U( J, I ) / VT( J, I )
   70    CONTINUE
         TEMP = DNRM2( K, U( 1, I ), 1 )
         TEMP = ONE / TEMP
         DO 80 J = 1, K
            U( J, I ) = U( J, I ) * TEMP
   80    CONTINUE
   90 CONTINUE
!       
!   Generate VT
!
      DO  J = 1, K
         DO I = 1, K
            VT( I,J ) = U( J,I )
         END DO
      END DO
!
    deallocate( U, Q )
!
      RETURN
!
!     End of Cauchy_VT
!
  END SUBROUTINE CAUCHY_VT

!!!!!!! 
  SUBROUTINE Cauchy_SVD1(N, D, DSIGMA, VT, Z, INFO )
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D(*), DSIGMA(*), VT(N, *), Z(*)
!     ..
!
!  Purpose
!  =======
!
!  Cauchy_SVD1 computes the singular values of M =[ DIAG(DSIGMA); Z] and the transpose of its right
!  singular vector matrix, VT. The difference with CAUCHY_VT is that now the computed singular values
!  D are in decreasing order. 
!
!  Cauchy_SVD1 is called from updating SVD test routines. The correctness needs to be checked later. 
!
!  Arguments
!  =========
!
!  N     (input) INTEGER
!         The column dimension of the matrix M, or the dimension of matrix 
!         A = DIAG( DSIGMA )* DIAG( DSIGMA ) + Z^T * Z. 
!
!  D      (output) DOUBLE PRECISION array, dimension(N)
!         On exit the square roots of the roots of the secular equation
!         in ascending order, the singular values of M.
!
!  DSIGMA (input) DOUBLE PRECISION array, dimension(N)
!         The singular values of previous matrix in ascending order. 
!         These are the poles of the secular equation.
!
!  VT     (output) DOUBLE PRECISION array, dimension (N, N)
!         It contains the transpose of right singular matrix of M, the singular values of 
!         which are in descending order. 
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         It contains the appended row vector and it will be updated after this routine. 
!
!  INFO   (output) INTEGER
!         = 0:  successful exit.
!         < 0:  if INFO = -i, the i-th argument had an illegal value.
!         > 0:  if INFO = 1, a singular value did not converge
!
!  Further Details
!  ===============
!   Modified from mdlasd3 in LAPACK by Shengguo Li, on Aug. 29th, 2012.
!
!   Modified on Dec. 13th, 2012.
!   It should be correct. But its correctness and usage are needed to 
!   be checked later. 
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, K, ierr
      DOUBLE PRECISION   RHO, TEMP
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION, ALLOCATABLE :: U(:,:), Q(:)
      INTEGER    JX(N)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLASCL, DLASD4
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      K  = N
      ALLOCATE( U(N,N), Q(N), stat=ierr )
      IF( ierr .ne. 0 ) THEN
         WRITE(*,*) 'Allocate failed in Cauchy_SVD1'
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
      END IF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2!DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 20 I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   20 CONTINUE
!
!     Keep a copy of Z.
!
      CALL DCOPY( K, Z, 1, Q, 1 )
!
!     Normalize Z.
!
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!
!     Find the new singular values.
!
      DO 30 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ), &
                     VT( 1, J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
   30 CONTINUE
!
!     Compute updated Z.
!
      DO 60 I = 1, K
         Z( I ) = U( I, K )*VT( I, K )
         DO 40 J = 1, I - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J ) ) / &
                    ( DSIGMA( I )+DSIGMA( J ) ) )
   40    CONTINUE
         DO 50 J = I, K - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / &
                    ( DSIGMA( I )-DSIGMA( J+1 ) ) /  &
                    ( DSIGMA( I )+DSIGMA( J+1 ) ) )
   50    CONTINUE
         Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I ) )
   60 CONTINUE
!
!     Generate the right singular vectors V
!
      DO 90 I = 1, K
         VT( 1, I ) = Z( 1 ) / U( 1, I ) / VT( 1, I )
         DO 70 J = 2, K
            VT( J, I ) = Z( J ) / U( J, I ) / VT( J, I )
   70    CONTINUE
         TEMP = DNRM2( K, VT( 1, I ), 1 )
         TEMP = ONE / TEMP
         DO 80 J = 1, K
            VT( J, I ) = VT( J, I ) * TEMP
   80    CONTINUE
   90 CONTINUE
!
!     Generate the right singular vectors VT, which must try 'transpose' later
!
!        VT(1:K,1:K) = transpose(Vt(1:K,1:K))

      DO 120 I = 1, K
         DO 110 J = 1, K
            U( I, J ) = VT( J, I )
  110    CONTINUE
  120 CONTINUE
            
      DO I = 1,N
         JX(I) = N-I+1
      END DO

      DO I = 1, N
         DO J = 1,N
            VT(JX(J),JX(I)) = U(J,I)
        END DO
      END DO

    DEALLOCATE( U, Q )
!
      RETURN
!
!     End of Cauchy_SVD1
!
  END SUBROUTINE CAUCHY_SVD1

!!!!!!!
    SUBROUTINE TCauchyDC( A,D,F,Z,APHAV,DIFL,DIFR,M,N, WORK )
!
! .. Scalar Arguments ..
    INTEGER M, N
! .. Array Arguments ..
    DOUBLE PRECISION A(M,*),D(*),F(*),Z(*),APHAV(*),DIFL(*),DIFR(*), &
                     WORK(*)
!
! Purpose
! =======
! It tests whether the Cauchy matrix defined by D, F, APHAV and Z is 
! orthogonal. This routine is designed for Structured DC algorithm.
! 

! .. Local Scalars
    INTEGER  LWORK, MN, info, ierr
    DOUBLE PRECISION  err, one
! ..
! .. Local Arrays 
    DOUBLE PRECISION, ALLOCATABLE :: Dt(:)

    MN = MAX( M, N )
    lwork = 5*MN
    ALLOCATE( Dt(MN), stat=ierr )
    IF( ierr .ne. 0 ) THEN
       WRITE(*,*) 'Allocate failed in TCauchyDC.'
    END IF

    call Cauchy_VT1( A,D,F,APHAV,Z,DIFL,DIFR,M,N )
    call dgesvd('N','N',N,N,A,N,Dt,A,M,A,N,WORK,lwork,info )
    
    one = 1.0E+0
    IF( info .eq. 0 ) THEN
       err = MAX( ABS( Dt(1)-one), ABS(Dt(MN)- one) )
       WRITE(*,*) 'Cauchy orthogonality is ', err
    ELSE
       WRITE(*,*) 'Error in Cauchy'
    END IF
    
    DEALLOCATE( Dt )

  END SUBROUTINE TCauchyDC

!!!!!!!
  SUBROUTINE CauchyDC_VT(A,D,F,U,Z,DIFL,DIFR,M,N)
! 
! .. Scalar Arguments ..
    INTEGER  M, N
! .. Array Arguments ..
    DOUBLE PRECISION A(M,*),D(*),F(*),U(*),Z(*),DIFL(*),DIFR(*)
!
! Purpose
! =======
! It computes the right singular vector matrix, VT, for bidiagonal DC algorithm.
! On entry, VT( i,j ) = U(i)*Z(j) / ( F(j)**2 - D(i)**2 ). 
! VT is the right singular vector matrix of M, M = [Z; diag(F)], 
! where F(1) = 0, and F(i)< F(i+1). 
!
! Parameters
! ==========
! A  (output)  DOUBLE PRECISION array, DIMENSION ( M, N )
!    It is the transpose of right singular vector matrix of matrix M. 
! 
! D  (input) DOUBLE PRECISION array, dimension ( N )
!    On entry, D are the new singular values of M in ascending order, M = [Z^T; diag(F)] where
!    Z is an appended vector. F(i) < D(i) < F(i+1). 
!
! F  (input) DOUBLE PRECISION array, dimension ( M )
!    The old singular values in ascending order. These are the poles of the secular equation.

! U  (input) DOUBLE PRECISION array, dimension ( M )
!    The row normalization scalar of matrix A, or the transpose of right singular vector matrix, VT.
!    
! Z  (input) DOUBLE PRECISION array, dimension ( N )
!    The appended row vector, or the column generator of A;
!
! DIFL (input) DOUBLE PRECISION array, dimension ( MAX(M,N) )
!      The difference between old singular values and new singular values
!      DIFL(J) = D(J) - F(J), positive value, the same as LAPACK. 
!
! DIFR (INPUT) DOUBLE PRECISION array, DIMENSION( MAX(M,N) )
!      The difference between old singular values and new singular values
!      DIFR(J) = D(J) - F(J+1), negative value, the same as LAPACK. DIFR( MAX(M,N) )
!      will not be referenced. 
!
! M   (input) INTEGER, row dimension of A
!
! N   (input) INTEGER, column dimension of A. 
!     N should equal to M. 
!      
! ============
! Written by S.-G. Li, in Changsha China, on Mar. 11th, 2013
! =========================================================
!
! .. Local Scalars ..
    INTEGER   I,J
    DOUBLE PRECISION  DIFLJ,ZJ,FJ
! ..
!     .. External Functions ..
    double precision   DLAMC3
    external           DLAMC3

    DO J = 1, N
         DIFLJ = DIFL( J )
         ZJ = Z(J)
         FJ = F( J )
         A(J, J ) = -U(J)*ZJ / DIFLJ / ( D( J )+FJ ) 
         DO I = 1, J - 1
            IF(I.EQ. 1120 .and. J .eq. 1121 ) Then
!               write(*,*) 'Look check it' 
            END IF
            A( I,J ) = U(I)*ZJ / ( DLAMC3( FJ, -F(I+1) )-DIFR(I) ) / ( D( I )+FJ ) 
         END DO
         DO I = J + 1, M
            A( I,J ) = U(I)*ZJ / ( DLAMC3( FJ, -F(I) )-DIFL(I) ) / ( D( I )+FJ ) 
         END DO
      END DO

    END SUBROUTINE CauchyDC_VT

!!!!!!!
    subroutine testOrthRight(A, D,F,ALPHA_L,Z,DIFL,DIFR,N)
!
! .. Scalar Arguments ..
      integer N
! .. Array Arguments ..
      double precision :: D(*),F(*),ALPHA_L(*),Z(*),DIFL(*),DIFR(*),A(N,*)
!
! .. Purpose ..
! =============      
! This routine is written for testing the orthogonality of left singular vector matrix.
! 
! .. Parameters ..
! F       (input)  double precision
!         the old singular values
! D       (input)  double precision
!         the updated singular values
! ALPHA_L (input)  double precision
!         row scalars of the transpose of right singular vector matrix
!
! Z       (input)  double precision, the vector of the numerator
!
! DIFL (input) DOUBLE PRECISION array, dimension ( MAX(M,N) )
!      The difference between old singular values and new singular values
!      DIFL(J) = D(J) - F(J), positive value, the same as LAPACK. 
!
! DIFR (INPUT) DOUBLE PRECISION array, DIMENSION( MAX(M,N) )
!      The difference between old singular values and new singular values
!      DIFR(J) = D(J) - F(J+1), negative value, the same as LAPACK. DIFR( MAX(M,N) )
!      will not be referenced. 

      integer lwork,info
      double precision err
      double precision, allocatable :: work(:),Dt(:), V(:,:)

      lwork = N*(N+1)
      allocate( V(N,N),work(lwork),Dt(N) )

!      call CauchylikesvdU(A,D1,D,ALPHA_L,ZZ,DIFL,DIFR,N,N+1) 
      call CauchyDc_VT(A,D,F,ALPHA_L,Z,DIFL,DIFR,N,N) 

      ! *****************************************************
      !           Check the orthogonality of V              *
      ! *****************************************************
      Dt = 0.0D0
      V(1:N,1:N) = A(1:N,1:N)
      call dgesvd('N','N',N,N,V,N,Dt,V,N,V,N,work,lwork,info)  !svd of V1  V1 would be destroyed
      err = max( ABS(Dt(1)-1.0D0), ABS(Dt(N)-1.0D0) )
      write(*,*) "Othogonality of exact: ", err

      deallocate( V,work,Dt )

    end subroutine testOrthRight

!!!!!!!!
    SUBROUTINE TzStruct( A, LDA, M, N ) 
! 
!   Test the structure of matrix A
!
      INTEGER  M, N, LDA, num
      DOUBLE PRECISION A( LDA,* )
!
      INTEGER i, j
      DOUBLE PRECISION  TOL

      TOL  = 1.0E-13
      num = 0
      write(*,*) ' Test structure '
      DO j = 1, N
         DO i = 1, M
            IF( abs(A(i,j)) > TOL ) then
               num = num + 1
               write(*,*) 'I=', i, 'J=', j, A(i,j), num
            END IF
         END DO
      END DO

    END SUBROUTINE TzStruct

!!!!!!!!
    SUBROUTINE Tzlarge( A, LDA, N, K ) 
! ..
! ..Scalar Arguments..
      INTEGER, INTENT(IN) :: LDA, N, K
! ..
! ..Array Arguments..
      DOUBLE PRECISION, INTENT(IN) :: A(LDA,*)
!
! Purpose
! =======
! 
! Test the structure of matrix A, and whether there is a big entry in
! the K-th row of A, larger than TOL in magnitude. 
!
! ==========
! Modified on April 16th, 2013
! ============================
!
      INTEGER j
      DOUBLE PRECISION  TOL

      TOL  = 1.0E-10
      Write(*,*) 'Large in ', K, '-th row'
      DO j = 1, N
         IF( abs(A(K,j)) > TOL ) &
                 write(*,*) 'I=', K, 'J=', j, A(K,j)
      END DO
         
    END SUBROUTINE Tzlarge

!!!!!!!!
  SUBROUTINE Cauchy_svd( A,X,Y,Z,N )
! 
! Construction an Cauchy matrix A_{ij} = z_j / (x_j^2 - y_i^2)
! 
! X column basis
! Y row basis
! Z nominator
!
! This procedure is design for updating SVD problem 
!===============
! Written by Shengguo Li, On Aug. 16th, 2012
!===============

    REAL(8)  :: A(N,*), X(*), Y(*), Z(*), temp
    INTEGER  :: N       ! Leading dim of A
!
    INTEGER  :: I,J
    
    DO J = 1, N
       DO I = 1, N
          temp = ( X(J)-Y(I) ) *( X(J) + Y(I) )
          A(I,J) = Z(J) / temp
       END DO
    END DO
    
  END SUBROUTINE Cauchy_svd 

!!!! NormlzM
  subroutine NormLzM(A, N, RC)
!
! .. Scalar Arguments ..
    INTEGER   :: N
    CHARACTER(LEN = 3) :: RC

! .. Array Arguments ..
    DOUBLE PRECISION  :: A(N,N)
!
! Purpose
! ========
! Normalize the column or row of matrix A. 
!
! .. Parameters .. 
! A  (INPUT/OUTPUT)  DOUBLE PRECISION Array, DIMENSION(N,N) 
! 
! N  (INPUT)  INTEGER
!
! RC (INPUT)  CHARACTER STRING, 'ROW' or 'COL'.
!
! =================
! Written by Shengguo Li, On Aug. 16th, 2012
! =================

    INTEGER :: I
    DOUBLE PRECISION  :: TEMP, dnrm2

!  .. External Functions
      LOGICAL     LSAME, GR, GC
      EXTERNAL    LSAME, dnrm2

    GR = LSAME(RC, 'ROW')
    GC = LSAME(RC, 'COL')

    IF( .not.(GR .or. GC) ) THEN
       print *, 'Normalization must specify row or col. '
    END IF
    
    IF( GC ) THEN    
       DO I =1, N
          temp = dnrm2(N,A(1:N, I),1 )
          temp = 1.0D0 / temp
          A(1:N,I) = temp * A(1:N,I)
       END DO
    ELSE
       DO I =1, N
          temp = dnrm2(N,A(I, 1:N),1 )
          temp = 1.0D0 / temp
          A(I, 1:N) = temp * A(I, 1:N)
       END DO
    END IF
       
  end subroutine NormLzM

!!!!!!
  subroutine searchMax(U,W,D2,D1,LDU,LDW,ii,jj,Nmax)
! search the maximum magnitude entry: first column then row
! The track looks like a snake. It can save flops by half but may not
! get the largest entry. 

! This is designed for Cauchy-HSS method, not used any more. 

    double precision U(*),W(*),D2(*),D1(*)
    integer LDU,LDW,ii,jj

    double precision zero,piv,junkL,junkU,Nmax
    integer pivot,jju,jjL,flg
    parameter(zero = 0.0D0)
    
    pivot = 0
    call CauchyMax(D2,D1(1),u,LDU,junkL,jjL)
    junkL = junkL*abs(w(1))       
    call CauchyMax(D1,D2(1),w,LDW,junkU,jju)
    junkU = junkU*abs(u(1))
    piv = abs(u(1)*w(1)/(D2(1)-D1(1)) )
    if(junkL .le. piv .and. junkU .le. piv ) then
       ii = 1
       jj = 1
       Nmax = piv
!       write(*,*) 'pivot = ', pivot
       return
    end if
    
    if(junkL .gt. junkU) then
       flg = 1
       jju = 1
    else
       flg = 0
       jjL = 1
    end if

    do while( 1 .lt. 2) ! ii: col, jj: row, ii < jj
       pivot = pivot + 1
       if(flg .eq. 1) then
          ii = jju
          call CauchyMax(D1,D2(jjL),w,LDW,junkU,jju)
          if(jju .eq. ii) then 
             jj = jjL
             Nmax = junkU*abs(u(jj))
 !            write(*,*) 'pivot = ', pivot
             return
          end if
          flg = 0
       end if
!
       jj = jjL
       call CauchyMax(D2,D1(jju),u,LDU,junkL,jjL)
       if(jjL .eq. jj) then
          ii = jju
          Nmax = junkL*abs(w(ii))
!          write(*,*) 'pivot = ', pivot
          return
       end if
       flg = 1

    end do ! while

  end subroutine searchMax

!!!!!!
  subroutine CauchyMax(D,F,u,N,junk,jj)

    double precision D(*),u(*),junk,F
    integer N,jj
    integer temp(1)

!  .. Intrinsic Functions ..
    intrinsic    Maxloc, ABS,MAXVAL

    double precision, allocatable :: LK(:) ! change it to LK(N)

    allocate( LK(N) )

    Lk = u(1:N)/ ( D(1:N)-F )
    Lk = abs(Lk)
    junk = maxval(Lk(1:N))
    temp = maxloc(Lk(1:N))
    jj = temp(1)
    
    deallocate(LK)
    
  end subroutine CauchyMax

!!!!!!
      subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************
!
!! R8MAT_PRINT prints a R8MAT.
!
!  Discussion:
!
!    A R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return

end subroutine r8mat_print

!!!!!!!!
   subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************
!
!! R8MAT_PRINT_SOME prints some of a R8MAT.
!
!  Discussion:
!
!    A R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end subroutine r8mat_print_some
   
!!!!!!!!!!!!
   subroutine r8vec_print ( n, a, title )

!*****************************************************************************
!
!! R8VEC_PRINT prints a R8VEC.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end subroutine r8vec_print

!!!!!!!!!!!
    subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to print.
!    The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i_hi
  integer   ( kind = 4 ) i_lo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
  end do

  return
end subroutine r8vec_print_some

!!!!!!!!
   subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end subroutine i4vec_print

!!!!!!!!!
  subroutine i4vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_SOME prints "some" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to print.
!    The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i_hi
  integer   ( kind = 4 ) i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,i8)' ) i, a(i)
  end do

  return
end subroutine i4vec_print_some

!!!!!!
subroutine init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
          
  CALL SYSTEM_CLOCK(COUNT=clock)
          
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
          
  DEALLOCATE(seed)

end subroutine init_random_seed


end module BasicMM
