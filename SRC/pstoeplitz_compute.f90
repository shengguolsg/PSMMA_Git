!
      SUBROUTINE PSTOEPLITZ_COMPUTE( N,K, ALPHA,A,LDA,DESCA,B,LDB,&
                  BETA, C, LDC, IMROW, IMCOL, WORK )
!
      IMPLICIT NONE                        
!
!  -- PSMMA Package routine (version 0.1) --
!     We assume that A is a general matrix and B is a structured matrix with O(N) parameters.
!     Here B is a Toeplitz matrix, which is defined by two vectors, stored in 
!     the matrix B(LDB,2). 
!     We assume that every processes has a local copy of matrix B, the generators. 
!
!     Written by Shengguo Li, 2022-06-26
!     National University of Defense Technology, Changsha 410073, China. 
!
!   =======================================                                                                      
!
!     .. Scalar Arguments ..
      INTEGER            IMROW, IMCOL, LDA, LDB, LDC
      INTEGER            N, K
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA(*)
      DOUBLE PRECISION             :: A( LDA, * ), C( LDC,* )
      DOUBLE PRECISION, INTENT(IN) :: B(*)
      DOUBLE PRECISION   WORK( * )
!     ..
!
!  Parameters
!  ==========
!
!  N      - (input) INTEGER
!           N specifies the (global) number of rows of the matrix A
!           and of the matrix C.  N >= 0.
!
!  K      - (input) INTEGER
!           K specifies the (global) number of columns of the matrix A
!           and the (global) number of rows of the matrix B.  K >= 0.
!
!  MB     - (input) INTEGER
!           MB specifies the row block size of the matrix A and of the
!           matrix C.  MB >= 1.
!
!  NB     - (input) INTEGER
!           NB specifies the column block size of the matrix B and of
!           the matrix C.  NB >= 1.
!
!  KB     - (input) INTEGER
!           KB specifies the column block size of the matrix A and the
!           row block size of the matrix B.  KB >= 1.
!
!  ALPHA  - (input) DOUBLE PRECISION
!           ALPHA specifies the scalar alpha.
!
!  A      - (input) DOUBLE PRECISION array of DIMENSION ( LDA, Kq ).
!           The leading Mp by Kq part of the array A must contain the
!           (local) matrix A.  Mp and Kq are local variables  (see de-
!           scription of local parameters).
!
! ICTXT     (input) INTEGER
!            The Process grid CONTEXT      
!            
!  LDA    - (input) INTEGER
!           The leading dimension of the (local) array A.
!           LDA >= max( 1, Mp ).
!
!  B      - (input) DOUBLE PRECISION array of DIMENSION ( LDB ).
!           The first K entries stores the first column entries of a Toeplitz 
!           matrix from bottom to top and B(K) stores the diagonal, and the 
!           first row entries of Toeplitz is stored in B(K+1:end)
!      
!  LDB    - (input) INTEGER
!           The length of the array , LDB = 2*K-1 when B is unsymmetric and 
!           LDB = K when B is symmetric. 
!      
!  BETA   - (input) DOUBLE PRECISION
!           BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!
!  C      - (input/output) DOUBLE PRECISION array of DIMENSION (LDC,Nq)
!           On entry, the leading Mp by Nq part of the array C must
!           contain the (local) matrix C except when  beta  is zero,
!           in which case C need not be set on entry.  Mp and Nq are
!           local variables (see description of local parameters).
!           On exit, the array C is overwritten by the Mp by Nq matrix
!           ( alpha*A*B + beta*C ).
!
!  LDC    - (input) INTEGER
!           The leading dimension of the (local) array C.
!           LDC >= max( 1, Mp ).
!
!  IMROW  - (input) INTEGER
!           IMROW specifies a row of the process template, which holds
!           the first block of the matrices.  0 <= IMROW < NPROW.
!
!  IMCOL  - (input) INTEGER
!           IMCOL specifies a column of the process template, which
!           holds the first block of the matrices.  0 <= IMCOL < NPCOL.
!
!  WORK   - (workspace) DOUBLE PRECISION array
!           See requirements.
!
!  IWORK  - (workspace) INTEGER array
!           See requirements.
!
! ===============
!  June 22, 2022. 
!
!  1) Currently we assume that VECTOR B is of length 2K-1. i.e., both the 
!  lower and upper triangular generators are stored. 
!
!  2) Low-rank approximation is implemented in this routine. 
!
! July 3th, 2022
! 
! 1) A : N-by-K, B: K-by-K, C: N-by-K. Essentially we assume N == K.
!       
! 2) When NPROW .NE. NPCOL, the MB, NB and KB for BDD may not be the same.        
!  
!  =================
!
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_, &
                         MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                         CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                         RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE,           ZERO
      PARAMETER          ( ONE  = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER              INFO, IPA, IPB, ICTXT
      INTEGER              K1, KAQ, KBP, Rk, KP, NB, KB, MB
      INTEGER              MP, MRCOL, MRROW, MYCOL, MYROW, LengthC
      INTEGER              NPCOL, NPROW, NROLL, NQ, ierr, LengthR, &
                           IPDW, IPB2
      DOUBLE PRECISION     TBETA
      LOGICAL              islowrank      
!     ..
!     .. Local Arrays ..      
      INTEGER, ALLOCATABLE :: RIndex(:), CIndex(:) 
!     Record the Row and Column indexes of generators
!     ..
!     .. External Subroutines ..
      EXTERNAL             PXERBLA, DGEMM, is_intersect
      EXTERNAL             BLACS_GRIDINFO, DGESD2D, DGERV2D
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC            MAX, MOD
!     ..
!     .. External Functions ..
      INTEGER, EXTERNAL    :: NUMROC 
!     ..
!     ..
!     .. Executable Statements ..
!
!     Get grid parameters
!
!     MB : row block of A and C
!     KB : col block of A, row block of B
!     NB : col block of B and C

      ICTXT = DESCA( CTXT_ )
      MB = DESCA( MB_ )  
      KB = DESCA( NB_ ) ! We assume MB == NB
      NB = KB

      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
!
!     Test for the input parameters.
!
      INFO = 0
      IF( N  .LT. 0                      ) THEN
         INFO = 2
      ELSE IF( K  .LT. 0                      ) THEN
         INFO = 3
      ELSE IF( NB .LT. 1                      ) THEN
         INFO = 5
      ELSE IF( KB .LT. 1                      ) THEN
         INFO = 6
      ELSE IF( LDA.LT. 1                      ) THEN
         INFO = 9
      ELSE IF( LDC.LT. 1                      ) THEN
         INFO = 14
      ELSE IF( IMROW.LT.0 .OR. IMROW.GE.NPROW ) THEN
         INFO = 15
      ELSE IF( IMCOL.LT.0 .OR. IMCOL.GE.NPCOL ) THEN
         INFO = 16
      END IF
!
   10 CONTINUE
      IF( INFO .NE. 0 ) THEN
          CALL PXERBLA( ICTXT, 'PSTOEPLITZ_COMPUTE ', INFO )
          RETURN
      END IF
!
      MRROW = MOD( NPROW+MYROW-IMROW, NPROW )
      MRCOL = MOD( NPCOL+MYCOL-IMCOL, NPCOL )
!
      ! ********* For some process, MP and NQ can be ZERO. *******
      MP = NUMROC( N, NB, MRROW, 0, NPROW )
      NQ = NUMROC( K, NB, MRCOL, 0, NPCOL )
      KP = NUMROC( K, KB, MRROW, 0, NPROW )

      IF( MP.EQ.0 ) GOTO 900
!
!     Test for the input parameters again with local parameters
!
      IF(      LDA .LT. MP ) THEN
         INFO = 9
      ELSE IF( LDC .LT. MP ) THEN
         INFO = 14
      END IF
      IF( INFO .NE. 0 ) GO TO 10
!
!     Quick return if possible.
!
      IF( ( N.EQ.0 ).OR. &
        ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
         RETURN
!
      TBETA  = BETA
      KAQ    = NUMROC( K, KB, MRCOL, 0, NPCOL )
!     Note that KAQ == NQ.       
!
!     Compute pointer addresses for buffer
!
      IPA  = 1
      IPB  = IPA +  MP*NQ
      IPDW = IPB +  2*MP*NQ  ! only for without low-rank
!
!     Copy all blocks of A to WORK(IPA) with column block presorting
!
!      ===================================       
      CALL DLACPY( 'A', MP, NQ, A, LDA, WORK(IPA), MP  )

!     Allocate the Row and Column index of generators
      ALLOCATE( RIndex(NQ+NB), CIndex(NQ),stat=ierr )     
      IF( ierr .ne. 0 ) THEN
            WRITE(*,*) 'Allocation failed in pdmnnbeig, and return'
            RETURN
      END IF       
!
      CALL ConstIndex( K,NB,MRCOL,NPCOL,CIndex,LengthC )
!
!     Main Multiplication Routine
!
      DO 40 K1 = 0, NPCOL-1
!
         NROLL = MOD( MRCOL+K1,NPCOL )
         CALL ConstIndex( K,KB,NROLL,NPCOL,RIndex,LengthR )
         KBP = LengthR
!
         IF( LengthC.EQ.0 .OR. LengthR.EQ.0 ) GOTO 80
!         
!        Check whether the intersections of RIndex and CIndex is empty.         
!        If empty, use low-rank approximation; Different from Cauchy-like matrix,
!        it need to construct the full matrix first and then use RRQR to construct
!        low-rank approximation. Then, perform low-rank matrix-matix multiplication. 
!         
!        Otherwise, treat it as a dense matrix.
         call is_intersect( LengthC,CIndex,LengthR,RIndex,islowrank )
!         islowrank = .false.
!
         IF ( islowrank ) THEN
            ! Form low-rank approximation
            CALL ConstToeplitzlowrank( KBP,NQ,B,LDB,WORK(IPB),RIndex,&
                  CIndex,Rk )
            !IF (MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            !     WRITE(*,*) "MYROW,MYCOL,Rk=", MYROW,MYCOL,KBP,NQ,Rk
            !END IF
!
             IF( Rk == MIN(KBP,NQ) ) THEN
                  CALL DGEMM( 'N', 'N', MP, NQ, KBP,ALPHA,WORK(IPA), &
                    MAX(1,MP),WORK(IPB),MAX(1,KBP),TBETA, C,LDC )
             ELSE                  
!                 Perform two low-rank matrix multiplications
                  CALL DGEMM('N','N',MP,Rk,KBP,ONE,WORK(IPA),MAX(1,MP),&
                        WORK(IPB),KBP,ZERO,WORK(IPDW),MP  )             
!                 MP == KBP
                  IPB2 = IPB+KBP*Rk 
                   CALL DGEMM( 'N','N',MP,NQ,Rk,ALPHA,WORK(IPDW),MP,&
                         WORK(IPB2), Rk, TBETA, C, LDC )
             END IF
             TBETA = ONE
!
         ELSE
            ! Construct full local submatrix B
            CALL ConstToeplitz( KBP,NQ,B,LDB,K,WORK(IPB),RIndex,CIndex )

            CALL DGEMM( 'N', 'N', MP, NQ, KBP,ALPHA,WORK(IPA), &
                    MAX(1,MP),WORK(IPB),MAX(1,KBP),TBETA, C,LDC )
            TBETA = ONE
!
         END IF

!        Shift A (WORK(IPA)) to the left
!
   80    IF( K1.LT.NPCOL-1 ) THEN
!
!            Rotate A (WORK(IPA))
!
             CALL DGESD2D( ICTXT, MP, KAQ, WORK(IPA), MP, MYROW, &
                           MOD(NPCOL+MYCOL-1, NPCOL) )
             KAQ = NUMROC( K, KB, MRCOL+K1+1, 0, NPCOL )
             CALL DGERV2D( ICTXT, MP, KAQ, WORK(IPA), MP, MYROW, &
                           MOD(MYCOL+1, NPCOL) )
         END IF
   40 CONTINUE

       DEALLOCATE(RIndex, CIndex)
!
  900  CONTINUE

      RETURN
!
!     End of PSTOEPLITZ_COMPUTE
!
    END SUBROUTINE PSTOEPLITZ_COMPUTE
