!
      SUBROUTINE PSCAUCHYEIG_COMPUTE(N,K,NB,KB,ALPHA,A,LDA,DESCA,F,D,U,&
                   V,DIFL,DIFR,BETA,C,LDC,IMROW,IMCOL,WORK )
!
      IMPLICIT NONE                        
!
!  -- PSMMA Package routine (version 0.1) --
!     This routine is designed for structured DC algorithm. We assume that
!     A is a general matrix and B is a structured matrix with O(N) parameters.
!     Here B is a Cauchy-like matrix, which is defined by SIX vectors, stored in 
!     F,D,U,V,DIFL,DIFR, and B(i,j) = U(i)*V(j) / (F(i)-D(j)). 
!     We assume that every processes has a local copy of matrix B, the generators. 
!
!     Written by Shengguo Li, 2020-06-26
!     National University of Defense Technology, Changsha 410073, China. 
!
!   =======================================                                                                      
!
!     .. Scalar Arguments ..
      INTEGER            IMROW, IMCOL, LDA, LDC
      INTEGER            N, K, NB,KB
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA(*)
      DOUBLE PRECISION   A( LDA, * ), C( LDC,* ), DIFR(*), D(*)
      DOUBLE PRECISION   WORK( * ), F(*), U(*), V(*), DIFL(*)
!     ..
!
!  Purpose
!  =======
!
!  PSCAUCHYEIG_COMPUTE routine is one of the PSMMA package based on block cyclic
!  data distribution (BCDD) on 2-D process configuration. This routine does not redistribute
!  matrix A from BCDD to BDD (block data distribution). 
!  This routine is designed for the tridiagonal DC algorithm, and assume that 
!  DESCC == DESCA and A is N-by-K, B is K-by-K, and C is N-by-K. 
!
!  PSCAUCHYEIG_COMPUTE performs the matrix-matrix operations
!
!     C := alpha*A*B + beta*C,
!
!  where alpha and beta are scalars, and A, B and C are matrices,
!  with A an M by K matrix (globally), B a K by N matrix (globally)
!  and C an M by N matrix (globally).
!
!  The basic principles are to simulate "Columnwise Broadcast B"
!                       and perform "Rowwise Shift A"
!
!  Since matrix B is stored by its generators and every process has a copy of B,
!  only matrix A is needed to be communicated. Therefore, we can simulate the matrix
!  multiplication easily even in block cyclic data distribution form. We do not 
!  need to perform PUMMA algorithm, which makes the algorithm too complicated. 
!  Comparing with PUMMA, we don't need to permute matrix A, but we need to shift
!  A to left too. We also need to workspace to construct the local submatrix of B. 
! 
!  Question: Can we use less workspace than PUMMA ? 
!
!  The main points are 1) the column indexes of B are fixed for each process;
!  2) the row indexes of B are updated from the column indexes of A after shift. 
!   
!  The whole algorithm executes NPCOL times and only matrix A shifts leftward. 
!  This routine can also deal with the cases that IMROW and IMCOL is not ZERO. 
!      
!  This is a function, and the input parameters are
!   1) K, N  : The row and column dimension of B, which is not assumed to be square; 
!   2) MB:  The row block size of matrix A and C;
!   3) KB:  The column block size of matrix A and row block size of B;
!   4) NB:  The column block size of matrix B and C;
!   3) NPROW, NPCOL: The shape of the process grid. 
!                     
!            
!  Parameters
!  ==========
!
!  M      - (input) INTEGER
!           M specifies the (global) number of rows of the matrix A and
!           of the matrix C.  M >= 0.
!
!  N      - (input) INTEGER
!           N specifies the (global) number of columns of the matrix B
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
!  B      - (input) DOUBLE PRECISION array of DIMENSION ( LDB, NCB ).
!           The leading N and NCB rows and columns store the generators
!           of the Cauchy-like matrix.
!
!  LDB    - (input) INTEGER
!           The leading dimension of the array B.
!
!  NCB    - (input) INTEGER
!           The number of columns of matrix B. It is Four.                   
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
!  Local  Parameters
!  =================
!
!  MpxKq = size of (local) matrix A in the process, iam
!  KpxNq = size of (local) matrix B in the process, iam
!  MpxNq = size of (local) matrix C in the process, iam
!  ITEM  =  temporal integer parameter
!
!  Three buffers for A and B
!     WORK(IPA) <== A
!     WORK(IPB) <== B
!
!  One interger buffer
!     IWORK(k) <== its starting point of column block of TA
!
!  Requirements (approximate)
!  ==========================
!
!   Size(IWORK) = LCM(P,Q) / Q
!               or P for the worst case
!   Size(WORK)  = Size(A) + Size(B)
!   where
!     MG = Ceil( M, MB )
!     NG = Ceil( N, NB )
!     KG = Ceil( K, KB )
!     Mp <= Mp0 = Ceil(MG,P)xMB
!     Nq <= Nq0 = Ceil(NG,Q)xNB
!     Size(A) = Mp x Ceil(KG,Q)xKB
!     Size(B) = Ceil(KG,LCM)xKB x Nq
!
! ======================================================================
!
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,&
                         MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                         CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,  &
                         RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE,           ZERO
      PARAMETER          ( ONE  = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER              INFO, IPA, IPB, ICTXT
      INTEGER              K1, KAQ, KBP, Rk, KP
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
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

!      IF(MYROW.EQ.0 .AND. MYCOL.EQ.0 ) & 
!         WRITE(*,*) "pscauchyeig_compute, K,KB,NB=", K,KB,NB,ICTXT
!
!     Test for the input parameters.
!
      INFO = 0
      IF( N  .LT. 0                      ) THEN
         INFO = 1
      ELSE IF( K  .LT. 0                      ) THEN
         INFO = 2
      ELSE IF( NB .LT. 1                      ) THEN
         INFO = 3
      ELSE IF( KB .LT. 1                      ) THEN
         INFO = 4
      ELSE IF( LDA.LT. 1                      ) THEN
         INFO = 7
      ELSE IF( LDC.LT. 1                      ) THEN
         INFO = 17
      ELSE IF( IMROW.LT.0 .OR. IMROW.GE.NPROW ) THEN
         INFO = 18
      ELSE IF( IMCOL.LT.0 .OR. IMCOL.GE.NPCOL ) THEN
         INFO = 19
      END IF
!
   10 CONTINUE
      IF( INFO .NE. 0 ) THEN
          CALL PXERBLA( ICTXT, 'PSCAUCHYEIG_COMPUTE ', INFO )
          RETURN
      END IF
!
!     Initialize parameters
!
!     MRROW is the logical process coordinate corresponding to the physical
!     process MYROW, in other words, if the matrices started in the physical
!     process (3,1), then for this process MRROW = 0, for the physical
!     processes (2,*), MRROW = NPROW-1,
!     MRCOL is the equivalent of MRROW for the process columns. 
!
!     For a PxQ process template, the logical process (i,j) contains the
!     submatrices (global indices)
!
!           A(    i,j) A(    i,j+Q) ... A(    i,j+k*Q) ...
!           A(  i+P,j) A(  i+P,j+Q) ... A(  i+P,j+k*Q) ...
!               .            .      ...        .       ...
!           A(i+l*P,j) A(i+l*P,j+Q) ... A(i+l*P,j+k*Q) ...
!               .            .      ...        .       ...
!
!     IGE is the smallest power of 2 greater or equal than IGD.
!
      MRROW = MOD( NPROW+MYROW-IMROW, NPROW )
      MRCOL = MOD( NPCOL+MYCOL-IMCOL, NPCOL )
!
!     MP is the number of local rows of A (and C) I own,
!     NQ is the number of local columns of B (and C) I own,
!     KP is the number of local rows of B I own,
!     KQ is the number of local columns of B I own.
!     These statements are equivalent to
!           MP = NUMROC( M, MB, MYROW, IMROW, NPROW )
!           NQ = NUMROC( N, NB, MYCOL, IMCOL, NPCOL )
!           KP = NUMROC( K, KB, MYROW, IMROW, NPROW )
!           KQ = NUMROC( K, KB, MYCOL, IMCOL, NPCOL )
!
      ! The local A submatrix is initially MP-by-NQ. Its row dimension is fixed, but
      ! its column dimension is varied since it is shift leftwards. 
      ! The local B submatrix is X-by-NQ. Its columns dimension is fixed, but its row
      ! dimension is determined by the corresponding submatrix A. 
      ! Only MP and NQ matters. KP is useless. 
      ! 
      ! ********* For some process, MP and NQ can be ZERO. *******
      MP = NUMROC( N, NB, MRROW, 0, NPROW )
      NQ = NUMROC( K, KB, MRCOL, 0, NPCOL )
      KP = NUMROC( K, KB, MRROW, 0, NPROW )

!      WRITE(*,*) "MP,NQ,KP=", MP,NQ,KP,MYROW,MYCOL,MRROW,MRCOL

      IF( MP.EQ.0 ) GOTO 900
!
!     Test for the input parameters again with local parameters

      IF(      LDA .LT. MP ) THEN
         INFO = 7
      ELSE IF( LDC .LT. MP ) THEN
         INFO = 17
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
      IPA = 1
      IPB = MP * (KAQ+KB) + 1
!      IPDW = IPB + 2*KP*NQ
      IPDW = IPB + (NQ+NB)*NQ  ! only for without low-rank
!
!     Copy all blocks of A to WORK(IPA) with column block presorting
!
!      ****** Beware NQ can be ZERO ******
!      ===================================       
      CALL DLACPY( 'A', MP, NQ, A, LDA, WORK(IPA), MP  )
!      CALL PDLACPY( 'A',N,K,A,DESCA,WORK(IPA),DESCC )
!      CALL PDtestOrth1( N,K,WORK(IPA),1,1,DESCC )

!     Allocate the Row and Column index of generators
      ALLOCATE( RIndex(NQ+NB), CIndex(NQ),stat=ierr )     
      IF( ierr .ne. 0 ) THEN
            WRITE(*,*) 'Allocation failed in pdmnnbeig, and return'
            RETURN
      END IF       
!
!     The column indexes of generators of B are static. B is K-by-N, and the
!     the local B is KBP-by-NQ. 
!     ============================================
!     ****** Beware LengthC==NQ can be ZERO ******
!     ============================================
      CALL ConstIndex( K,KB, MRCOL,NPCOL,CIndex,LengthC )
!      WRITE(*,*) "After ConstIndex, CIndex ", LengthC, MYROW, MYCOL
!      WRITE(*,*) "After ConstIndex, CIndex()=", MYROW,MYCOL,CIndex(1:LengthC)
!
!     Main Multiplication Routine
!
      DO 40 K1 = 0, NPCOL-1
!
!        Construct the local row indexes of B based on the column
!        indexes of current matrix A
         NROLL = MOD( MRCOL+K1,NPCOL )
         CALL ConstIndex( K,KB,NROLL,NPCOL,RIndex,LengthR )
         KBP = LengthR
!         WRITE(*,*) "After RIndex(1:LengthR)", LengthR, MYROW, MYCOL, NROLL
         IF( LengthC.EQ.0 .OR. LengthR.EQ. 0) GOTO 80
!         
!        Check whether the intersections of RIndex and CIndex is empty.         
!        If empty, use low-rank approximation; Otherwise, treat it as a dense matrix.
         call is_intersect( LengthC,CIndex,LengthR,RIndex,islowrank )
!         WRITE(*,*) "Is lowrank:", MYROW,MYCOL, islowrank
!         islowrank = .false.

         IF ( islowrank ) THEN
            ! Form low-rank approximation
            CALL ConstCauchyEiglowrank( KBP,NQ,WORK(IPB),F,D,U,V,DIFL,&
                   DIFR, RIndex,CIndex,Rk )
!             WRITE(*,*) "MYROW,MYCOL,Rk=", MYROW,MYCOL,KBP,NQ,LengthC,Rk
!            IF( MYROW.EQ.0 .AND. MYCOL.EQ.1 .AND. K1.EQ.1 ) THEN
!                  WRITE(*,*) "RIndex=", RIndex(1:LengthR)
!                  WRITE(*,*) "CIndex=", CIndex(1:LengthC)
!            END IF
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
            CALL ConstCauchyEig( KBP,NQ,WORK(IPB),F,U,V,DIFL,DIFR,&
                  RIndex,CIndex )
!
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
!     End of PSCAUCHYEIG_COMPUTE
!
      END SUBROUTINE PSCAUCHYEIG_COMPUTE
