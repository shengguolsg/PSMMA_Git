      SUBROUTINE PZMDFTL1( N, K, NB, KB, ALPHA, B, LDB, BETA, C, LDC, &
                           IMROW, IMCOL, WORK, RWORK )
!
!        use MPI
        IMPLICIT NONE
        
        include  'mpif.h'
!
!  -- PDMDFTM 
!
!     .. Scalar Arguments ..
      INTEGER            IMROW, IMCOL, K, KB, LDB, LDC
      INTEGER            N, NB
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         B(LDB, *), WORK( * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  PDMDFTL routine is based on block data distribution on 2-D process 
!  configuration. It is designed for computing Fourier Transform, C=F*B,
!  where F is DFT matrix with dimension N, square matrix.  
!
!  PDMNNA performs the matrix-matrix operations
!
!     C := alpha*F*B + beta*C,
!
!  where alpha and beta are scalars, and B and C are matrices,
!  with F an N by N matrix (globally), B a N by K matrix (globally)
!  and C an N by K matrix (globally).
!
!  The basic principles are "Rowwise Broadcast A"
!                       and "Columnwise Shift  B".
!
!  Compare with PZMDFTL, this routine uses RRQR to compute the low-rank
!  approximation, and vectorized initialization. Currently, we use the
!  LAPACK routine ZGEQP3 and it can be optimized later. 
!
!  Parameters
!  ==========
!
!  N      - (input) INTEGER
!           N specifies the (global) number of rows of the matrix A and
!           of the matrix C.  N >= 0.
!
!  K      - (input) INTEGER
!           N specifies the (global) number of columns of the matrix B
!           and of the matrix C.  N >= 0.
!
!  NB     - (input) INTEGER
!           NB specifies the row block size of the matrix A and of the
!           matrix C.  NB >= 1. It also specifies the column block size
!           of the matrix A and the row block size of the matrix B. 
!
!  KB     - (input) INTEGER
!           NB specifies the column block size of the matrix B and of
!           the matrix C.  KB >= 1.
!
!  ALPHA  - (input) DOUBLE PRECISION
!           ALPHA specifies the scalar alpha.
!
!  B      - (input) DOUBLE PRECISION array of DIMENSION ( LDB, Nq ).
!           The leading Kp by Nq part of the array B must contain the
!           (local) matrix B.  Kp and Nq are local variables  (see de-
!           scription of local parameters).
!
!  LDB    - (input) INTEGER
!           The leading dimension of the (local) array B.
!           LDB >= max( 1, Kp ).
!
!  BETA   - (input) DOUBLE PRECISION
!           BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!
!  C      - (input/output) DOUBLE PRECISION array of DIMENSION (LDC,Nq)
!           On entry, the leading Mp by Nq part of the array C must
!           contain the (local) matrix C except when  beta  is zero,
!           in which case C need not be set on entry.  Mp and Nq are
!           local variables  (see description of  local parameters).
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
!  Local  Parameters
!  =================
!
!  LCM   =  the lowest common multiple of P and Q
!  LCMP  =  LCM/P = number of template rows in LCM block
!  LCMQ  =  LCM/Q = number of template columns in LCM block
!  IGCD  =  the Greatest Common Divisor of P and Q
!  NROLL =  it specifies rolling index of B
!  MpxKq = size of (local) matrix A in the process, iam
!  KpxNq = size of (local) matrix B in the process, iam
!  MpxNq = size of (local) matrix C in the process, iam
!  ITEM  =  temporal integer parameter
!
!  Three buffers for A and B
!     WORK(IPA) <== A
!     WORK(IPB) <== B
!
!  One integer buffer
!     IWORK(k) <== its starting point of row block of TB
!
!  Requirements (approximate)
!  ==========================
!
!   Size(IWORK) = LCM(P,Q) / P,
!                 or Q for the worst case
!   Size(WORK)  = Size(A) + Size(B)
!   where
!     MG = Ceil( M, MB )
!     NG = Ceil( N, NB )
!     KG = Ceil( K, KB )
!     Mp <= Mp0 = Ceil(MG,P)xMB
!     Nq <= Nq0 = Ceil(NG,Q)xNB
!     Size(A) = Mp x Ceil(KG,LCM)xKB
!     Size(B) = Ceil(KG,P)xKB x Nq
!
! ======================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION    ONE, NEONE, ZERO, TWO 
      PARAMETER          ( ONE=1.0E+0, ZERO=0.0E+0, NEONE=-1.0E+0, TWO=2.0E+0 )
      COMPLEX*16, PARAMETER :: CI=(ZERO,ONE), CONE=(ONE,ZERO), CZERO=(ZERO,ZERO)
      DOUBLE PRECISION      :: PI 
!     ..
!     .. Local Scalars ..
      INTEGER              IBCOL, ICTXT, IGCD, INFO, IPB
      INTEGER              K1, KBP, KG, KITER, II, I2
      INTEGER              KP, KT, LCM, LCMP, LCMQ, LMPB, IINFO
      INTEGER              LMQB, MP, MRCOL, MRROW, MYCOL, MYROW, IFTOP
      INTEGER              NPCOL, NPROW, NROLL, NQ, MB, M, PowK, &
                           GAPR, GAPC, ICOL, IDL, IDR, IPR, IUC, IWAVT,&
                           IWA,IWTEMP,IWRT, IWRT1, J2, RK, &
                           LWORK, NB1, ITAU, IWORK
      DOUBLE PRECISION     TBETA, TOL, time, time1, time2
      LOGICAL              FORWRD
      COMPLEX*16           Wn
!     ..
!     .. Local Arrays ..
      INTEGER, ALLOCATABLE    :: JPVT(:)
      COMPLEX*16, ALLOCATABLE :: ZDWORK(:)
!     ..
!     .. External Subroutines ..
      EXTERNAL             BLACS_GRIDINFO, DGEMM, DPRESRT, DPRESS
      EXTERNAL             DBCSTGCD, DGESD2D, DGERV2D, &
                           PXERBLA, ZLASCL2, ZLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC            MAX, MOD
!     ..
!     .. External Functions ..
      INTEGER              ICEIL, ILCM, IEXPV, NUMROC, NPREROC
      EXTERNAL             ICEIL, ILCM, IEXPV, NUMROC, NPREROC
      DOUBLE PRECISION     DZNRM2
      EXTERNAL             DZNRM2
      
!     ..
!     .. Common blocks ..
      COMMON             / CONTEXT / ICTXT
!
!     .. Executable Statements ..
!
!     Get grid parameters
!
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
!
!     Test the input arguments
!
      M = N
      MB = NB
      INFO = 0
      IF( N  .LT. 0                      ) THEN
         INFO = 1
      ELSE IF( K  .LT. 0                      ) THEN
         INFO = 2
      ELSE IF( NB .LT. 1                      ) THEN
         INFO = 3
      ELSE IF( KB .LT. 1                      ) THEN
         INFO = 4
      ELSE IF( LDB.LT. 1                      ) THEN
         INFO = 7 
      ELSE IF( LDC.LT. 1                      ) THEN
         INFO = 10
      ELSE IF( IMROW.LT.0 .OR. IMROW.GE.NPROW ) THEN
         INFO = 11
      ELSE IF( IMCOL.LT.0 .OR. IMCOL.GE.NPCOL ) THEN
         INFO = 12
      END IF
!
   10 CONTINUE
      IF( INFO .NE. 0 ) THEN
          CALL PXERBLA( ICTXT, 'PDMDFTL ', INFO )
          RETURN
      END IF
!
!     Initialize parameters
!
!     MRROW is the logical process coordinate corresponding to the physical
!     process MYROW, in other words, if the matrices started in the physical
!     process (3,1), then for this process MRROW = 0, for the physical
!     processes (2,*), MRROW = NPROW-1,
!     MRCOL is the equivalent of MRROW for the process columns,
!     multiplication can be combined,
!     LCM is the minimal global block stride separating submatrices, which
!     LCMP is the corresponding minimal local block stride along a matrix row,
!     LCMP is also the local number of presorted subblocks to be arranged
!     for rowwise presorting,
!     LCMQ is the corresponding minimal local block stride along a matrix
!     column, is also the local number of presorted subblocks to be arranged
!     for columnwise presorting.
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
!     The blocks which can be combined in a row are A(*,j+k*LCM) and
!     A(i+l*LCM,*) in a column. Locally, A(i,j) = A(0,0),  and
!     A(i+l*P,j+k*Q) = A(l,k). Let A(m,n) be a block in logical process
!     (i,j), where m, n are local indexes, then the blocks, that can be
!     combined with it in a row are given locally by A(m,n+k*LCMQ)
!     (resp. A(m+l*LCMP,n) in a column).
!
!     IGCD is the Greatest Common Diviser (GCD) of (P,Q) LCM*IGCD = NPROW*NPCOL
!     IGE is the smallest power of 2 greater or equal than IGCD.
!
      MRROW = MOD( NPROW+MYROW-IMROW, NPROW )
      MRCOL = MOD( NPCOL+MYCOL-IMCOL, NPCOL )
      LCM   = ILCM( NPROW, NPCOL )
      LCMP  = LCM / NPROW
      LCMQ  = LCM / NPCOL
      IGCD  = NPROW  / LCMQ
!
!     MP is the number of local rows of A (and C) I own,
!     NQ is the number of local columns of B (and C) I own,
!     KP is the number of local rows of B.
!     These statements are equivalent to
!           MP = NUMROC( M, MB, MYROW, IMROW, NPROW )
!           NQ = NUMROC( N, NB, MYCOL, IMCOL, NPCOL )
!           KP = NUMROC( K, KB, MYROW, IMROW, NPROW )
!
      MP = NUMROC( M, MB, MRROW, 0, NPROW )
      NQ = NUMROC( N, NB, MRCOL, 0, NPCOL )
      KP = NUMROC( K, KB, MRROW, 0, NPROW )
!
!     Test the input arguments
!
      IF( LDB .LT. KP ) THEN
         INFO = 7 
      ELSE IF( LDC .LT. MP ) THEN
         INFO = 10
      END IF
      IF( INFO .NE. 0 ) GO TO 10
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
       ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
        RETURN
!
!     LMPB is the number of local rows between two blocks own by one process,
!     which multiplication can be combined,
!     LMQB is the number of local columns between two blocks own by one
!     process, which multiplication can be combined,
!     KG is the total number of blocks used to store K rows (or columns)
!     in KB-blocks,
!     KITER is the maximum number of submatrices, which can be combined at
!     each step of the algorithm,
!     KT is the number rows or columns per combined subblock.
!
!     KBP is the number of local rows of B, KBP = KP.
!
      LMPB   = LCMP * KB
      LMQB   = LCMQ * KB
      KG     = ICEIL( K, KB )
      KITER  = ICEIL( KG, LCM )
      KT     = KITER * KB
      IBCOL  = IGCD * ( MRCOL / IGCD )
      TBETA  = BETA
      KBP    = NUMROC( K, KB, MRROW, 0, NPROW )
!
!     Compute pointer addresses for buffer
!
      IPB  = 1
!
!     Copy all blocks of B to WORK(IPB) with row block presorting
!
!      CALL DPRESRT( 'Row', KP, NQ, KB, B, LDB, RWORK(IPB), KBP, &
!                    KITER, LMPB, LCMP )
!
!    *********************************************************************
!     Compute the low-rank approximation for F_{11}. First, construct 
!     the entries of F_{11}; Then, compute its truncated SVD, based on
!     the accuracy parameter TOL. 
!    *********************************************************************
!
!     The size of WORK is required to be larger than 3*MP*MP+MP+3*MP
!     The size of RWORK requires MP*MP+5*MP+MP. 
!
      NB1   = 16
      IFTOP = 1
      ITAU  = IFTOP+MP*MP
      IWORK = ITAU+MP
  
      PI  = 4.0E+0 * ATAN( ONE )     
      Wn  = EXP( NEONE*TWO*PI*CI/N )
!      
      time = MPI_Wtime() 
      DO I2 = 1, MP
         DO J2 = 1, MP   
           PowK = (I2-1)*(J2-1)
           WORK( MP*(I2-1)+J2 ) = Wn**PowK
         END DO
      END DO
      time2 = MPI_Wtime() 
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
          WRITE(*,*) 'Matrix initialization costs ', time2 - time 
      ENDIF
!
!     Compute RRQR of ZWORK and construct its low-rank approximation
!     For simplicity, we assume that N is dividable by NB, and N= NB*NPROW=NB*NPCOL, 
!     in Version-V0.  
! 	  
      ALLOCATE( JPVT(MP) )
      LWORK = (MP+1)*NB1

      JPVT( 1:MP ) = 0

      time1 = MPI_Wtime()
      CALL ZGEQP3( MP,MP,WORK(IFTOP),MP,JPVT,WORK(ITAU),WORK(IWORK),&
                   LWORK,RWORK,INFO )
      time2 = MPI_Wtime()
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
          WRITE(*,*) 'ZGEQP3 costs ', time2 - time1
      ENDIF

      ! Determine the rank of A
      TOL = 1E-13
      DO II = 1, MP
         IUC = (II-1)*MP+II 
         IF ( ABS( WORK(IUC) ) <= TOL*ABS( WORK(1) ) ) THEN
              EXIT
         ENDIF
      END DO 
      RK = II- 1
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) WRITE(*,*) 'GEQP3 rank is', RK

      ! Construct the orthogonal matrix Q
      IWRT  = IWORK+MP*RK
      LWORK = MP*RK       
      CALL ZLASET( 'A',MP,RK,CZERO,CONE,WORK(IWORK),MP ) 

      CALL ZUNMQR( 'Left','N',MP,RK,RK,WORK(IFTOP),MP,WORK(ITAU),WORK(IWORK),MP,&
                    WORK(IWRT),LWORK,INFO ) 

      ! Construct the R part
      CALL ZLASET( 'A',RK,MP,CZERO,CZERO,WORK(IWRT),RK ) 
      CALL ZLACPY( 'U',RK,MP,WORK(IFTOP),MP,WORK(IWRT),RK )
      ! Update T:=T*P
      FORWRD = .FALSE.
      CALL ZLAPMT( FORWRD,RK,MP,WORK(IWRT),RK,JPVT )
! 
!     Check the correctness of the low-rank approximation computed by
!     RRQR. The following is wrong. WORK(1:MP*MP) has been destroyed. 
!      IWEND = IWRT+MP*RK
!      CALL ZGEMM( 'N','N',MP,MP,RK,ONE,WORK(IWORK),MP,WORK(IWRT),RK,&
!                  ZERO,WORK(IWEND),MP )
!      IF(MYROW.EQ.0 .AND. MYCOL.EQ.0) THEN
!         IWERR = IWEND+MP*MP
!         WORK( IWERR:IWERR+MP*MP-1 ) = WORK(1:MP*MP) - &
!                   WORK( IWEND:IWEND+MP*MP-1 )  
!         Err = dznrm2( MP*MP,WORK(IWERR),1 )
!         WRITE(*,*) 'Low-rank error by RRQR is ', Err 
!      ENDIF 
!	
      ! Copy the low-rank approximation to the front of WORK
      ! Copy Q
      CALL ZLACPY('A',MP,RK,WORK(IWORK),MP,WORK(IFTOP),MP)
  
      ! Copy RT
      IWRT1 = IFTOP+MP*RK
      CALL ZLACPY('A',RK,MP,WORK(IWRT),RK,WORK(IWRT1),RK )
      time2 = MPI_Wtime() 
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
          WRITE(*,*) 'Low-rank totally costs ', time2 - time
      ENDIF

      DEALLOCATE( JPVT )
!
!   *********************************************************************
!	  Compute diagonal scaling matrix DL and DR for each processor 
!   *********************************************************************
      ALLOCATE( ZDWORK(NPCOL*MP+MP) )
  
      ! The DR for all processes in the row are same
          IDR  = 1
          GAPC = MP*MYROW
          ! Store the DR terms
          DO I2 =1, MP
             ZDWORK( IDR +I2 -1) = Wn**( GAPC*(I2-1) )
          END DO
  
          IDL = IDR+MP
          DO II = 0, NPCOL-1
             ICOL = MOD( MYROW+II, NPCOL)
             IPR  = (MYROW*MP)*(ICOL*KP)
             GAPR = MP*ICOL
 
             ! Store the DL terms
             DO I2 = 1, MP
                 ZDWORK( IDL +I2 -1) = Wn**( IPR+ GAPR*(I2-1) )
             END DO
             IDL = IDL +MP
          END DO
          IDL = IDR + MP

!          IF(MYROW.EQ.0 .AND. MYCOL.EQ.0 ) write(*,*) 'DR0 and DL0=', ZDWORK(1:2*MP) 
!
!     Main Multiplication Routine
!
      IWA = IFTOP+MP*RK*2
      IWAVT = IWA+MP*RK
      IWTEMP = IWAVT+MP*RK
      DO K1 = 0, NPROW-1
         NROLL = (K1 + MRROW) / NPROW
!
!   *********************************************************************
!        Copy blocks of A to WORK(IPA) and perform diagonal scaling
!   *********************************************************************
!
         CALL ZLACPY( 'A',MP,RK*2,WORK(IFTOP),MP,WORK(IWA),MP ) 
         ICOL = MOD( MYROW+K1,NPCOL )

         IF( MYROW .EQ. 0 ) THEN   ! First row, ONLY DL
             IPR = IDL + K1*MP
             CALL ZZLASCL2( MP,RK,ZDWORK(IPR),WORK(IWA),MP )
         ELSE IF( ICOL .EQ. 0 )  THEN                            ! First column,  ONLY DR 
             DO II =1, MP
                IUC=IWA+MP*RK+(II-1)*RK
                CALL ZZLASCL('G',0,0,ONE,ZDWORK(IDR+II-1),RK,1,WORK(IUC),RK,IINFO)
             END DO
         ELSE     ! Requires both DL and DR
             ! DL
             IPR = IDL + K1*MP
             CALL ZZLASCL2( MP,RK,ZDWORK(IPR),WORK(IWA),MP )
             ! DR
             DO II =1, MP
                IUC=IWA+MP*RK+(II-1)*RK
                CALL ZZLASCL('G',0,0,ONE,ZDWORK(IDR+II-1),RK,1,WORK(IUC),RK,IINFO)
             END DO
             ! Print F_{22}
!             IF( MYROW.EQ.1 .AND. MYCOL.EQ.1 ) THEN
!                 WRITE(*,*) 'K1=', K1, 'DL=', ZDWORK(IPR:IPR+MP-1) 
!                 CALL ZGEMM( 'N','N',MP,MP,RK,ONE,WORK(IWA),MP,WORK(IWAVT),RK, &
!                              ZERO,WORK(IWTEMP),MP )
!                 WRITE(*,*) 'F22=', REAL( WORK( IWTEMP:IWTEMP+MP*MP-1 ) )
!             ENDIF
         END IF

!
!   *********************************************************************
!        Perform two low-rank matrix-matrix multiplication
!   *********************************************************************
!
         ! VT*B
         !WRITE(*,*) 'MP,KBP,RK,NQ,KP=',MP,KBP,RK,NQ,KP,IWTEMP 
         CALL ZGEMM( 'No', 'No', RK, NQ, KP, ALPHA, WORK(IWAVT), &
                      RK, B, MAX(1, KBP), CZERO, &
                      WORK(IWTEMP), RK )
         ! U*(VT*B)
         CALL ZGEMM( 'No', 'No', MP, NQ, RK, ONE, WORK(IWA), &
                    MP, WORK(IWTEMP), RK, TBETA, &
                    C, LDC )
!         WRITE(*,*) 'GEMM2: MP,KBP,RK=',MP,KBP,RK 
         TBETA = CONE
!
!        Shift B (WORK(IPB)) to the top (Upwards)
!
         IF( K1.LT.NPROW-1 ) THEN
!
!            Rotate B (WORK(IPB)) based on simultaneous shifting
!
!             WRITE(*,*) 'The communication steps: KBP,NQ=',KBP,NQ
             CALL ZGESD2D( ICTXT, KBP, NQ, B, KBP, &
                          MOD(NPROW+MYROW-1, NPROW), MYCOL )
             KBP = NUMROC( K, KB, MRROW+K1+1, 0, NPROW )
             CALL ZGERV2D( ICTXT, KBP, NQ, B, KBP, & 
                           MOD(MYROW+1, NPROW), MYCOL )
         END IF

      END DO
   
   
      DEALLOCATE( ZDWORK )
!
      RETURN
!
!     End of PDMDFTL
!
      END SUBROUTINE PZMDFTL1
