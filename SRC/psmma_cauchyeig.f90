!
   SUBROUTINE PSMMA_CAUCHYEIG( N,K,ID,ALPHA,A,IA,JA,DESCA,F,D,U,V,DIFL,&
                           DIFR,BETA,Q,IQ,JQ,DESCQ,WORK,Redist )
!
!      USE MPI
      !USE Scal_Aux
      USE scal_aux
      implicit none
    
      include 'mpif.h'
!
!  -- PSMMA computational routine (version 1.0.0) --
!
!     .. Scalar Arguments ..
      INTEGER            N, K, ID, IA, JA, IQ, JQ
      DOUBLE PRECISION   ALPHA, BETA
      LOGICAL            Redist
!     ..
!     .. Array Arguments ..
      INTEGER          ::   DESCA( * ), DESCQ( * )
      DOUBLE PRECISION ::   A(*),Q(*),F(*),D(*),U(*),V(*),DIFL(*),&
                            DIFR(*), WORK(*)
!     ..
!
!  Purpose
!  =======
!
!  PSMMA_CauchyEig computes the multiplication of a matrix A with an HSS matrix, 
!  Q = alpha* A * H + beta * Q, where A is N-by-K, and H is K-by-K,
!  H(i,j) = u(i)*z(j) / ( d(i)-lambda(j) ).
!  It uses structured Cauchy-like parallel matrix-matrix multiplication  
!  algorithms, PSMMA, which is introduced in TPDS 2020, by
!  Xia Liao, Shengguo Li, Yutong Lu and Jose E. Roman.   
!  
!  The matrix A is initially distributed in the block-cyclic form (BCDD), and the 
!  generators of H are stored globally, all the entries can be
!  constructed locally. Based on the parameter 'Redist', this routine decides
!  whether redistribute matrix A from BCDD to BDD or not. If it requires to redistribute
!  matrix A, it calls PDGEMR2D and then call PSCAUCHYEIG_COMPUTE to perform PSMMA. 
!
!  PUMMA works for block-cyclic from but difficult to implement. While, PSMMA further takes
!  advantage with the off-diagonal property and several advantages. 
! 
!  We assume that A is N-by-K, and H is K-by-K. 
! 
! ===========================
! Written by S.G. Li, NUDT, Oct. 2019
! =======================================================
!
!     .. Parameters ..
!
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_, &
                         MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                         CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                         RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            ICTXT, IIQ, IQCOL,IQROW, JJQ, info, &
                         LDAT,LDQT,MYCOL,MYROW,NB,NP, NPCOL, NPROW,&
                         MBT, NBT, NQT, NPT, NQ, MB, IAT, IQT, IWEND2,&
                         LDB, LDA, LDQ
      REAL(8)            time1, time2
!     ..
!     .. Local Arrays ..
      INTEGER       ::   DESCAT(DLEN_),DESCQT(DLEN_)
!     ..
!     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DESCINIT, INFOG1L,   &
                         INFOG2L
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_* &
          RSRC_.LT.0 ) RETURN
!
!     Test the input parameters.
!
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!
      ICTXT = DESCQ( CTXT_ )
      LDA = DESCA( LLD_ )
      LDQ = DESCQ( LLD_ )
      MB    = DESCQ( MB_ )
      NB    = DESCQ( NB_ )
       
      CALL INFOG2L( ID, ID, DESCQ, NPROW, NPCOL, MYROW, MYCOL, &
                    IIQ, JJQ, IQROW, IQCOL )
!                    
      NP = NUMROC( N, DESCQ( MB_ ), MYROW, IQROW, NPROW )
      NQ = NUMROC( K, DESCQ( NB_ ), MYCOL, IQCOL, NPCOL )
!      
!      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!            WRITE(*,*) 'In psmma_cauchyeig, IQROW,IQCOL= ', IQROW,IQCOL
!      ENDIF

      IF( .NOT. Redist ) THEN
            IQT = 1
            LDQT = NP
            IWEND2 = IQT+NP*NQ
            CALL DESCINIT( DESCQT, N, K, MB, NB, IQROW, IQCOL, ICTXT, LDQT, &
                           INFO )

            !CALL PDtestOrth( N,K,A,1,1,DESCA )

            !time1 = MPI_Wtime()
            CALL PSCAUCHYEIG_COMPUTE( N,K,MB,NB,ALPHA,A,LDA,DESCA,F,D,U,V,DIFL, &
                              DIFR,BETA,WORK(IQT),LDQT,IQROW,IQCOL,WORK(IWEND2) )
            !time2 = MPI_Wtime()
            !IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            !      write(*,*) 'The PSMMA without redistribution costs', time2-time1
            !END IF
            !CALL PDtestOrth( N,K,WORK(IQT),1,1,DESCQT )

            ! Copy back the computed Q to Q
            CALL PDLACPY( 'A',N,K,WORK(IQT),1,1,DESCQT,Q,IQ,JQ,DESCQ )

            GOTO 90  
      END IF
!
!     Construct the descriptor for the reshaped matrix A and Q
      MBT = N / NPROW 
      NBT = K / NPCOL
      IF( MBT*NPROW .NE. N )  MBT = MBT+1
      IF( NBT*NPCOL .NE. K )  NBT = NBT+1 
      NPT = NUMROC( N, MBT, MYROW, 0, NPROW )
      NQT = NUMROC( K, NBT, MYCOL, 0, NPCOL )

      !IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
      !  WRITE(*,*) 'In PDHSSEVC1,MB,NB,MBT,NBT,NPT,NQT=',MB,NB,MBT,NBT,&
      !             NPT,NQT, 'MYROW,MYCOL=',MYROW,MYCOL 
      !ENDIF

      !LDAT = MAX( NPT,MBT )   ! The locate leading dimension of AT
      LDAT = NPT
      LDQT = LDAT 
      LDB  = MAX( NQT,NBT )
!
!     Here we assume N == K
      CALL DESCINIT( DESCAT, N, K, MBT, NBT, 0, 0, ICTXT, LDAT, &
                     INFO )
      CALL DESCINIT( DESCQT, N, K, MBT, NBT, 0, 0, ICTXT, LDQT, &
                     INFO )
!
!      ALLOCATE( AT(NPT*NQT), QT(NPT*NQT), stat=ierr )
!      ALLOCATE( AT(MBT*NBT), QT(MBT*NBT), stat=ierr )
!      IF( ierr .ne. 0 ) THEN
!         WRITE(*,*) 'Allocation failed in pdhssevc1, and return'
!         RETURN
!      END IF

       IAT = 1
       IQT = 1+ NPT*NQT
       IWEND2= IQT+NPT*NQT ! IWEND2 must left at 3*NPT*NQT
    
      ! Redistribute matrix A to blocked form.  No matter A(IA,JA) is located in which
      ! process, we always redistribute it AT(1,1) to process (0,0). Then, we can use
      ! SPUMMA via generators. Otherwise, we can only use SPUMMA at the top level of DC. 
!$      time1 = MPI_Wtime()
      CALL PDGEMR2D( N,K,A,1,1,DESCA,WORK(IAT),1,1,DESCAT,ICTXT ) 
!$      time2 = MPI_Wtime()
!$      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!$           write(*,*) 'The first redistribution costs', time2-time1
!$      END IF 

!      CALL PDtestOrth( N,K,WORK(IAT),1,1,DESCAT )
!      CALL PDtestOrth( N,K,A,1,1,DESCA )
!      IF(MYROW.EQ.0 .AND. MYCOL.EQ.0) THEN
!          WRITE(*,*) 'We are testing Orthogonality!', N, K
!          CALL TestOrthCauchylike( F,U,V,DIFL,DIFR,K,K )
!      END IF

      ! Call structured PSMMA matrix-matrix multiplication
      ! If using Redistribution, all the matrices are starting from (0,0)
      IQROW = 0
      IQCOL = 0
!$      time1 = MPI_Wtime()
      CALL PSCAUCHYEIG_COMPUTE( N,K,MBT,NBT,ALPHA,WORK(IAT),LDAT,DESCAT,F,D,U,V,DIFL, &
                              DIFR,BETA,WORK(IQT),LDQT,IQROW,IQCOL,WORK(IWEND2) )  
!$      time2 = MPI_Wtime()
!$      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!$           write(*,*) 'The PSMMA costs', time2-time1, MYROW,MYCOL
!$      END IF 
!      CALL PDtestOrth( N,K,WORK(IQT),1,1,DESCQT )

!      CALL BLACS_BARRIER(ictxt,'All')
      ! Redistribute the matrix Q back
!$      time1 = MPI_Wtime()
      CALL PDGEMR2D( N,K,WORK(IQT),1,1,DESCQT,Q,IQ,JQ,DESCQ,ICTXT ) 
!$      time2 = MPI_Wtime()
!$      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!$           write(*,*) 'The second redistribution costs', time2-time1
!$      END IF
       
      ! Clean-up
      !DEALLOCATE( AT,QT )

  90    RETURN

      END SUBROUTINE PSMMA_CAUCHYEIG
