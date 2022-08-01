Module Scal_Aux

  use CauchyHssEig
  use BasicMM
!
  implicit none

Contains

!!!!!!
   Subroutine PDtestOrth( M,N, A,IA,JA,DESCA)
!
!    use MPI
!
    include 'mpif.h'
!
!  This routine tests the orthogonality of matrix A, which is a distributed !  matrix. 
!
   INTEGER  ::  M,N,IA,JA
!
   INTEGER  ::  DESCA(*)
   DOUBLE PRECISION :: A(*)
! ..
! .. Parameters ..
   INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_, &
                      MB_, NB_, RSRC_, CSRC_, LLD_
   PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                      CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                      RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
   DOUBLE PRECISION   ZERO, ONE
   PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
! ..
! .. Local Scalars .. 
   INTEGER            ICTXT, info, LDB, MYCOL,MYROW,&
                      MB,NB,NPCOL, NPROW, IQROW, IQCOL,mpierr, &
                      NP,NQ,IIA,JJA
   REAL(8)            err,errmax
!  ..
!  .. Local Arrays ..
   INTEGER    ::   DESCB(DLEN_)
   DOUBLE PRECISION, ALLOCATABLE ::  B(:,:), D(:,:)
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!
!  
   CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!
   ICTXT = DESCA( CTXT_ )
   MB    = DESCA( MB_ )
   NB    = DESCA( NB_ )
!
   CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, &
                 IIA, JJA, IQROW, IQCOL )

   NP = NUMROC( N, DESCA( NB_ ), MYROW, IQROW, NPROW )
   NQ = NUMROC( N, DESCA( NB_ ), MYCOL, IQCOL, NPCOL )
   
   !LDB = MAX(NB,NP)
   LDB = NP

!   WRITE(*,*) 'NB,NQ=',NB,NQ,MYROW,MYCOL
!  B is the transpose of A. Here, we assume NPROW == NPCOL
   CALL DESCINIT( DESCB, N, N, NB, NB, IQROW, IQCOL, ICTXT, LDB, &
                  INFO )
!
   ALLOCATE( B(LDB,NQ),D(NP,NQ) )

!  B = AT*A
   CALL PDGEMM( 'T','N',N,N,M,ONE,A,1,1,DESCA,A,1,1,DESCA,&
                ZERO,B,1,1,DESCB  )  
   D = ZERO
   CALL PDLASET( 'A',N,N,ZERO,ONE,D,1,1,DESCB ) 

   D(:,:) = B(:,:)-D(:,:)
  
   err = maxval(abs(D))
   CALL MPI_ALLREDUCE(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD, &
                         mpierr)
   IF(MYROW.EQ.0 .AND. MYCOL.EQ.0 ) print *,'Error Orthogonality:',errmax

   DEALLOCATE( B,D )

  End Subroutine PDtestOrth
        
!!!!!!
  SUBROUTINE TestOrthCauchylike( F,U,V,DIFL,DIFR,M,N )
!
! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  ::  M, N
! .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN)  :: F(*),U(*),V(*),DIFL(*),DIFR(*)
!
! Purpose 
! =======
! It computes an Cauchy matrix A with dimension M-by-N, 
! A(i,j) = U(i)*V(j) / ( F(i) - D(j) ). The entries of F are the old eigenvalues,
! and D contains the updated eigenvalues. The computed eigenvector
! matrix satisfies F*A-A*D = U*V^T. 
!
! Parameters
! F   (input) DOUBLE PRECISION array, DIMENSION ( LN )
!     The original old eigenvalues which are used directly, the row generators.
!     Its entries are referenced by PL.
!
! D   (input)  DOUBLE PRECISION array, DIMENSION ( N )
!     D contains the updated eigenvalues for current HSS block,
!     D corresponds to PU, and their relation is D(i) = D0( PU(i) );
!
! U   (input) DOUBLE PRECISION array, DIMENSION ( M )
!     Row Generator of Cauchy matrix, which are referenced directly.
! 
! V   (input) DOUBLE PRECISON array, DIMENSION ( N )
!     Column Generator of Cauchy matrix, which are referenced directly.
!
! DIFL (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFL(i) = D0(i)-F(i), the order of which will not change and its elements
!      are also referenced by using PL or PU. Positive
!
! DIFR (input) DOULBE PRECISION array, DIMENSION ( LN )
!      DIFR(i) = D0(i)-F(i+1), the order of which will not change and its elements
!      are also referenced by using PL or PU. Negative
!
! PL  (input) INTEGER array, DIMENSION ( N )
!     The row permutations of generator.  It relates with F.
!
! PU  (input) INTEGER array, DIMENSION ( N )
!     The column permutations of generator. It relates with D.
!
!  M  (input) INTEGER, the row dimension of A
!
!  N  (input) INTEGER, the column dimension of A
! 
! nflops (output) INTEGER 
!        Floating point operations of generating A.
!
! =============
!  Written by S.-G. Li, on April 17th, 2013
!  for tridiagonal eigenvalue problems
! =======================================
!
! .. Local Scalars ..
    INTEGER           i,j,PFII,PDJJ
    DOUBLE PRECISION  FJ, VJ
!
! .. Local Arrays ..
    DOUBLE PRECISION, ALLOCATABLE :: A(:,:)

    Allocate( A(M,N) ) 

!$OMP PARALLEL PRIVATE(PDJJ,FJ,VJ,PFII)
!$OMP DO SCHEDULE(dynamic)
    DO j = 1,N
       PDJJ = j
       FJ = F(PDJJ)
       VJ = V(j)
       DO i = 1, M
          PFII = i
          IF( PFII .le. PDJJ ) THEN
             A(i,j) = U(i)*VJ / ( (F(PFII)-FJ )-DIFL(PDJJ) )
          ELSE
             A(i,j) = U(i)*VJ / ( ( F(PFII)-F(PDJJ+1) )-DIFR(PDJJ) )
          END IF
       END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CALL TestOrth( A,M,N ) 

  DEALLOCATE( A )
    
  END SUBROUTINE TestOrthCauchylike


End Module Scal_Aux 
