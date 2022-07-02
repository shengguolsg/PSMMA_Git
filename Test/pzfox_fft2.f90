
  program pzfox_fft2

!  use MPI

  implicit none
  include  'mpif.h'
  
  integer, parameter :: BLACSCTXTSIZE=9

  integer :: descA(BLACSCTXTSIZE), descXB(BLACSCTXTSIZE)
  integer :: n, i, j, jj, ii
  integer :: nb
  integer :: locr, locc
  integer :: ierr, dummy
  integer :: myid, np
  integer :: myrow, mycol, nprow, npcol
  integer :: ictxt, lwork, LWORK1
  double precision :: xnorm, xnorm0, time, time1
  double precision, allocatable :: RWORK(:)
  complex*16, allocatable :: A(:), B(:), WORK(:),X(:),X1(:)
  
  double precision,parameter :: one=1.0D+0, zero=0.0D+0
  
  integer, parameter :: IONE = 1, INONE=-1, IZERO=0
  !DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793E+0
  DOUBLE PRECISION :: PI 
  complex*16, parameter :: CI=( zero, one )

  complex*16  Wn

  integer, external :: indxl2g, numroc

  EXTERNAL  PDGEMM

  COMMON /CONTEXT / ICTXT

  n  = 1024   ! Size of the problem
  nb = 512   ! Blocksize of the 2D block-cyclic distribution
  PI  = 4.0E+0 * ATAN( ONE )
  Wn  = EXP( -2.0E+0*PI*CI/N )

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

  ! Initialize the BLACS grid
  do npcol = NINT(SQRT(REAL(np))),2,-1
      if(mod(np,npcol) == 0 ) exit
  enddo
  ! at the end of the above loop, nprocs is always divisible by np_cols

   nprow = np/npcol

   call blacs_get( IZERO,IZERO,ictxt )
   call blacs_gridinit( ictxt,'R',nprow,npcol )
   call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)

   IF( myid .eq. 0 ) write(*,*) 'N, NB, nprow, npcol=', N,NB,nprow,npcol

  ! A is a dense n x n distributed FFT matrix
  if(myid<nprow*npcol) then
    locr=numroc(n,nb,myrow,IZERO,nprow)
    locc=numroc(n,nb,mycol,IZERO,npcol)
    allocate(A(locr*locc))
    A=0
    dummy=max(1,locr)
    call descinit(descA,n,n,nb,nb,IZERO,IZERO,ictxt,dummy,ierr)

    do i=1,locr
      do j=1,locc
        ii=indxl2g(i,nb,myrow,IZERO,nprow)
        jj=indxl2g(j,nb,mycol,IZERO,npcol)
        ! DFT matrix 
          A(locr*(j-1)+i)= Wn**( (ii-1)*(jj-1) )
      end do
    end do
  else
    call descset(descA,n,n,nb,nb,IZERO,IZERO,INONE,IONE)
  end if

!  IF( myid .eq. 0 ) write(*,*) 'A=', REAL ( A(1:locr*locc) )

  ! Set the random matrix
  if(myid<nprow*npcol) then
    locr=numroc(n,nb,myrow,IZERO,nprow)
    locc=numroc(n,nb,mycol,IZERO,npcol)
    allocate( B(locr*locc) )
    dummy=max(1,locr)
    call descinit(descXB,n,n,nb,nb,IZERO,IZERO,ictxt,dummy,ierr)
    B=1.0E+0
    allocate( X(locr*locc) )
    allocate( X1(locr*locc) )
  else
    call descset(descXB,n,n,nb,nb,IZERO,IZERO,INONE,IONE)
  end if

  ! /* Original matrix-matrix mulitplication */
  time1 = MPI_Wtime()
  call PZGEMM( 'N','N',n,n,n,one,A,1,1,descA,B,1,1,descA,zero,X1,1,1,descA )
  time = MPI_Wtime()-time1
  IF( myid .eq. 0 ) write(*,*) 'PZGEMM costs ', time


  LWORK = (5*NB+6)*NB
  LWORK1= 7*NB
  allocate( WORK(LWORK), RWORK(LWORK1) )
  IF(myid .eq. 0 ) WRITE(*,*) 'LWORK=', LWORK

  ! /* Matrix-matrix multiplication */
  time1 = MPI_Wtime()
  call pzmdftL( n,n,nb,nb,one,B,locr,zero,X,locr,0,0,WORK,RWORK )
  time = MPI_Wtime()-time1
  IF( myid .eq. 0 ) write(*,*) 'PZMDFTL costs ', time

  ! Check the accuracy of computed results 
  if(myid < nprow*npcol) then
     call pdnrm2( n, xnorm0, X, 1, 1, descXB, 1)
     locr=numroc(n,nb,myrow,IZERO,nprow)
     locc=numroc(n,nb,mycol,IZERO,npcol)
     do i =1, locr*locc
        X1(i) = X1(i) -X(i)
     end do
     call pdnrm2( n, xnorm, X1, 1, 1, descXB, 1)
  end if
  if(myid .eq. 0 ) write(*,*) 'Error norm is ', xnorm/xnorm0, nprow, npcol

  
  ! Clean-up
  if( myid < nprow*npcol ) then
    deallocate(A)
    deallocate(B)
    deallocate(X)
    deallocate(X1)
    deallocate(WORK)
    deallocate(RWORK)
  end if

  ! The end
  call MPI_Finalize(ierr)

end program pzfox_fft2
