
program psmma_testhankel

  implicit none
  include  'mpif.h'
  
  integer, parameter :: BLACSCTXTSIZE=9

  integer :: descA(BLACSCTXTSIZE), descXB(BLACSCTXTSIZE)
  integer :: n
  integer :: nb, ldta
  integer :: locr, locc
  integer :: i, j, ii, jj
  integer :: ierr, dummy
  integer :: myid, np
  integer :: myrow, mycol, nprow, npcol
  integer :: ictxt, lwork, NTAC
  logical :: REDIST
  double precision :: xnorm, xnorm0, time1,time2
!
  double precision, allocatable :: A(:), B(:), X(:), X1(:), WORK(:)
  double precision, allocatable :: D(:),TA(:)
  double precision,parameter :: one=1.0D+0, zero=0.0D+0
  
  integer, parameter :: IONE = 1, INONE=-1, IZERO=0

  integer, external :: indxl2g, numroc


  !n= 8192  ! Size of the problem
  n= 16384  ! Size of the problem
  !n= 32768   ! Size of the problem
  !n = 768 

  allocate( D(2*N-1),TA(2*N-1) )
  ldta = 2*N -1

  call random_number(D)
  TA = D

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

  ! Initialize the BLACS grid
  !nprow=floor(sqrt(real(np)))
  !npcol=np/nprow

  do nprow = NINT(SQRT(REAL(np))),2,-1
      if(mod(np,nprow) == 0 ) exit
  enddo
  ! at the end of the above loop, nprocs is always divisible by np_cols
   npcol = np/nprow

!  do npcol = NINT(SQRT(REAL(np))),2,-1
!      if(mod(np,npcol) == 0 ) exit
!  enddo
!  ! at the end of the above loop, nprocs is always divisible by np_cols
!   nprow = np/npcol

   call blacs_get(IZERO,IZERO,ictxt)
   call blacs_gridinit(ictxt,'R',nprow,npcol)
   call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)
  
   !nb = n/nprow 
   nb= 64   ! Blocksize of the 2D block-cyclic distribution
   if(myrow.eq.0 .and. mycol.eq.0) write(*,*) "n,nprow,npcol,nb=", n, nprow,npcol,nb

  ! A is a dense n x n distributed Toeplitz matrix
  if(myid<nprow*npcol) then
    locr=numroc(n,nb,myrow,IZERO,nprow)
    locc=numroc(n,nb,mycol,IZERO,npcol)
    allocate(A(locr*locc))
    A=0
    dummy=max(1,locr)
    call descinit(descA,n,n,nb,nb,IZERO,IZERO,ictxt,dummy,ierr)

!    write(*,*) 'myid=', myid, myrow, mycol, locr

    do j=1,locc
      jj=indxl2g(j,nb,mycol,IZERO,npcol)
      do i=1,locr
        ii=indxl2g(i,nb,myrow,IZERO,nprow)
        ! Toeplitz matrix
          A(locr*(j-1)+i)= D( ii+jj-1 )
      end do
    end do
  else
    call descset(descA,n,n,nb,nb,IZERO,IZERO,INONE,IONE)
  end if

  ! Set the random matrix
  if(myid<nprow*npcol) then
    locr=numroc(n,nb,myrow,IZERO,nprow)
    locc=numroc(n,nb,mycol,IZERO,npcol)
    allocate( B(locr*locc) )
    dummy=max(1,locr)
    call descinit(descXB,n,n,nb,nb,IZERO,IZERO,ictxt,dummy,ierr)
!    B=1D0
    allocate( X(locr*locc) )
    allocate( X1(locr*locc) )
  else
    call descset(descXB,n,n,nb,nb,IZERO,IZERO,INONE,IONE)
  end if

!  write(*,*) "Before call pdmma", N, nb, nprow, npcol
!  write(*,*) "myrow, mycol",myrow, mycol, locr, locc, ictxt

  call random_number(B)
  ! /* Original matrix-matrix mulitplication */
  ! X1 = B*A, A is Toeplitz
  time1 = MPI_Wtime()
  call PDGEMM( 'N','N',n,n,n,one,B,1,1,descA,A,1,1,descA,zero,X1,1,1,descA )
  time1 = MPI_Wtime()-time1
  IF(myid .eq. 0) WRITE(*,*) "PDGEMM finishes and costs ", time1

  time1 = MPI_Wtime()
  call PDGEMM( 'N','N',n,n,n,one,B,1,1,descA,A,1,1,descA,zero,X1,1,1,descA )
  time1 = MPI_Wtime()-time1
  IF(myid .eq. 0) WRITE(*,*) "PDGEMM finishes and costs ", time1

!  WRITE(*,*) "PDGEMM finishes"

  !REDIST = .true.
  REDIST = .false.
  LWORK = 6*locr*locc
  allocate( WORK(LWORK) )
  ! /* Matrix-matrix multiplication */
  time2 = MPI_Wtime()
  call psmma_hankel( n,n,one,B,1,1,desca,TA,ldta,zero,X,1,1,desca,WORK,REDIST )
  time2 = MPI_Wtime()-time2
  IF(myid .eq. 0) WRITE(*,*) "PDMMMA finishes and costs ", time2, time1/time2

  time2 = MPI_Wtime()
  call psmma_hankel( n,n,one,B,1,1,desca,TA,ldta,zero,X,1,1,desca,WORK,REDIST )
  time2 = MPI_Wtime()-time2
  IF(myid .eq. 0) WRITE(*,*) "PDMMMA finishes and costs ", time2, time1/time2

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
  if(myid .eq. 0 ) write(*,*) 'Error norm is ', xnorm/xnorm0, xnorm0
  
  ! Clean-up
  if( myid < nprow*npcol ) then
    deallocate(A,TA,D)
    deallocate(B)
    deallocate(X)
    deallocate(X1)
    deallocate(WORK)
  end if

  ! The end
  call MPI_Finalize(ierr)

end program psmma_testhankel
