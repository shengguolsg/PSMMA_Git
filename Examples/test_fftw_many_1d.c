
#include <stdlib.h>
#include <fftw3-mpi.h>

int main(int argc, char **argv){
    const ptrdiff_t N[1]={4096};
    ptrdiff_t howmany=4096;
    fftw_plan plan;
    fftw_complex *in;
    fftw_complex *out; 
    ptrdiff_t alloc_local, local_n0, local_0_start, block;
    int me,np,err;
    int rank; 
    double startwtime, endwtime;

    rank =1; 

    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD,&np);

    block=FFTW_MPI_DEFAULT_BLOCK;
    alloc_local = fftw_mpi_local_size_many(rank,N,howmany,block, MPI_COMM_WORLD,&local_n0, &local_0_start);
    in  = fftw_alloc_complex(alloc_local);
    //out = fftw_alloc_complex(alloc_local); 

    if (me==0){
        printf("Test is running at %d processes \n",np);
        printf("Before plan calculation, alloc_local=%ld, block=%ld\n",local_n0,block );
    }

    /* initialize data to some function my_function(x,y) */
    for (int i = 0; i < local_n0; ++i) 
    {
        in[i][0] =rand() / (double)RAND_MAX;
        in[i][1] =rand() / (double)RAND_MAX;
    }

    if(me ==0) {
        startwtime = MPI_Wtime();
    }
    plan = fftw_mpi_plan_many_dft(rank,N,howmany,block,block, in, in, MPI_COMM_WORLD, 
                      FFTW_FORWARD, FFTW_ESTIMATE); // FFTW_ESTIMATE); // FFTW_PATIENT);
    if(me ==0) {
        endwtime = MPI_Wtime();
        printf("wall clock time = %f\n", endwtime-startwtime);
    }

    if (me==0){
        printf("After plan calculation\n" );
    }

    fftw_destroy_plan(plan);
    MPI_Finalize();
}