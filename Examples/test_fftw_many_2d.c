#include <mpi.h>
#include <fftw3-mpi.h>

int main(int argc, char **argv){
    const ptrdiff_t N[2]={448,352};
    ptrdiff_t Nz=354;
    fftw_plan plan;
    fftw_complex *data;
    ptrdiff_t alloc_local, local_n0, local_0_start, block;
    int me,np,err;

    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD,&np);

    block=FFTW_MPI_DEFAULT_BLOCK;
    alloc_local = fftw_mpi_local_size_many(2,N,Nz,block, MPI_COMM_WORLD,&local_n0, &local_0_start);
    data = fftw_alloc_complex(alloc_local);

    if (me==0){
        printf("Test is running at %d processes \n",np);
        printf("Before plan calculation\n" );
    }
    plan = fftw_mpi_plan_many_dft(2,N,Nz,block,block, data, data, MPI_COMM_WORLD, 
                      FFTW_FORWARD, FFTW_PATIENT);
    if (me==0){
        printf("After plan calculation\n" );
    }

    fftw_destroy_plan(plan);
    MPI_Finalize();
}