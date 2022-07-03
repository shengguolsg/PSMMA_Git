    #include <fftw3-mpi.h>
    # include <stdlib.h>
    # include <stdio.h>
    #include <sys/stat.h>
    #include <fcntl.h>
    # include <time.h>
    #include <math.h>

    int main(int argc, char **argv)
    {
       const ptrdiff_t N0 = 4096;
       //const ptrdiff_t N0 = 4194304 ;
       //const ptrdiff_t N0=16777216;
       //const ptrdiff_t N0 = 8388608;
       fftw_plan planForw,planBack;
       fftw_complex *data,*dataOut,*data2;
       ptrdiff_t alloc_local, local_ni, local_i_start, i, j,local_no, local_o_start;
       int index,size;
       double startwtime, endwtime;
       MPI_Init(&argc, &argv);
       fftw_mpi_init();
       MPI_Comm_rank(MPI_COMM_WORLD,&index);
       MPI_Comm_size(MPI_COMM_WORLD,&size);

       /* get local data size and allocate */
       alloc_local = fftw_mpi_local_size_1d(N0, MPI_COMM_WORLD,FFTW_FORWARD, FFTW_ESTIMATE,
                                                  &local_ni, &local_i_start,&local_no, &local_o_start);
       data = fftw_alloc_complex(alloc_local);
       dataOut = fftw_alloc_complex(alloc_local);
       data2 = fftw_alloc_complex(alloc_local);
             /* create plan  */
       planForw = fftw_mpi_plan_dft_1d(N0, data, data2, MPI_COMM_WORLD,
                                         FFTW_FORWARD, FFTW_ESTIMATE);
       planBack = fftw_mpi_plan_dft_1d(N0, data2, dataOut, MPI_COMM_WORLD,
                                         FFTW_BACKWARD, FFTW_ESTIMATE);
       /* initialize data to some function my_function(x,y) */
       for (i = 0; i < local_ni; ++i) 
       {
        data[i][0] =rand() / (double)RAND_MAX;
        data[i][1] =rand() / (double)RAND_MAX;
       }
       if(index==0){
         startwtime = MPI_Wtime();
       }
        for (i=0; i<N0; i++) {
            fftw_execute(planForw);
        }
        if(index==0){
            endwtime = MPI_Wtime();
            printf("wall clock time = %f\n",
                           endwtime-startwtime);
        }
        if(index==0){
         startwtime = MPI_Wtime();
        }
        fftw_execute(planBack);
        if(index==0){
            endwtime = MPI_Wtime();
            printf("wall clock time = %f\n",
                           endwtime-startwtime);
        }

        fftw_destroy_plan(planForw);
        fftw_destroy_plan(planBack);
        MPI_Finalize();
}