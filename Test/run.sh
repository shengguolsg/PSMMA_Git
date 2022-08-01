OMP_NUM_THREADS=1 yhrun -N 171 -n 4096 ./test2  2>&1 | tee -a Results_4800.dft
#OMP_NUM_THREADS=1 yhrun -N 96 -n 2304 ./test2  2>&1 | tee -a Results_4800.dft
#OMP_NUM_THREADS=1 yhrun -N 43 -n 1024 ./test2  2>&1 | tee -a Results_2400.dft
#OMP_NUM_THREADS=1 yhrun -N 11 -n 256 ./test2  2>&1 | tee -a Results_2400.dft
#OMP_NUM_THREADS=1 yhrun -N 3 -n 64 ./test2  2>&1 | tee -a Results_2400.dft
