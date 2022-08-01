2022-07-02

1) To add low-rank approximation of Toeplitz matrix.
   Only used when T is stored in BDD format, and T is off-diagonally low-rank.  
   a) Using FFT and convert Toeplitz matrix into Cauchy-like matrix. 
      Since this approach requires to write a complex version of SRRSC, this 
      approach is "abandoned".
      
   b) Another approach is to use RRQR to construct a low-rank approximation. 
      This approach does not require the local matrix to be Toeplitz, and it 
      can be block Toeplitz. Its disadvantage is that it requires construct
      the local matrix explicitly in O(m^2) flops where m is the dimension of
      local matrix. The cost of RRQR is also O(m^2 r), and therefore when $r$
      is \emph{much} smaller than $m$, this approach can be faster.
      The traditional approach requires $O(m^3)$ flops.

      For this we need to modify the dgeqp3 routine, and use dgemqr or dorgqr to construct
      the orthogonal matrix Q. 

2) To re-test structured DFT matrix multiplication


To be done

1) Add PUMMA's original codes (Modified Fox's algorithm)

2) Add 2.5D matrix multiplication routines. The Japanese's scsumma routine
   can be faster than PDGEMM for the case that Multi(A^T, B^T) or multi(A, B^T).

 
2022-07-03

1) Cauchy and Toeplitz matrix multiplication requires MB and NB to be equal. 
   This constraint can be removed. Therefore, they can work for rectangular 
   process grid.   

2) Cauchy and Toeplitz only works for square process grid when using redistribution. 
   Without using redistribution, they works for rectangular process grids. 

3) Work on DFT matrix multiplication. 


4) /home/lsg/MyCodes/Ready_To_Release/YH-ParaNUT/Valid_Test/MyTest 
   There are many routines to refer. 

2022-07-06

1) The Tianhe-2 compiler environment is to use
module load MPI/mpich/3.2.1-gcc-4.8.5-dynamic
source /BIGDATA1/app/intelcompiler/18.0.0/compilers_and_libraries_2018/linux/bin/compilervars.sh intel64



