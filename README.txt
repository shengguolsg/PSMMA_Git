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

      For this we need to modify the dgeqp3 routine, and use dgemqr to construct
      the orthogonal matrix Q

2) To re-test structured DFT matrix multiplication


To be done

1) Add PUMMA's original codes (Modified Fox's algorithm)

2) Add 2.5D matrix multiplication routines. The Japanese's scsumma routine
   can be faster than PDGEMM for the case that Multi(A^T, B^T) or multi(A, B^T).

3) 
