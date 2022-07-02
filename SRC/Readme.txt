
1) pzmdftL.f90 uses SVD to compute low-rank approximation; pzmdftL1.f90 uses RRQR to comptue 
   low-rank approximation. 

   Right now DFT codes assume that 1) N=NP*NB; 2) NP=NQ; The process grid must be square. 
