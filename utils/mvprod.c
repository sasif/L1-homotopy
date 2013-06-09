/* mvprod.c 
 * Usage: lhs = mvprod(A,rhs,Gamma,flag);
 * y = A_Gamma * x_Gamma
 * x_Gamma = A_Gamma' * y
 Matrix vector multiplication using indices of the matrix.
 Written by: Salman Asif, Georgia Tech
*/

#include <stdlib.h>
#include "mex.h"

/* The gateway routine. */
/* lhs = mvprod(A, rhs, Gamma, flag); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *rhs, *lhs, *A, *Gamma, rhs_v;
  int n, m, flag, T, i, j, ij, gj;
  
  /* Check for the proper number of arguments. */
  if (nrhs != 4) {
    mexErrMsgTxt("Exactly four inputs required");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }
  if (mxGetN(prhs[1]) > 1) {
    mexErrMsgTxt("Input x must be a column vector");
  }
  if (mxIsComplex(prhs[1])) {
    mexErrMsgTxt("Input x must be real");    
  }
  
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  T = mxGetM(prhs[2]);
  flag = mxGetScalar(prhs[3]);

  if (mxGetN(prhs[2]) > mxGetM(prhs[2])){
      mexErrMsgTxt("Gamma must be a column vector");
  }
  if (flag > 1) {
    mexErrMsgTxt("Flag must be 0 or 1");
  }
  if (flag == 0) {
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    if (mxGetN(prhs[0]) != mxGetM(prhs[1])) {
        mexErrMsgTxt("Dimensions of A and x do not match");
    }
  }
  if (flag == 1) {
      plhs[0] = mxCreateDoubleMatrix(T, 1, mxREAL);
      if (mxGetM(prhs[0]) != mxGetM(prhs[1])) {
          mexErrMsgTxt("Dimensions of A' and y do not match");
      }
  }

  
  /* Assign pointers to each input and output. */
  A = mxGetPr(prhs[0]);
  rhs = mxGetPr(prhs[1]);
  Gamma = mxGetPr(prhs[2]);
  lhs = mxGetPr(plhs[0]);
  
  /* Do matrix-vector multiplication */ 
  if (flag == 0){
    for (i=0; i < m; i++) {
        lhs[i] = 0;
    }
    for (j=0; j < T; j++){
        gj = Gamma[j]-1;
        rhs_v = rhs[gj];
        gj = gj*m;
        for (i=0; i < m; i++){
            ij = gj+i;
            lhs[i] += A[ij]*rhs_v;
            /* mexPrintf("\t lhs[j] = %f, A[ij]=%f, rhs[i]=%f. \n", lhs[j],A[ij],rhs[i]); */
        }
    }
  }
  else{
    for (j=0; j < T; j++){
        lhs[j] = 0;
        gj = (Gamma[j]-1)*m;
        for (i=0; i < m; i++) {
            ij = gj+i;
            lhs[j] += A[ij]*rhs[i];
            /* mexPrintf("\t lhs[j] = %f, A[ij]=%f, rhs[i]=%f. \n", lhs[j],A[ij],rhs[i]); */
        }
    }
  }  
}