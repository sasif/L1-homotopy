/* fwt.c 
 * Usage: w = fwt(x, h0, h1, J, sym);
 Fast wavelet transform.
 Written by: Justin Romberg, Georgia Tech
*/

#include <stdlib.h>
#include "mex.h"

/* these are in fwt_level.c */
void wavelet_one_level(double*, double*, int, double*, double*, int, int);
int checkPowerTwo(int,int);



void wavelet_multi_level(double *w, double *x, int n, double *h0, double *h1, int l, int J, int sym)
{ 
  int j, nj, m;
  double *tmp;
  
  tmp = (double*) mxCalloc(n/2, sizeof(double));
  
  nj = n;
  for (j=0; j < J; j++) {
    
    if (j == 0)
      wavelet_one_level(w, x, nj, h0, h1, l, sym);
    else
      wavelet_one_level(w, tmp, nj, h0, h1, l, sym);
    
    nj >>= 1;
    
    if (j < J-1) 
      for (m=0; m < nj; m++) 
        tmp[m] = w[m];
      
    
  }
   
  mxFree(tmp);
}


/* The gateway routine. */
/* w = fwt(x, h0, h1, J, sym); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *x, *w, *h0, *h1;
  int n, l, l1, sym, J, Jmax;
  
  /* Check for the proper number of arguments. */
  if (nrhs != 5) {
    mexErrMsgTxt("Exactly five inputs required");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }
  
  if (mxGetN(prhs[0]) > 1) {
    mexErrMsgTxt("Input x must be a column vector");
  }
  if (mxIsComplex(prhs[0])) {
    mexErrMsgTxt("Input x must be real");    
  }
  n = mxGetM(prhs[0]);
  J = mxGetScalar(prhs[3]);
  Jmax = checkPowerTwo(n, J);
  if ((J < 0) || (J > Jmax)) {
    mexErrMsgTxt("Input J must be an integer between 0 and log2(n)");
  }
  l = (mxGetM(prhs[1]) > mxGetN(prhs[1])) ? mxGetM(prhs[1]) : mxGetN(prhs[1]);
  l1 = (mxGetM(prhs[2]) > mxGetN(prhs[2])) ? mxGetM(prhs[2]) : mxGetN(prhs[2]);
  if (l != l1) {
    mexErrMsgTxt("Filters must be the same length");
  }
  sym = mxGetScalar(prhs[4]);
  if (sym > 2) {
    mexErrMsgTxt("Symmetry flag must be 0, 1, or 2");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  
  /* Assign pointers to each input and output. */
  x = mxGetPr(prhs[0]);
  h0 = mxGetPr(prhs[1]);
  h1 = mxGetPr(prhs[2]);
  w = mxGetPr(plhs[0]);
  
  /* Call the C subroutine. */
  wavelet_multi_level(w, x, n, h0, h1, l, J, sym);
  return;
}