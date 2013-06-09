/* afwt.c
 * Usage: u = afwt(w, h0, h1,J,sym);
 Adjoint of Fast wavelet transform.
 Written by: M. Salman Asif, Georgia Tech
 December 2007
 */

#include <stdlib.h>
#include "mex.h"

/* these are in afwt_level.c */
void adj_wavelet_one_level(double*, double*, int, double*, double*, int, int);
int checkPowerTwo(int,int);

void adj_wavelet_multi_level(double *u, double *w, int n, double *h0, double *h1, int l, int J, int sym)
{
    int j, nj, m;
    double *u_tmp;
    
    u_tmp = (double*) mxCalloc(n, sizeof(double));
   /* u_tmp = (double*) calloc(n, sizeof(double)); */

    nj = n;
	for (m=0; m < n; m++)
	    u[m] = w[m];
    
	for (j=J-1; j >= 0; j--) {
        nj = n>>j;
        /*if (j == J-1){ /*although i realized now we don't need any if else here because u is already w.*/
        /*    adj_wavelet_one_level(u_tmp, w, nj, h0, h1, l, sym);
        }
        else*/
        adj_wavelet_one_level(u_tmp, u, nj, h0, h1, l, sym);
        
        for (m=0; m < nj; m++)
            u[m] = u_tmp[m];
    }
    mxFree(u_tmp);
}


/* The gateway routine. */
/* u = afwt(w, h0, h1, J, sym); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *u, *w, *h0, *h1;
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
        mexErrMsgTxt("Input w must be real");
    }
    n = mxGetM(prhs[0]);
    J = mxGetScalar(prhs[3]);
    Jmax = checkPowerTwo(n, J);
    if ((J < 0) || (J > Jmax)) {
        mxErrMsgTxt("Input J must be an integer between 0 and log2(n)");
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
    w = mxGetPr(prhs[0]);
    h0 = mxGetPr(prhs[1]);
    h1 = mxGetPr(prhs[2]);
    u = mxGetPr(plhs[0]);
    
  /* Call the C subroutine. */
    adj_wavelet_multi_level(u, w, n, h0, h1, l, J, sym);
    return;
}