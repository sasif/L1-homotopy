/* ifwt.c
 *Usage: x = ifwt(w, g0, g1, J, sym);
 * Inverse fast wavelet transform.
 * Written by: Justin Romberg, Georgia Tech
 */

#include <stdlib.h>
#include "mex.h"

/* these are in fwt_level.c */
void inv_wavelet_one_level(double*, double*, int, double*, double*, int, int);
int checkPowerTwo(int, int);




void inv_wavelet_multi_level(double *x, double *w, int n, double *g0, double *g1, int l, int J, int sym) {
    int j, nj, m;
    double *tmp;
    
    tmp = (double*) mxCalloc(n, sizeof(double));
    
    /* start with the correct nj */
    nj = n;
    for (j=0; j < J; j++)
        nj >>= 1;
    
    for (j=J-1; j >= 0; j--) {
        
        nj <<= 1;
        
        if (j == J-1)
            inv_wavelet_one_level(x, w, nj, g0, g1, l, sym);
        else
            inv_wavelet_one_level(x, tmp, nj, g0, g1, l, sym);
        
        if (j > 0)
            for (m=0; m < nj; m++) {
            tmp[m] = x[m];
            tmp[nj+m] = w[nj+m];
            }
    }
    
    mxFree(tmp);
}


/* The gateway routine. */
/* x = ifwt(w, g0, g1, J, sym); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *x, *w, *g0, *g1;
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
    g0 = mxGetPr(prhs[1]);
    g1 = mxGetPr(prhs[2]);
    x = mxGetPr(plhs[0]);
    
    /* Call the C subroutine. */
    inv_wavelet_multi_level(x, w, n, g0, g1, l, J, sym);
    return;
}