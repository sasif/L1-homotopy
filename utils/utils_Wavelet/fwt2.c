/* Fast 2D wavelet transform
 * Wriiten by: Justin Romberg, Georgia Tech
 * Modified by Salman Asif, Georgia Tech
 * 
 * Revisions: 
 * December 2011:   Added the ability to work with non-square images that 
 *                  are dyadic upto scale J in both dimensions
 */

#include <stdlib.h>
#include "mex.h"

/* these are in fwt_level.c */
void wavelet_one_level(double*, double*, int, double*, double*, int, int);
int checkPowerTwo(int, int);


void wavelet_multi_level2D(double *w, double *x, int nR, int nC, double *h0, double *h1, int l, int J, int sym) {
    int j, nRj, nCj, n, m1, m2;
    double *tmp, *tmpo;
    
    n = (nR > nC) ? nR: nC;
    /* Arrays for row/col processing */
    tmp = (double*) mxCalloc(n, sizeof(double));
    tmpo = (double*) mxCalloc(n, sizeof(double));
    
    nRj = nR;
    nCj = nC;
    for (j=0; j < J; j++) {
        /* columns */
        for (m2=0; m2 < nCj; m2++) {
            if (j==0)
                wavelet_one_level(&w[m2*nR], &x[m2*nR], nR, h0, h1, l, sym);
            else {
                for (m1=0; m1 < nRj; m1++)
                    tmp[m1] = w[nR*m2+m1];
                wavelet_one_level(&w[m2*nR], tmp, nRj, h0, h1, l, sym);
            }
        }
        /* rows */
        for (m1=0; m1 < nRj; m1++) {
            for (m2=0; m2 < nCj; m2++)
                tmp[m2] = w[nR*m2+m1];
            wavelet_one_level(tmpo, tmp, nCj, h0, h1, l, sym);
            for (m2=0; m2 < nCj; m2++)
                w[nR*m2+m1] = tmpo[m2];
        }
        
        nRj >>= 1;
        nCj >>= 1;
    }
    
    mxFree(tmp);
    mxFree(tmpo);
}

/* The gateway routine. */
/* w = fwt2(x, h0, h1, J, sym); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *x, *w, *h0, *h1;
    int nR, nC, l, l1, sym, J, JmaxR, JmaxC;
    
    /* Check for the proper number of arguments. */
    if (nrhs != 5) {
        mexErrMsgTxt("Exactly five inputs required");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments");
    }
    if (mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input x must be real");
    }
    
    nR = mxGetM(prhs[0]);
    nC = mxGetN(prhs[0]);
    /*
     * if (n != n1) {
     * mexErrMsgTxt("Input x must be a square image.");
     * }
     */
    
    J = mxGetScalar(prhs[3]);
    JmaxR = checkPowerTwo(nR, J);
    JmaxC = checkPowerTwo(nC, J);
    if ((J < 0) || (J > JmaxR) || (J > JmaxC)) {
        mxErrMsgTxt("Input J must be an integer between 0 and log2(n), and dyadic for ROW and COL---use smaller J.");
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
    plhs[0] = mxCreateDoubleMatrix(nR, nC, mxREAL);
    
    /* Assign pointers to each input and output. */
    x = mxGetPr(prhs[0]);
    h0 = mxGetPr(prhs[1]);
    h1 = mxGetPr(prhs[2]);
    w = mxGetPr(plhs[0]);
    
    /* Call the C subroutine. */
    wavelet_multi_level2D(w, x, nR, nC, h0, h1, l, J, sym);
    return;
}