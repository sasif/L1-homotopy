
#include <stdlib.h>
#include "mex.h"
#include <math.h>

/* check that the vector length is a power of 2,
   just using bitshifting instead of log */
int checkPowerTwo(int m, int J)
{
    int scl=0;
    
  /* check that it's not a degenerate 0 by 1 vector or singleton */
    if (m <= 1) {
        mexErrMsgTxt("Vector length must be greater than 1.");
    }
  /* keep dividing by two until result is odd */
    while( (m & 1) == 0){
        m >>= 1;
        scl++;
    }
  /* check that m is not a multiple of an odd number greater than 1 */
  /*  if (m > 1) {
        mexErrMsgTxt("Vector length must be power of 2.");
    } */ 
  /* check that scl is less than Jmax */
    if (scl < J){
        mexErrMsgTxt("Image is not dyadic upto the specified scale---use smaller J.");
    }
    return(scl);
}

/*
void adj_wavelet_one_level_old(double *u, double *w, int nj, double *h0, double *h1, int l, int sym) {
    int m, d, L, shift, ii, total_taps, start, shift_end, m_even_odd;
    double sum_h0, sum_h1;
    int *h_indices;
    
    L = l/2;
    d = nj/2;
    h_indices = (int*) mxCalloc(L, sizeof(int));
    
    for (m = 0; m < nj; m++) {
        if (m <= L && m < d){
            u[m] = 0;
            shift_end = (nj-2 < l)? (nj-2):l;
            for (shift = 0; shift<=shift_end; shift=shift+2){
                m_even_odd = ((m+L)%2);
                total_taps = filter_taps(h_indices, L, nj, m, shift, sym, m_even_odd);
                
                sum_h0 = 0;
                for (ii=0; ii<total_taps; ii++){
                    sum_h0 += h0[h_indices[ii]];
                }
                sum_h1 = 0;
                for (ii=0; ii<total_taps; ii++){
                    sum_h1 += h1[h_indices[ii]];
                }
                u[m] += (sum_h0*w[shift/2] + sum_h1*w[shift/2+d]);
            }
        }
        else if (m >= nj-(L+1)){
            u[m] = 0;
            shift_end = (nj-1 < l)? (nj-1):l;
            for (shift = nj-2; shift>=(nj-2-shift_end); shift=shift-2){
                
                m_even_odd = ((m+L)%2);
                total_taps = filter_taps(h_indices, L, nj, m, shift, sym, m_even_odd);
                sum_h0 = 0;
                for (ii=0; ii<total_taps; ii++){
                    sum_h0 += h0[h_indices[ii]];
                }
                sum_h1 = 0;
                for (ii=0; ii<total_taps; ii++){
                    sum_h1 += h1[h_indices[ii]];
                }
                u[m] += (sum_h0*w[shift/2] + sum_h1*w[shift/2+d]);
            }
            
        }
        else {
            // for indices where we have no foldings :)
            //start = (int)(ceil(m/2)-(L+1)/2+1);
            
            start = (ceil((m+1-L)/2)); // works fine
            
            u[m] = 0;
            
            if ( (m+L-1) % 2 == 0){
                for (ii=0; ii<L; ii++){
                    u[m] += (h0[2*ii+1]*w[start+ii] + h1[2*ii+1]*w[d+start+ii]);
                }
            }
            else{
                for (ii=0; ii<L; ii++){
                    u[m] = u[m] + h0[2*ii]*w[start+ii] + h1[2*ii]*w[d+start+ii];
                }
            }
        }
        
    }
    mxFree(h_indices);
}
*/

/* The key idea is to look at symmetrization and convolution as two separate steps */ 
void adj_wavelet_one_level(double *u, double *w, int nj, double *h0, double *h1, int l, int sym) {
    int m, d, k, L, Ly, m_even_odd, kl, a, cm, cn;
    double *y;
    
    L = l/2;
    d = nj/2;
    Ly = nj+l-1;
    y = (double*) mxCalloc(Ly, sizeof(double));
    
    /* Decompose symmetric extension and filtering into two steps */
    
    /* adjoint for linear convolution */
    for (m = 0; m < Ly; m++) {
        y[m] = 0;
        m_even_odd = ((m)%2)+1;
        if ( m < (l-1) && m < nj){
            for (k = 0; k <= m; k = k+2){
                y[m] += h0[l-k-m_even_odd]*w[m/2-k/2] + h1[l-k-m_even_odd]*w[d+m/2-k/2];
                /* mexPrintf("\n\t First: filter index = %i, wave index = %i, m_even_odd = %i ",l-k-m_even_odd, m/2-k/2,m_even_odd);*/
            }
        }
        else if (m >= nj){
            kl = ((Ly-m-2) < (nj-2))? (Ly-m-2) : (nj-2);
            for (k = 0; k <= kl; k = k+2){
                y[m] += h0[Ly-m-k-2]*w[nj/2-k/2-1] + h1[Ly-m-k-2]*w[d+nj/2-k/2-1];
                /* mexPrintf("\n\t Last: filter index = %i, wave index = %i, m_even_odd = %i ",Ly-m-k-2, nj/2-k/2-1,m_even_odd); */
            }
        }
        else{
            for (k = 0; k <= l-1; k = k+2){
                y[m] += h0[l-k-m_even_odd]*w[m/2-k/2] + h1[l-k-m_even_odd]*w[d+m/2-k/2];
                /* mexPrintf("\n\t Middle: filter index = %i, wave index = %i, m_even_odd = %i ",l-k-m_even_odd, m/2-k/2,m_even_odd); */
            }
        }
    }
    
    /* clear memory for output */
    for (cn=0; cn<nj; cn++) {
        u[cn] = 0;
    }
    
    /* adjoint for symmetric extension */
    if (sym == 0){
        for (m=0; m<Ly; m++) {
            cm = m+1-(l/2);
            if (cm<0) {
                a = cm%nj;
                cn = (a < 0) ? a+nj: a;
            }
            else if (cm > nj-1){
                cn = cm % (nj);
            }
            else{
                cn = cm;
            }
            u[cn] = u[cn]+y[m];
            /* mexPrintf("\n\t nj = %i, periodic indices are = %i, convolution index = %i",nj, cn, cm); */
        }
    }
    if (sym == 1){
        for (m = 0; m < Ly; m++) {
            cm = m+1-(l/2);
            if (cm < 0){
                a = (-cm) % (2*(nj-1));
                cn = (a > (nj-1)) ? (2*(nj-1)-a) : (a);
            }
            else if (cm > (nj-1)){
                a = (cm) % (2*(nj-1));
                cn = (a > (nj-1)) ? (2*(nj-1)-a) : (a);
            }
            else {
                cn = cm;
            }
            /* mexPrintf("\n\t cn = %i, Value u[cn] = %lf, y[m] = %lf ",cn, u[cn],y[m]); */
            u[cn] = u[cn]+y[m];
        }
    }
    if (sym == 2){
        for (m = 0; m < Ly; m++) {
            cm = m+1-(l/2);
            if (cm < 0){
                a = (-cm-1) % (2*nj);
                cn = (a > (nj-1)) ? (2*nj-a-1) : (a);
            }
            else if (cm > (nj-1)){
                a = (cm) % (2*nj);
                cn = (a > (nj-1)) ? (2*nj-1-a) : (a);
            }
            else {
                cn = cm;
            }
            u[cn] = u[cn]+y[m];
        }
    }
    mxFree(y);
}