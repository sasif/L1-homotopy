
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
        /* mexPrintf("Hello, world!\n");*/
    }
    return(scl);
}



void wavelet_one_level(double *w, double *x, int n, double *h0, double *h1, int l, int sym)
{
    int m, L, ln, k, kl, kx, kxp, rmkx;
    
    L = l/2;
    ln = l-L-1;
    
    for (m = 0; m < n/2; m++) {
        /* printf("m is : %d \n", m); */
        w[m] = 0;
        w[n/2+m] = 0;
        kl = 2*m-ln;
        for (k=0; k < l; k++){
            kx = kl+k;
            rmkx = 2;
            if (kx < 0) {
                if (sym == 0) {
                    kx = kx % n;
                    kxp = (kx < 0) ? n+kx : kx;
                }
                else if (sym == 1) {
                    rmkx = (-kx-1)/(n-1);
                    kxp = ((rmkx % 2) == 0) ? (((-kx-1) % (n-1))+1) : ((n-2) - ((-kx-1)%(n-1)));
                }
                else if (sym == 2) {
                    rmkx = (kx+1)/n;
                    kxp = ((rmkx % 2) == 0) ? ((-kx-1) % n) : (n - 1 - ((-kx-1) % n ));
         /* kxp = -kx - 1; */
                }
            }  else if (kx >= n) {
                if (sym == 0) {
                    kxp = kx % n;
                }
                else if (sym == 1) {
                    rmkx = (kx-1)/(n-1);
                    kxp = ((rmkx % 2) == 0) ? (((kx-1) % (n-1))+1) : ((n-2) - ((kx-1)%(n-1)));
         /* kxp = 2*n - kx - 2; */
                }
                else if (sym == 2) {
                    rmkx = kx/n;
                    kxp = ((rmkx % 2) == 0) ? (kx % n) : (n - 1 - (kx%n));
         /*kxp = 2*n - kx - 1;*/
                }
            } else {
                kxp = kx;
            }
            /*  printf("kxp is : %d l-k-1 is : %d  rmkx is: %d h0[l-k-1] is : %f   h1[l-k-1] is: %f \n", kxp, l-k-1, rmkx, h0[l-k-1], h1[l-k-1]);*/
            w[m] += x[kxp]*h0[l-k-1];
            w[n/2+m] += x[kxp]*h1[l-k-1];
        }
        
    }
}



void inv_wavelet_one_level(double *x, double *w, int n, double *g0, double *g1, int l, int sym)
{
    int m, L, ln, k, kl, kw, kwp, n2;
    double wl, wh;
    
    L = l/2-1;
    ln = l-L-1;
    n2 = n/2;
    
    for (m=0; m < n; m++) {
        
        x[m]=0;
        kl = m-ln;
        /*  printf("m is : %d \n", m); */
        for (k=0; k < l; k++) {
            
            kw = kl+k;
            if ((kw & 1) == 0) {
                kw = kw/2;
                if (sym == 0) {
                    kw = kw % n2;
                    kwp = (kw < 0) ? kw + n2 : kw;
                    wl = w[kwp];
                    wh = w[n2+kwp];
                    /*  printf("kw is : %d, kwp is : %d n2+kwp is : %d l-k-1 is : %d \n", kw, kwp, n2+kwp, l-k-1); */
                }
                else if (sym == 1) {
                    /* lowpass is 'I-II', highpass is 'II-I' */
                    kw = kw % (2*n2-1);
                    kw = (kw < 0) ? kw + 2*n2-1 : kw;
                    if (kw >= n2) {
                        wl = w[2*n2-kw-1];
                        wh = w[n2+2*n2-kw-2];
                        /*              printf("kw is : %d, wl_index is : %d, wh_index is : %d l-k-1 is : %d \n ", kw,2*n2-kw-1, n2+2*n2-kw-2, l-k-1); */
                        
                    } else {
                        wl = w[kw];
                        wh = w[n2+kw];
                        /*                      printf("kw is : %d, wl_index is : %d, wh_index is : %d l-k-1 is : %d \n ", kw,kw, n2+kw, l-k-1);*/
                        
                    }
                }
                else if (sym == 2) {
                    /* lowpass is II symmetric, highpass is II anti-symmetric */
                    kwp = kw % (2*n2);
                    kwp = (kwp < 0) ? kwp + 2*n2 : kwp;
                    if (kwp >= n2) {
                        wl = w[2*n2-kwp-1];
                        wh = -w[n2+2*n2-kwp-1];
                        /*         printf("2*n2-kwp-1 is : %d l-k-1 is : %d wl is: %d   wh is: %d \n", 2*n2-kwp-1, l-k-1, wl, wh);*/
                    } else {
                        wl = w[kwp];
                        wh = w[n2+kwp];
                        /*       printf("kwp is : %d l-k-1 is : %d  wl is: %d  wh is: %d \n", kwp, l-k-1, wl, wh);*/
                    }
                    
                    
                }
                x[m] += wl*g0[l-k-1] + wh*g1[l-k-1];
            }
            
        }
    }
}

void inv_wavelet_one_level_gen(double *x, double *w, int n, double *g0, double *g1, int l, int sym, int oddSymFlip)
{
    int m, L, ln, k, kl, kw, kwp, n2;
    double wl, wh;
    
    L = l/2-1;
    ln = l-L-1;
    n2 = n/2;
    
    for (m=0; m < n; m++) {
        
        x[m]=0;
        kl = m-ln;
        for (k=0; k < l; k++) {
            
            kw = kl+k;
            if ((kw & 1) == 0) {
                kw = kw/2;
                if (sym == 0) {
                    kw = kw % n2;
                    kwp = (kw < 0) ? kw + n2 : kw;
                    wl = w[kwp];
                    wh = w[n2+kwp];
                }
                else if (sym == 1) {
                    /* lowpass is 'I-II', highpass is 'II-I' */
                    kw = kw % (2*n2-1);
                    kw = (kw < 0) ? kw + 2*n2-1 : kw;
                    if (kw >= n2) {
                        wl = w[2*n2-kw-1];
                        wh = w[n2+2*n2-kw-2];
                        if (oddSymFlip == 1){
                            wl = w[2*n2-kw-2];
                            wh = w[n2+2*n2-kw-1];
                            /* printf("oddsymflip %d \n", oddSymFlip);*/
                        }

                    } else {
                        wl = w[kw];
                        wh = w[n2+kw];
                    }
                }
                else if (sym == 2) {
                    /* lowpass is II symmetric, highpass is II anti-symmetric */
                    kwp = kw % (2*n2);
                    kwp = (kwp < 0) ? kwp + 2*n2 : kwp;
                    if (kwp >= n2) {
                        wl = w[2*n2-kwp-1];
                        wh = -w[n2+2*n2-kwp-1];
                    } else {
                        wl = w[kwp];
                        wh = w[n2+kwp];
                    }
                }
                x[m] += wl*g0[l-k-1] + wh*g1[l-k-1];
            }
            
        }
    }
}

/* NOT USED ANYMORE */ 
void adjmod_inverse_wavelet_one_level(double *u, double *v, int nj, double *g0, double *g1, int l, int sym)
{
    int m, d, L, shift, start_midpart, start_end, ii, shift_end,scale_wave,length;
    double sum_g0, sum_g1;
    double k1, k2;
    int k1_int, k2_int, ind, check, ind1,ind2;
    /*  int *g_indices; */
    
    L = l/2-1;
    d = nj/2;
    /* g_indices = (int*) calloc(L+1, sizeof(int));
        g_indices = (int*) mxCalloc(l/2, sizeof(int)); */
    
    /* Scaling part of wavelet matrix */
    start_midpart = 0;
    scale_wave = 0;
    if (nj <= 8){
        /* For fine sacles */
        for (m = 0; m < d; m++) {
            if (m < (L+2)/2 && m < d/2){
                u[m] = 0;
                shift_end = (nj-1 < l)? (nj-1):l;
                for (shift = 0; shift<=shift_end; shift++){
                    length = 0;
                    start_end =shift%2;
                    sum_g0 = 0;
                    for (ind = start_end; ind <= 2*(L+1)-1; ind=ind+2){
                        k1 = (double)(ind-shift+2*m-L)/(double)(2*(nj-1));
                        k2 = (double)(ind-shift-2*m-L)/(double)(2*(nj-1));
                        k1_int = (floor(k1)==k1);
                        k2_int = (floor(k2)==k2);
                        check = k1_int || k2_int;
                        sum_g0 += (check*g0[ind]);
                  /*  if (k1_int || k2_int)
                    {
                        //g_indices[length] = ind;
                        //length = length+1;
                        sum_g0 += g0[ind];
                    }*/
                    }
                    /* total_taps = filter_taps_iwt(g_indices, L, nj, m, shift, start_end,scale_wave); */
                    
/*                for (ii=0; ii<length; ii++){
                    sum_g0 += g0[g_indices[ii]];
                }*/
                    u[m] += sum_g0*v[shift];
                }
            }
            else if (m >= (d-(L+2)/2)+1){
                /* for indices near nj */
                u[m] = 0;
                shift_end = (nj-1 < l)? (nj-1):l;
                for (shift = 0; shift<=shift_end; shift++){
                    start_end = (shift%2)+1;
                    length = 0;
                    sum_g0 = 0;
                    for (ind = start_end; ind <= 2*(L+1)-1; ind=ind+2){
                        k1 = (double)(ind+shift+(nj-1-2*m)-L)/(double)(2*(nj-1));
                        k2 = (double)(ind+shift-(nj-1-2*m)-L)/(double)(2*(nj-1));
                        k1_int = (floor(k1)==k1);
                        k2_int = (floor(k2)==k2);
                        check = k1_int || k2_int;
                        sum_g0 += (check*g0[ind]);
                    /*if (k1_int || k2_int)
                    {
                        //g_indices[length] = ind;
                        //length = length+1;
                        sum_g0 +=g0[ind];
                    }*/
                    }
                    /*total_taps = filter_taps_iwt(g_indices, L, nj, m, shift, start_end,scale_wave); */
                    
/*                sum_g0 = 0;
                for (ii=0; ii<length; ii++){
                    sum_g0 += g0[g_indices[ii]];
                }*/
                    u[m] += sum_g0*v[nj-shift-1];
                }
            }
            else {
                /* for indices where we have no foldings :)*/
                start_midpart = start_midpart+2;
                u[m] = 0;
                for (ii=0; ii<l; ii++){
                    u[m] += (g0[ii]*v[start_midpart+ii]);
                }
                
            }
        }
        
        /* Wavelet part of adjoint matrix */
        start_midpart = 0;
        scale_wave = 1;
        for (m = d; m < nj; m++) {
            if ((m-d) < ((L+2)/2-1) && (m-d) < d/2){
                u[m] = 0;
                shift_end = (nj-1 < l)? (nj-1):l;
                for (shift = 0; shift<=shift_end; shift++){
                    start_end = shift%2;
                    length = 0;
                    sum_g1 = 0;
                    for (ind = start_end; ind <= 2*(L+1)-1; ind=ind+2){
                        k1 = (double)(ind-shift-1+(2*(m-nj/2)+1)-L)/(double)(2*(nj-1));
                        k2 = (double)(ind-shift-1-(2*(m-nj/2)+1)-L)/(double)(2*(nj-1));
                        k1_int = (floor(k1)==k1);
                        k2_int = (floor(k2)==k2);
                        check = k1_int || k2_int;
                        sum_g1 += (g1[ind]*check);
                    /*if (k1_int || k2_int)
                    {
                        g_indices[length] = ind;
                        length = length+1;
                        sum_g1 += g1[ind];
                    }*/
                    }
                    /*   total_taps = filter_taps_iwt(g_indices, L, nj, m, shift, start_end,scale_wave); */
                    
/*                for (ii=0; ii<length; ii++){
                    sum_g1 += g1[g_indices[ii]];
                }*/
                    u[m] += sum_g1*v[shift];
                }
            }
            else if ((m-d) >= (d-(L+2)/2)){
                /* for indices near nj */
                u[m] = 0;
                shift_end = (nj-1 < l)? (nj-1):l;
                for (shift = 0; shift<=shift_end; shift++){
                    start_end = (shift%2)+1;
                    length = 0;
                    sum_g1 = 0;
                    for (ind = start_end; ind <= 2*(L+1)-1; ind=ind+2){
                        k1 = (double)(ind+shift-1+(nj-1-(2*(m-nj/2)+1))-L)/(double)(2*(nj-1));
                        k2 = (double)(ind+shift-1-(nj-1-(2*(m-nj/2)+1))-L)/(double)(2*(nj-1));
                        k1_int = (floor(k1)==k1);
                        k2_int = (floor(k2)==k2);
                        check = k1_int || k2_int;
                        sum_g1 += (g1[ind]*check);
                    /*if (k1_int || k2_int)
                    {
                        g_indices[length] = ind;
                        length = length+1;
                        sum_g1 += g1[ind];
                    }*/
                    }
                    /*total_taps = filter_taps_iwt(g_indices, L, nj, m, shift, start_end,scale_wave);*/
                    
/*                for (ii=0; ii<length; ii++){
                    sum_g1 += g1[g_indices[ii]];
                }*/
                    u[m] += sum_g1*v[nj-shift-1];
                }
            }
            else {
                /* for indices where we have no foldings :) */
                
                u[m] = 0;
                for (ii=0; ii<l; ii++){
                    u[m] += (g1[ii]*v[start_midpart+ii]);
                }
                start_midpart = start_midpart+2;
            }
        }
    }
    
    /* For coarse scales */
    else{
        for (m = 0; m < d; m++) {
            u[m] = 0;
            if (m==0){
                for(shift = 0; shift<=l/2; shift++){
                    u[m] += g0[shift+L]*v[shift];
                }
            }
            else if (m < (L+2)/2 && m >0){
                shift_end = l/2+2*m;
                for (shift = 0; shift<=shift_end; shift++){
                    ind1 = L-2*m+shift;
                    ind2 = L+2*m+shift;
                    if(ind2 > l-1)
                        u[m] += g0[ind1]*v[shift];
                    else
                        u[m] += (g0[ind1]+g0[ind2])*v[shift];
                }
            }
            else if (m >= (d-(L+2)/2)+1){
                /* for indices near nj */
                shift_end = l/2+2*(d-m-1);
                for (shift = 0; shift<=shift_end; shift++){
                    ind1 = L+(nj-1-2*m)-shift;
                    ind2 = L-(nj-1-2*m)-shift;
                    sum_g0 = (ind2 < 0) ? (g0[ind1]):(g0[ind1]+g0[ind2]);
                    u[m] += sum_g0*v[nj-shift-1];
                }
            }
            else {
                /* for indices where we have no foldings :) */
                start_midpart = start_midpart+2;
                /* u[m] = 0; */
                for (ii=0; ii<l; ii++){
                    u[m] += (g0[ii]*v[start_midpart+ii]);
                }
            }
        }
        
        /* Wavelet part of adjoint matrix */
        start_midpart = 0;
        for (m = d; m < nj; m++) {
            u[m] = 0;
            if (m==(nj-1)){
                for(shift = 0; shift<=l/2; shift++)
                    u[m] += g1[l/2-shift]*v[nj-1-shift];
            }
            else if ((m-d) < L/2 && (m-d) < d/2){
                shift_end = l/2+2*(m-d);
                for (shift = 0; shift<=shift_end; shift++){
                    ind1 = L-(2*(m-d)+1)+shift+1;
                    ind2 = L+(2*(m-d)+1)+shift+1;
                    
                    if(ind2 > l-1)
                        /*sum_g1 = g1[ind1]; */
                        u[m] += g1[ind1]*v[shift];
                    else
                        /*sum_g1 = g1[ind1]+g1[ind2];*/
                        u[m] += (g1[ind1]+g1[ind2])*v[shift];
                    /*u[m] += sum_g1*v[shift]; */
                }
            }
            else if ((m-d) >= (d-(L+2)/2) && m < nj-1){
                /* for indices near nj */
                shift_end = l/2+2*(nj-m-1);
                for (shift = 0; shift<=shift_end; shift++){
                    ind1 = L+(nj-1-2*(m-nj/2)-1)-shift+1;
                    ind2 = L-(nj-1-2*(m-nj/2)-1)-shift+1;
                    sum_g1 = (ind2 < 0) ? (g1[ind1]):(g1[ind1]+g1[ind2]);
                    u[m] += sum_g1*v[nj-shift-1];
                }
            }
            else {
                /* for indices where we have no foldings :) */
                
                u[m] = 0;
                for (ii=0; ii<l; ii++){
                    u[m] += (g1[ii]*v[start_midpart+ii]);
                }
                start_midpart = start_midpart+2;
            }
        }
    }
    /*free(g_indices); */
    /*    mxFree(g_indices); */
}
