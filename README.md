L1-homotopy
===========

Codes related to L1-norm minimization using homotopy

%----------L1 Homotopy Package-------------+

Created by: Salman Asif @ Georgia Tech.  
Email: sasif@gatech.edu  

  
Code website: http://users.ece.gatech.edu/~sasif/homotopy  
and  
Github repository: https://github.com/sasif/L1-homotopy  
  
  (let's say this is)  
  	Version 2.0 released June 2013  
  	
	Previous versions:   
  	v1.1 release:	July 2012
  	v1.0 release: 	April 2009
	
%------------------------------------------+

References: 

	- M. Salman Asif and Justin Romberg, "Sparse recovery of streaming signals using L1-homotopy,"   
	  			preprint available at http://users.ece.gatech.edu/~sasif/ and ???  

	- M. Salman Asif, "Dynamic compressive sensing: Sparse recovery algorithms for streaming signals and video." 
	 			Doctoral Thesis, Georgia Institute of Technology, 2013. 

%------------------------------------------------------------

L1-homotopy is a highly versatile homotopy program that can solve a variety of L1-norm minimization problems using a warm start.   

	We consider the following linear model of observations: 
		
		y = Ax+e
		
	where x is a sparse vector of interest, A is the system matrix, y is the measurement vector, and e is the noise.   
	We want to solve the following weighted L-norm minimization program:  
	
	minimize_x  ||W x||_1 + 1/2*||Ax-y||_2^2,  
	
	where W is a diagonal matrix with positive weights.
	
	Instead of solving the optimization program from scratch, we use a vector xh_old as the starting point   
	and solve the following homotopy program:  
	
	minimize_x  ||W x||_1 + 1/2*||Ax-y||_2^2 + (1-epsilon)u'x,  

	u is defined as u = -W*sign(xh_old)-A'*(A*xh_old-y)   
	xh_old is an arbitrary warm-start vector (or a zero vector if no warm-start is available)  

	
	l1homotopy.m is the main function that solves the following homotopy program by changing epsilon from 0 to 1, 
	the solution of homotopy program changes from the warm-start vector (xh_old) to the desired solution. 
	Further details in Section IV of the paper or Chapter III of the thesis.  


Scripts for different problems are also included in this package to demonstrate the use of l1homotopy:   
 
demo_BPDN -- solves LASSO/BPDN problem from scratch  
demo_posBPDN -- solves BPDN problem with positive sign constraint on the estimate  
demo_dynamicX -- updates the solution for BPDN as the signal changes  
demo_dynamicSeq -- updates the signal as sequential measurements are added  
demo_rwtL1-- solves iterative reweighting for BPDN  
demo_dynamicRWT -- iteratively updates weights in the L1-norm cost while estimating a time-varying signal  
demo_streamingLOT -- iteratively estimates a streaming signal using lapped orthogonal transform as the representation basis  
demo_KalmanRWT -- iteratively estimates a streaming signal that follows a linear dynamic model  
and more...   


You may need to compile mex codes for  
1. matrix-vector product of the form A_Gamma x_Gamma and A_Gamma^T y  
2. realnoiselet  
3. Wavelets
  
See compile.m for further details.   
Add all the folders in MATLAB path or only those that are required for each solver.   

This code is in development stage; any comments or bug reports are very welcome.  

%------------------------------------------------------------  

This file is part of L1 homotopy toolbox.  
Copyright (C) 2013, M. Salman Asif, all rights reserved.  

Redistribution and use of this code, with or without modification, are permitted provided that the following conditions are met:  

The software is provided under the terms of this license strictly for academic, non-commercial, not-for-profit purposes. Redistributions of source code must retain the above copyright notice, this list of conditions (license) and the following disclaimer. The name of the author may not be used to endorse or promote products derived from this software without specific prior written permission.   

This software is being provided "as is", without any express or implied warranty. In particular, the authors do not make any representation or warranty of any kind concerning the merchantability of this software or its fitness for any particular purpose.  

%------------------------------------------------------------  



%------------------------------------------------------------

Parts of old releases

%------------------------------------------------------------


Other than L1 decoding and adaptive reweighting methods, these homotopy programs can also be solved 
using l1homotopy.m (for the LASSO/BPDN formulation).


%---------------  
% Pursuits_Homotopy (Standard homotopy solvers)  
%---------------  

% Basis pursuit denoising (BPDN) homotopy  
BPDN_homotopy_function.m  
BPDN_homotopy_demo.m  

% Dantzig selector (DS) homotopy based on primal-dual pursuit  
DS_homotopy_function.m  
DS_homotopy_demo.m  
 
%---------------  
% DynamicX (Homotopy update for time-varying sparse signals)  
%---------------  

% BPDN  
DynamicX_BPDN_function.m  
DynamicX_BPDN_demo.m  
DynamicX_BPDN_Visual_demo.m  

% DS  
DynamicX_DS_function.m  
DynamicX_DS_demo.m  

% Simulations (use functions from GPSR and FPC_AS)  
Simulation_DynamicX_BPDN.m  
Simulation_DynamicX_BPDN_Pathological.m  
Simulation_DynamicX_BPDN_Wavelet.m % This script uses WaveLab functions.  

%---------------  
% DynamicSeq (Homotopy update for sequential measurements)  
%---------------  

% BPDN   
DynamicSeq_BPDN_function.m  
DynamicSeq_BPDN_demo.m  

% DS  
DynamicSeq_DS_function.m  
DynamicSeq_DS_demo.m  
 
% Simulations (use functions from GPSR and FPC_AS)  
Simulation_DynamicSeq_BPDN.m  
 
%---------------  
% Decoding (Correct sparse errors in an encoded signals)  
%---------------  

% L1 decoding  
l1Decode_homotopy_fast.m  
l1Decode_homotopy_qr.m  
Simulation_l1Decode.m   

% Robust error correction  
DynamicSeq_REC_function.m  
DynamicSeq_REC_demo.m  
Simulation_DynamicSeq_REC.m  
  

%-----------------  
% Iterative and adaptive reweighting (WeightedBPDN folder)  
%-----------------  
The codes are based on iterative and adaptive reweighting methods described in the following paper:  
"Fast and accurate algorithms for re-weighted L1 norm minimization," by M. Salman Asif and Justin Romberg  
  
More details in the readme.txt file inside the folder  
