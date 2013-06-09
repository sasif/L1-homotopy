%
%----------L1 Homotopy Package-------------+
%
%  Created by: Salman Asif @ Georgia Tech.
%  Email: sasif@ece.gatech.edu
%
%  First release: 	April 2009
%  Modified: 		June 2012
%  
%  Code website: http://users.ece.gatech.edu/~sasif/homotopy
%  
%  (let's say) 
%  Version 1.1 released July 2012
%
%------------------------------------------+
%
%------------------------------------------------------------

Copyright (2012) M. Salman Asif and Justin Romberg

This file is part of L1 homotopy toolbox.

The code distributed under the terms of the GNU General Public
License 3.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose.

%------------------------------------------------------------

Modifications: July 2012

I have modified the code for BPDN_homotopy_function from the previous release. 
I have added QR and Cholesky factorization for rank-one update and created 
common functions that are used by all the homotopy problems. In the future, 
I will try to clean up all the old codes by replacing redundant scripts with 
common function calls. 

This new release contains codes for iterative and adaptive reweighting.
The codes for reweighting are the most recent. 

In this new release, you may need to compile mex codes for
1. matrix-vector product of the form A_Gamma x_Gamma and A_Gamma^T y
2. realnoiselet
3. Wavelets

See compile.m for further details. 
Add all the folders in MATLAB path or only those that are required for each
solver. 

This code is in development stage; any comments or bug reports
are very welcome.

%--------------------------------------------------------------------------
% Part of old release.
%--------------------------------------------------------------------------

First Release: April 2009

% Add L1_homotopy folder alongwith all the subfolders in your MATLAB path.
% Or only the compoenents that are needed for particular problem. 
% Either do it yourself or run 'setup_path.m' for PC

% Each folder contains functions for respective homotopy algorithms and the 
% demos which demonstrate the usage of those functions and in some case
% working of the algorithm (e.g., piecewise linear path structure)
% 
% Below we have listed main function names contained in each folder. 
% We are not writing the description of each function separately here. 
% However, the names are hopefully self-descriptive! But if you have trouble 
% in figuring things out, play with the files and suggest some changes :)

% Note that some of the simulation scripts use functions from GPSR and FPC_AS packages. 
% You must have these packages in your MATLAB path in order to run those scripts.
%
% The wavelet functions in this package use WaveLab, which must also be in
% your MATLAB path. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Pursuits_Homotopy (Standard homotopy)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basis pursuit denoising (BPDN) homotopy
BPDN_homotopy_function.m
BPDN_homotopy_demo.m

% Dantzig selector (DS) homotopy based on Primal dual pursuit
DS_homotopy_function.m
DS_homotopy_demo.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  DynamicX (Homotopy update for time-varying sparse signals)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  DynamicSeq (Homotopy update for sequential measurements)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BPDN 
DynamicSeq_BPDN_function.m
DynamicSeq_BPDN_demo.m

% DS
DynamicSeq_DS_function.m
DynamicSeq_DS_demo.m

% Simulations (use functions from GPSR and FPC_AS)
Simulation_DynamicSeq_BPDN.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Decoding (Correct sparse errors in an encoded signals)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L1 decoding
l1Decode_homotopy_fast.m
l1Decode_homotopy_qr.m
Simulation_l1Decode.m

% Robust error correction
DynamicSeq_REC_function.m
DynamicSeq_REC_demo.m
Simulation_DynamicSeq_REC.m
