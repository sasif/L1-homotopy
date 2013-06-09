%
%----------L1 Homotopy Package-------------+
%
%  Created by: Salman Asif @ Georgia Tech.
%  Email: sasif@ece.gatech.edu
%
%------------------------------------------+
%

% Add L1_homotopy folder alongwith all the subfolders in your MATLAB path.
% Either do it yourself or run 'set_path.m' for PC

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
