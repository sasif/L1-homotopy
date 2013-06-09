% Comparison of various solvers for adaptive reweighting
% with largescale examples
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% while adaptively selecting the weights w_i
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: June 16, 2011
%
% Reference: 
% "Fast and accurate algorithms for re-weighted L1 norm minimization," by 
% M. Salman Asif and Justin Romberg
% 
% To reproduce experiments in the paper, use this scrip with 
% 
% rseed = 2012;
% rand('state',rseed);
% randn('state',rseed);
%  
% before running each simulation in 
% script_simulation_adpWBPDN     

clear
% close all force

%% parfor setup
% numCores = str2double(getenv('NUMBER_OF_PROCESSORS'));
% mpSize = numCores-1;
% if matlabpool('size') ~= mpSize
%     if matlabpool('size')~=0
%         matlabpool close;
%     else
%         matlabpool('open', mpSize);
%     end
% end

%% Setup path
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath ../Pursuits_Homotopy/
addpath ../utils/
addpath ../utils/utils_Wavelet/
addpath ../utils/utils_meas/
addpath ../solvers/
addpath src/

disp(['--------------------',datestr(now),'-------------------------'])

% simulation parameters
rwt_mode = 5;
lambda = 0;

SAVE_RESULTS = false; 

% SAVE_RESULTS = true;  diary(sprintf('%s-reproduce.txt',mname));

largescale = 1;

for SNR = [40]
    for M = [25000 30000] % [2 3]        
        IMG_LIST = {'barbara','boats', 'cameraman','house','peppers','shapes','lena','airplane','baboon','sailboat','tiffany'};
        for img = 1:5; % length(IMG_LIST)
            sType = IMG_LIST{img};
            
            if largescale
                mType = 'noiselets';
                N = (256)^2;
                % M = round(N/R);
                T = N;
                str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M) = %d, %d, rwt_mode-%d, lambda%3.4g.', mType, sType, SNR, N, M, rwt_mode, lambda);
                disp(str0);
            else
                mType = 'randn';
                sType = 'HeaviSine';
                N = 512;   % signal length
                M = round(N/2);    % no. of measurements
                T = round(M/3);    % sparsity level
                str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d, rwt_mode-%d, lambda%3.4g.', mType, sType, SNR, N, M, T, rwt_mode, lambda);
                disp(str0);
            end
            
            % rank-1 update mode
            delx_mode = 'mil'; % mil or qr
            
            %% Simulation
            maxsim = 10;
            script_simulation_adpWBPDN            
            
        end
    end
end
diary off;