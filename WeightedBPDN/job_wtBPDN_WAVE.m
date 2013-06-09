% Comparison of various solvers for iterative reweighting and adaptive reweighting
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% while updating the weights w_i
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
% script_simulation_wtBPDN
%
% The paper has results for Blocks and HeaviSine at SNR = 40
% Run publish_results.m in the misc folder to generate pdfs presented int the paper.
%
% NOTE: While preparing the code for release, I realized that I used the average
% of SER in db when I compiled results for publication, which is probably 
% not the correct way to compute the average. 
%
% The average SER in the paper was computed as mean(-10log10(normalized error))
% The standard way to compute average is -10log10(mean(normalized squared error))
% where mean is over all the experiments. 
%
% To reproduce results in the paper, use publish_results.m
% Otherwise use publish_results_new.m
%
% This changes the values of average SER for all the solvers. However, the 
% comparison among the solvers remains almost same. 


% function job_wtBPDN_WAVE(mT, sT, snr, rwt_mode, lam_mode)
% mT = 1, sT = 2, snr = 3, rwt_mode = 5, lam_mode = 1;

clear

EXP_LIST = [1 1 2 5 1; 1 1 3 5 1; 1 1 4 5 1; 1 2 2 5 1; 1 2 3 5 1; 1 2 4 5 1];
maxNumCompThreads(1);
% numCores = str2double(getenv('NUMBER_OF_PROCESSORS'));
% mpSize = 1;% size(EXP_LIST,1);
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
addpath ../utils/utils_Wavelet
addpath ../solvers/
addpath src/

%% Simulation parameters
% reweighted setup
rwt = 5;        % number of reweighting iterations
rwt_adp = 0;    % number of reweighting iterations after adaptive reweighting

% rank-1 update mode
delx_mode = 'mil'; % mil or qr

% simulation setup
maxsim = 100;

mType_list = {'randn','orth'};
sType_list = {'blocks','HeaviSine','pcwPoly'};
SNR_list = [20:10:40 inf];

lambda_list = [0, 1e-1, 1e-2, 1e-4]; 
% multiplication factor for tau = lambda *\|A'y\|_\infty
% 0 --> sigma*log(N)

N_list = [256 512 1024];
M_ratio = [2:.5:4];

for pf = 1:length(EXP_LIST)
    
    mT = EXP_LIST(pf,1);
    sT = EXP_LIST(pf,2);
    snr = EXP_LIST(pf,3);
    rwt_mode = EXP_LIST(pf,4);
    lam_mode = EXP_LIST(pf,5);
    
    EXP_stack = {};
    estack = 1;
    EXP_stack{1,1} = 'mType';
    EXP_stack{1,2} = 'sType';
    EXP_stack{1,3} = 'SNR';
    EXP_stack{1,4} = 'rwt_mode';
    EXP_stack{1,5} = 'lambda_mode';
    EXP_stack{1,6} = 'N';
    EXP_stack{1,7} = 'M';
    EXP_stack{1,8} = 'T';
    EXP_stack{1,9} = 'str0';
    EXP_stack{1,10} = 'str2';
    EXP_stack{1,11} = 'avg SIM_stack';
    EXP_stack{1,12} = 'SIM_stack';
    EXP_stack{1,13} = 'SIM_memory';
    
    % for mT = 1:length(mType_list);
    %     for sT = 1:length(sType_list);
    %         for snr = 1:length(SNR_list);    
    %                 for lam_mode = 1:length(lambda_list)
    
    for nl = 1:length(N_list)
        for mr = 1:length(M_ratio)
            
            N = N_list(nl);     % signal length
            M = round(N/M_ratio(mr));   % no. of measurements
            T = N;   % sparsity level
            
            sType = char(sType_list{sT});
            mType = char(mType_list{mT});
            SNR = SNR_list(snr);
            
            lambda = lambda_list(lam_mode);
            
            str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d, rwt_mode-%d, lambda%3.4g.', mType, sType, SNR, N, M, T, rwt_mode, lambda);
            disp(str0);
            
            % filename_save = sprintf('results_comparison_ALL/comparison_wtBPDN_mType%s_sType%s_SNR%d_rwt_mode%d_lambda%3.4g_N%d_M%d_T%d.mat',mType,sType,SNR,rwt_mode,lambda,N,M,T);
            % filename_save = sprintf('results_comparison_ALL/comparison_wtBPDN_mT-%s_sT-%s_SNR%d_rwt_mode%d_lam_mode%d.mat',mType,sType,SNR,rwt_mode,lam_mode);
            
            filename_save = sprintf('results_comparison_ALL/comparison_wtBPDN_mT-%s_sT-%s_SNR%d-reproduce-Trwt.mat',mType,sType,SNR);
            % filename_save = sprintf('results_comparison_ALL/comparison_wtBPDN_mT-%s_sT-%s_SNR%d-reproduce-rwt.mat',mType,sType,SNR);
            
            % load RandomStates
            %
            % Experiments in the earlier draft...
            % rseed = 0;
            % rseed = sum(100*clock);
            % RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));            
            %
            % Experiments that reproduce results in the paper...
            rseed = 2012;
            rand('state',rseed);
            randn('state',rseed);
            
            % rng(rseed);
            % rseed = rng(2012);
            
            %% Run simulation
            script_simulation_wtBPDN
            
            %% save results
            estack = estack+1;
            
            EXP_stack{estack,1} = mType;
            EXP_stack{estack,2} = sType;
            EXP_stack{estack,3} = SNR;
            EXP_stack{estack,4} = rwt_mode;
            EXP_stack{estack,5} = lambda;
            EXP_stack{estack,6} = N;
            EXP_stack{estack,7} = M;
            EXP_stack{estack,8} = T;
            EXP_stack{estack,9} = str0;
            EXP_stack{estack,10} = str2;
            EXP_stack{estack,11} = mean(cell2mat(SIM_stack),1);
            EXP_stack{estack,12} = SIM_stack;
            EXP_stack{estack,13} = SIM_memory;
            
            eval(sprintf('save %s', filename_save));
        end
    end
end