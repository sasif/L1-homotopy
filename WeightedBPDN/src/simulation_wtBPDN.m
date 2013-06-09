% Comparison of various re-weighted schemes
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% and dynamically update the weights w_i
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
% Created: April 16, 2011

clear
close all

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
addpath ../utils/qr/
addpath ../utils/utils_Wavelet
addpath(genpath('C:\Users\asif\MATLAB\Sparsity_solvers'));

% load RandomStates
rseed = sum(100*clock);
rseed = 0;
% rand('state',rseed);
% randn('state',rseed);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));

%% Sim parameters
disp(['--------------------',datestr(now),'-------------------------'])

% reweighted setup
rwt = 5; % number of iterative reweighted iterations
rwt_adp = 2; % number of adaptive reweighted iterations

% rank-1 update mode 
delx_mode = 'mil'; % mil or qr

% simulation setup
maxsim = 200;

N_list = [256 512 1024 2048];
M_ratio = [2:5];
T_ratio = [3:5];
sType_list = {'randn','sign'}; % 'blocks','pcwPoly', 'HeaviSine'
mType_list = {'randn','orth'};
SNR_list = 30; % [20:10:40 1e16];
rwt_list = 2; % [1 2 3 4];
lambda_list = 0; % [0, 1e-1, 1e-2, 1e-4]; % multiplication factor for tau = lambda *\|A'y\|_\infty

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

for mT = 1:length(mType_list);
    for sT = 1:length(sType_list);
        for snr = 1:length(SNR_list);
            for rwt_mode = rwt_list
                for lam = 1:length(lambda_list)
                    for nl = 1:length(N_list)
                        for mr = 1:length(M_ratio)
                            for tr = 1:length(T_ratio)
                                
                                N = N_list(nl);     % signal length
                                M = round(N/M_ratio(mr));   % no. of measurements
                                T = round(M/T_ratio(tr));   % sparsity level
                                
                                sType = char(sType_list{sT});
                                mType = char(mType_list{mT});
                                SNR = SNR_list(snr);
                                
                                lambda = lambda_list(lam);
                                
                                str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d, rwt_mode-%d, lambda%3.4g.', mType, sType, SNR, N, M, T, rwt_mode, lambda);
                                disp(str0);
                                
                                %% Run simulation
                                script_simulation_wtBPDN_ALL;
                                
                                %                         %% Plot results...
                                %                         EXP = {'adp','rwt','yall1','sparsa','spgl1'};
                                %                         exp_plot = [1 2 3 4 5];
                                %                         figure(101);
                                %                         for exp_no = 1:length(exp_plot)
                                %                             exp = exp_plot(exp_no);
                                %                             sim_err = [];
                                %                             sim_itr = [];
                                %                             for ii = 1:maxsim
                                %                                 sim_err = [sim_err; [SIM_memory{ii}{exp,3}]];
                                %                                 sim_itr = [sim_itr; [SIM_memory{ii}{exp,2}]];
                                %                             end
                                %                             eval(sprintf('err_%s = sim_err;',char(EXP{exp})));
                                %                             eval(sprintf('iter_%s = sim_itr;',char(EXP{exp})));
                                %                             subplot(2,length(exp_plot),exp_no); plot(sim_err); title(sprintf('rel. err for %s',char(EXP{exp})));
                                %                             subplot(2,length(exp_plot),length(exp_plot)+exp_no); plot(sim_itr); title(sprintf('iter count for %s',char(EXP{exp})));
                                %                         end
                                
                                
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
                                
                                % eval(sprintf('save Results_%s',mfilename));
                            end
                        end
                    end
                end
            end
        end
    end
end