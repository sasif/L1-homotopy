% demo_posBPDN
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \tau \|x\|_1 + 1/2*||y-Ax||_2^2
% with a positivity constraint... i.e., x >= 0
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: February 2013

clear
close all force

%% Setup path
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath utils/ 

disp(['--------------------',datestr(now),'-------------------------'])

% load RandomStates
%
% rseed = 0;
% rseed = sum(100*clock);
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));

rseed = 2012;
rand('state',rseed);
randn('state',rseed);

% simulation parameters
mType = 'randn'; % {'randn','orth','rdct'};
sType = 'randn'; % {'randn','sign','highD', 'blocks','pcwPoly'}
SNR = 40;       % additive Gaussian noise

N = 256;   % signal length
M = round(N/4);    % no. of measurements
T = round(M/4);    % sparsity level

% rank-1 update mode 
delx_mode = 'mil'; % mil or qr 
fprintf('Standard LASSO/BPDN homotopy..\n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d.', mType, sType, SNR, N, M, T);
disp(str0);

maxsim = 10;
SIM_stack = cell(maxsim,1);
                              
 
for sim = 1:maxsim
    
    % Generate a random signal
    in = []; in.type = sType; in.T = T; in.randgen = 1;
    x = genSignal(N,in); 
    x = abs(x);
    [val ind] = sort(abs(x),'descend');
    ind_pos = ind(val>1e-1);
    gamma_orig = ind_pos(1:min(length(ind_pos),M-1));
    

    % measurement matrix
    in = []; in.type = mType;
    A = genAmat(M,N,in);
    
    % measurements
    sigma = sqrt(norm(A*x)^2/10^(SNR/10)/M);    
    e = randn(M,1)*sigma;
    y = A*x+e;
    
    %     % orthogonalize rows of A
    %     [Q, R] = qr(A',0);
    %     A = Q'; y = R' \ y;
    
    % parameter selection 
    % tau = sigma*sqrt(log(N));
    tau = max(1e-4*max(abs(A'*y)),sigma*sqrt(log(N)));
    
    err_fun = @(z) (norm(x-z)/norm(x))^2;
    maxiter = 4*N;
    
    %% BPDN using l1homotopy
    in = [];   
    in.tau = tau;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
    in.err_fun = err_fun;
    out = l1homotopy(A,y,in);
    xh = out.x_out;
    iter_bpdn = out.iter;
    time_bpdn = out.time;
    gamma_bpdn = out.gamma;
    err_bpdn = out.error_table;
        
    %% BPDN homotopy with positivity constraint
    in = [];
    in.tau = tau;
    in.nonneg = 1; % Positivity constraint flag
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
    in.err_fun = err_fun;
    out = l1homotopy(A,y,in);
    xh_nonneg = out.x_out;
    iter_nonneg = out.iter;
    time_nonneg = out.time;
    gamma_nonneg = out.gamma;
    err_nonneg = out.error_table;
    
    SIM_stack{sim} = [sim, tau, ...
        norm(x-xh)^2/norm(x)^2, sum(iter_bpdn,2), sum(time_bpdn,2), ...
        norm(x-xh_nonneg)^2/norm(x)^2, sum(iter_nonneg,2), sum(time_nonneg,2)];  
    
    fprintf('sim %d. tau = %3.4g, (err,iter,time): l1homotopy-%3.4g,%3.4g,%3.4g; nonneg homotopy-%3.4g,%3.4g,%3.4g. \n', ...
        SIM_stack{sim});

    %% plot recovered signals
    figure(1); plot([x xh xh_nonneg]);    
    title('Comparison betweeen the new and nonneg homotopy code')
end
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: maxsim %d. tau = %3.4g, (err,iter,time): l1homotopy-%3.4g,%3.4g,%3.4g; nonneg homotopy-%3.4g,%3.4g,%3.4g. \n', maxsim, mS(2:end));
