% demo_dynamicX
%
% Solves the following BPDN problem
% min_x  tau\|x\|_1 + 1/2*||y-Ax||_2^2
% 
% and updates the solution using updated measurements 
% for a slighly modified x
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
addpath utils/qr/

disp(['--------------------',datestr(now),'-------------------------'])

% load RandomStates
% 
rseed = 2012;
rand('state',rseed);
randn('state',rseed);

% simulation parameters
mType = 'randn'; % {'randn','orth','rdct'};
sType = 'randn'; % {'randn','sign','highD', 'blocks','pcwPoly'}
SNR = inf;       % additive Gaussian noise

% reweighted setup
rwt = 5;        % number of reweighting iterations

N = 512;   % signal length
M = round(N/2);    % no. of measurements
T = round(M/5);    % sparsity level
Tn = round(T/20);  % innovations in the signal

% rank-1 update mode
delx_mode = 'mil'; % mil or qr 
fprintf('Update time-varying signal..\n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d; Tn = %d.', mType, sType, SNR, N, M, T,Tn);
disp(str0);
 
maxsim = 10;
SIM_stack = cell(maxsim,1);
SIM_memory = cell(maxsim,1);

for sim = 1:maxsim
    
    % Generate a random signal
    in = []; in.type = sType; in.T = T; in.randgen = 1;
    x = genSignal(N,in);
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
    
    maxiter = 4*N;
    err_fun = @(z) (norm(x-z)/norm(x))^2;
    
    %% BPDN using l1homotopy
    in = [];
    in.tau = tau;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
    in.err_fun = err_fun;
    tic
    out = l1homotopy(A,y,in);
    xh = out.x_out;
    iter_bpdn = out.iter;
    time_bpdn = toc;
    gamma_bpdn = out.gamma;
    err_bpdn = out.error_table;
    
    %% Update the solution for a time-varying signal (DynamicX)    
    xh_old = xh; tau_old = tau;
    x_old = x; y_old = y; 
    
    % Update the sparse signal
    in = []; in.type = sType; in.T = Tn; in.randgen = 1;
    dx = genSignal(N,in);
    x = x_old + dx;
    e = randn(M,1)*sigma;
    y = A*x+e;
    
    % parameter selection 
    % tau = sigma*sqrt(log(N));
    tau = max(1e-4*max(abs(A'*y)),sigma*sqrt(log(N)));

    
    % create a dummy variable...
    % use homotopy on the measurements...
    % in principle, we can start with any xh_old with this formulation
    % and any starting value of tau or W...
    W = tau;
    Atr = A'*(A*xh_old-y);
    u =  -W.*sign(xh_old)-Atr;
    pk_old = Atr+u;
    
    in = out;
    in.xh_old = xh_old;
    in.pk_old = pk_old;
    in.u = u;
    in.W = W;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
    in.err_fun = @(z) (norm(x-z)/norm(x))^2;
    tic
    out = l1homotopy(A,y,in);
    xh_dynX = out.x_out;
    gamma_dynX = out.gamma;
    iter_dynX = out.iter;
    time_dynX = toc;
         
    
    SIM_stack{sim} = [sim, tau, ...
        norm(x_old-xh)^2/norm(x)^2, sum(iter_bpdn,2), sum(time_bpdn,2), ...
        norm(x-xh_dynX)^2/norm(x)^2, sum(iter_dynX,2), sum(time_dynX,2)];  
    
    fprintf('sim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; dynX l1homotopy-%3.4g,%3.4g,%3.4g. \n', ...
        SIM_stack{sim});
   
    %% plot recovered signals
    figure(1); plot([x xh xh_dynX]);
    title('Comparison betweeen the new and old homotopy code')
end
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: maxsim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; dynX l1homotopy-%3.4g,%3.4g,%3.4g. \n', maxsim, mS(2:end));
