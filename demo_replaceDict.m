% demo_replaceDict
%
% Solves the following BPDN problem
% min_x  tau\|x\|_1 + 1/2*||y-Ax||_2^2
% 
% and updates the solution as the A matrix is updated... 
% Infact, we can modify all the terms: y, A, and x simultaneously
%
% Applications:
%   	in dictionary learning where A is iteratively updated
%       in tracking where y, A, and x, are updated..
% 
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: August 2012

clear
close all force
 
% Limit the number of computational threads (for profiling)
% maxNumCompThreads(1);

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
T = round(M/4);    % sparsity level
Tn = 0; % round(T/10);  % update in the signal...

% rank-1 update mode
delx_mode = 'mil'; % mil or qr  
fprintf('Replace matrix A and/or signal x..\n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d.', mType, sType, SNR, N, M, T);
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
    
    %% Update the solution after updating the measurement matrix and/or the
    % sparse signal 
    
    xh_old = xh; tau_old = tau;
    x_old = x; y_old = y; A_old = A;
    
    % Update the measurement matrix
    in = []; in.type = mType;
    B = genAmat(M,N,in);
    
    A = A_old+1/5*B;
    
    % Update the sparse signal.
    in = []; in.type = sType; in.T = Tn; in.randgen = 1;
    dx = genSignal(N,in);
    x = x_old + dx;
    % e = randn(M,1)*sigma;    
    y = A*x+e;
           
    % parameter selection    
    % tau = sigma*sqrt(log(N));
    tau = max(1e-4*max(abs(A'*y)),sigma*sqrt(log(N)));
    
    homotopy_mode = 'dummy';
    switch homotopy_mode
        case 'dummy'
            % create a dummy variable...
            % use homotopy on the measurements...
            % in principle, we can start with any xh_old with this formulation      
            % and any starting value of tau or W...
            W = tau; 
            W_old = tau_old;
            % yh = A*xh_old;
            % u =  -sign(xh_old).*W_old;
            % Atdyu = A'*(yh-y)-u;
            yh = A*xh_old;
            Atr = A'*(A*xh_old-y);
            u =  -W.*sign(xh_old)-Atr;
            pk_old = Atr+u;
        otherwise            
            didp('Go back ... no escape');
    end 
    
    in = out;
    gamma_old = in.gamma;
    switch delx_mode
        case 'mil';
            % in.delx_mode = 'mil';
            % The following gram matrix and its inverse can be used from the
            % previous homotopy. Too lazy to include that right now...
            % wt BPDN homotopy update
            AtAgx = A(:,gamma_old)'*A(:,gamma_old);
            iAtAgx = inv(AtAgx);
            in.iAtA = iAtAgx;
        case {'qr','chol'};
            % in.delx_mode = 'qr';
            [Q R] = qr(A(:,gamma_old),0);
            in.Q = Q; in.R = R;
        case 'qrM'
            % in.delx_mode = 'qrM';
            [Q0 R0] = qr(A(:,gamma_old));
            in.Q0 = Q0; in.R0 = R0;
    end

    in.xh_old = xh_old;
    in.pk_old = pk_old;
    % in.Atdyu = Atdyu;
    in.u = u;
    in.W = W;
    in.W_old = W_old;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
    in.err_fun = @(z) (norm(x-z)/norm(x))^2;
    tic
    out = l1homotopy(A,y,in);
    xh_dynDict = out.x_out;
    gamma_dynDict = out.gamma;
    iter_dynDict = out.iter;
    time_dynDict = toc; 
     
    SIM_stack{sim} = [sim, tau, ...
        norm(x_old-xh)^2/norm(x)^2, sum(iter_bpdn,2), sum(time_bpdn,2), ...
        norm(x-xh_dynDict)^2/norm(x)^2, sum(iter_dynDict,2), sum(time_dynDict,2)];  
    
    fprintf('sim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; dynDict l1homotopy-%3.4g,%3.4g,%3.4g. \n', ...
        SIM_stack{sim});
   
    %% plot recovered signals
    figure(1); plot([x xh xh_dynDict]);
    title('Comparison betweeen the new and old homotopy code')
end
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: maxsim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; dynDict l1homotopy-%3.4g,%3.4g,%3.4g. \n', maxsim, mS(2:end));
