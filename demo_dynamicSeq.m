% demo_dynamicSeq
%
% Solves the following BPDN problem
% min_x  tau\|x\|_1 + 1/2*||y-Ax||_2^2
% 
% and updates the solution as new measurements w = Bx are received. 
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: August 2012
%
% Modified: May 2013 (using l1homotopy for update)
% 
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
  
N = 512;   % signal length
M = round(N/2);    % no. of measurements
T = round(M/4);    % sparsity level
Mn = round(M/10);  % new measurements

% rank-1 update mode
delx_mode = 'mil'; % mil or qr 
fprintf('Add sequential measurements..\n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d; Mn = %d.', mType, sType, SNR, N, M, T,Mn);
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
    A = A*sqrt(M); % Because genAmat generates matrices with unit (expected) norm columns
    
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
    
    %% Update the solution after adding measurements (DynamicSeq)
    xh_old = xh; tau_old = tau;
    x_old = x; y_old = y; A_old = A;
    
    % Additional measurements
    in = []; in.type = mType;
    B = genAmat(Mn,N,in);
    B = B*sqrt(Mn); % Because genAmat generates matrices with unit (expected) norm columns
    
    % noise in the new measurements
    sigma = sqrt(norm(B*x)^2/10^(SNR/10)/Mn);    
    e = randn(Mn,1)*sigma;
    w = B*x+e;
    
    A = [A_old; B];
    y = [y_old; w];
    
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
            Atr = A'*(A*xh_old-y);
            u =  -W.*sign(xh_old)-Atr;
            pk_old = Atr+u;
        otherwise
            % Use direct information from the old solution
            W = tau; 
            W_old = tau_old;
            %             u = out.pk;
            %             yh = [y_old; B*xh_old];
            %             Atdyu = A'*(yh-y);
            Atr = A'*(A*xh_old-y);
            u =  -W.*sign(xh_old)-Atr;
            pk_old = Atr+u;
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
    xh_dynSeq = out.x_out;
    gamma_dynSeq = out.gamma;
    iter_dynSeq = out.iter;
    time_dynSeq = toc;
         
    SIM_stack{sim} = [sim, tau, ...
        norm(x_old-xh)^2/norm(x)^2, sum(iter_bpdn,2), sum(time_bpdn,2), ...
        norm(x-xh_dynSeq)^2/norm(x)^2, sum(iter_dynSeq,2), sum(time_dynSeq,2)];  
    
    fprintf('sim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; dynSeq l1homotopy-%3.4g,%3.4g,%3.4g. \n', ...
        SIM_stack{sim});
   
    %% plot recovered signals
    figure(1); plot([x xh xh_dynSeq]);
    title('Comparison betweeen the new and old homotopy code')
end
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: maxsim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; dynSeq l1homotopy-%3.4g,%3.4g,%3.4g. \n', maxsim, mS(2:end));
