% demo_dynamicRWT
%
% Solves the following dynamic BPDN problem
% min_x  \|Wx\|_1 + 1/2*||y-Ax||_2^2  
%
% which updates W and the solution as the signal changes 
%
%   
% Applications:
%   	tracking where y, A, and x, are updated..
%       predict an estimate of the solution and
%       update weights according to the predicted solution
%
% We can also add a difference operator ||x-xh||_2^2 in the reconstruction
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: August 2012

clear
close all force
 
%% Setup path
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath utils/ 
addpath utils/utils_wavelet/
addpath solvers/
% addpath WeightedBPDN/src/

disp(['--------------------',datestr(now),'-------------------------'])

% load RandomStates 
%
rseed = 2012;
rand('state',rseed);
randn('state',rseed);

% simulation parameters
mType = 'randn'; % {'randn','orth','rdct'};
sType = 'pcwpoly'; % {'randn','sign','highD', 'blocks','pcwPoly'}
SNR = 40;       % additive Gaussian noise

N = 256;   % signal length
M = round(N/4);    % no. of measurements
T = round(M/3);    % sparsity level

% signal dynamics
cshift = -2;
rshift_h = @(z) (rand-0.5)/20;

% reg. parameter for temporal difference..
lambda_diff = sqrt(1);

% rank-1 update mode
delx_mode = 'mil'; % mil or qr
lambda = 0;
fprintf('Tracking a dynamical signal using reweighting..\n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d.', mType, sType, SNR, N, M, T);
disp(str0);
 
maxsim = 200;
SIM_stack = cell(maxsim,1);
SIM_memory = cell(maxsim,1);

x_vec = zeros(N*maxsim,1);
xh_vec = zeros(N*maxsim,1);

%% Setup signal and system 

% Generate a random signal
in = []; in.type = sType; in.T = T; in.randgen = 0; in.wType = 'daub79'; in.take_fwt = 1;
[x sig wave_struct] = genSignal(N,in);
[val ind] = sort(abs(x),'descend');
ind_pos = ind(val>1e-1);
gamma_orig = ind_pos(1:min(length(ind_pos),M-1));

% Setup dynamical model
x_orig = x; % Initial signal

% At every time instance, add to the original/previous signal 
% an integer circular shift that is known
% a random drift that is unknown
if isempty(wave_struct)
    % cshfits
    F_h = @(x,cshift,rshift) interp1(1:N,circshift(x,cshift),[1:N]+rshift,'linear','extrap')';
    W_h = @(z) z; iW_h = @(z) z; 
else
      
    W_h = wave_struct.W_h;
    iW_h = wave_struct.iW_h;    
    
    F_h = @(x,cshift,rshift) W_h(interp1(1:N,circshift(iW_h(x),cshift),[1:N]+rshift,'linear','extrap')');
end

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

%% initialize...! 
% 
% Assume a previous signal is known.. 
% [val_sort ind_sort] = sort(abs(x),'descend');
% xh_old = x;
% xh_old(ind_sort(M+1:end)) = 0;
% 
% lambda_diff = sqrt(1e-0);
% 
% A = [A; lambda_diff*eye(N)];
% y = [A(1:M,:)*x+e; lambda_diff*xh_old];


% rwt L1 using l1homotopy
in = [];
in.tau = tau; W = tau;
in.delx_mode = delx_mode;
in.debias = 0;
in.verbose = 0;
in.plots = 0;
in.record = 1;
in.err_fun = err_fun;
tic 
for wt_itr = 1:5 
    
    W_old = W;
    
    out = l1homotopy(A,y,in);
    xh = out.x_out;
    iter_bpdn = out.iter;
    time_bpdn = toc;
    gamma_bpdn = out.gamma;
    err_bpdn = out.error_table;
    
    % Update weights
    xh_old = xh;
    
    alpha = 1; epsilon = 1;
    beta = M*(norm(xh_old,2)/norm(xh_old,1))^2;
    W = tau/alpha./(beta*abs(xh_old)+epsilon);
    
    % yh = A*xh_old;
    % u = -sign(xh_old).*W_old; % out.pk;
    % Atdyu = A'*(yh-y)-u;
    
    yh = A*xh_old;
    Atr = A'*(A*xh_old-y);
    u =  -W.*sign(xh_old)-Atr;
    pk_old = Atr+u;
    
    in = out;
    in.xh_old = xh_old;
    in.pk_old = pk_old;
    in.u = u;
    in.W_old = W_old;
    in.W = W;
    % in.Atdyu = Atdyu;
end
W = W_old;
 

% in = [];
% in.tau = tau;
% in.delx_mode = delx_mode;
% in.debias = 0;
% in.verbose = 0;
% in.plots = 0;
% in.record = 1;
% in.err_fun = err_fun;
% tic
% out = l1homotopy(A,y,in);
% xh = out.x_out;
% iter_bpdn = out.iter;
% time_bpdn = toc;
% gamma_bpdn = out.gamma;
% err_bpdn = out.error_table;
% 
% W = tau;
%
% x_sparsa = 0;
% tic;
% [x_sparsa,x_debias_SpaRSA_m,obj_SpaRSA_m_cont,...
%     times_SpaRSA_m_cont,debias_start_SpaRSA_m,mse_SpaRSA_m,taus_SpaRSA_m, numA, numAt]= ...
%     SpaRSA_adpW(y,A,tau,...
%     'Monotone',0,...
%     'adp_wt',0,...
%     'Debias',0,...
%     'StopCriterion',2,...
%     'ToleranceA',1e-4,...
%     'Safeguard',1,...
%     'MaxiterA',maxiter,...
%     'Verbose',0,...
%     'True_x',x,...
%     'Continuation',1,...
%     'Continuationsteps',-1);
% 
% time_sparsa = toc;
% iter_sparsa = (numA+numAt)/2;
% error_sparsa = norm(x-x_sparsa)/norm(x);


for sim = 1:maxsim
    
    % Update the solution after updating the measurement matrix and/or the
    % sparse signal
    x_old = x; sig_old = sig;
    y_old = y; A_old = A(1:M,:);
    
    % Time-varying signal setup
    dyn_ref = 'previous';
    switch dyn_ref
        case 'previous'
            % In this model, x may become dense after some time
            
            % Update the sparse signal.
            rshift = rshift_h(1);
            x = F_h(x_old,cshift, rshift);
            % Threshold the original signal...
            % (to stop it from becoming dense) ???
            % [x_sort ind_sort] = sort(abs(x),'descend');
            % x(ind_sort(M+1:end)) = 0; 
            
            % Update the signal estimate
            xh_old = F_h(xh,cshift,0);
            % xh_old(abs(xh_old)<max(xh_old)/1e3) = 0;
            xh_old(abs(xh_old)<tau/sqrt(log(N))) = 0;
            
        case 'initial'
            % In this model, x will probably remain sparse for all times
            
            % Update the sparse signal.
            cshift = -sim; rshift = rshift_h(1);
            x = F_h(x_orig,cshift, rshift);            
            
            % Update the signal estimate
            xh_old = F_h(xh,cshift-sim+1,0);
            % xh_old(abs(xh_old)<max(xh_old)/1e3) = 0;
            xh_old(abs(xh_old)<tau/sqrt(log(N))) = 0;
    end
    
    % Update weights
    W_old = W;
    
    alpha = 1; epsilon = 1;
    beta = M*(norm(xh_old,2)/norm(xh_old,1))^1;
    W = tau/alpha./(beta*abs(xh_old)+epsilon);
    
    % Update the measurements
    in = []; in.type = mType;
    B = genAmat(M,N,in);
    e = randn(M,1)*sigma;

    % add difference regularization    
    A = [0*A_old+1*B; lambda_diff*eye(N)];
    y = [A(1:M,:)*x+e; lambda_diff*xh_old];
    
    % To remove diff. regularization, set lambda_diff = 0;
    if lambda_diff == 0
        A = A(1:M,:);
        y = y(1:M,:);
    end
    
    % parameter selection (you may change it at every iteration...)
    %         % tau = sigma*sqrt(log(N));
    %         tau = max(1e-3*max(abs(A'*y)),sigma*sqrt(log(N)));
    
    homotopy_mode = 'dummy';
    switch homotopy_mode
        case 'dummy'
            % create a dummy variable...
            % use homotopy on the measurements...
            % in principle, we can start with any xh_old with this formulation
            % and any starting value of tau or W...
            W_old = W;
            gamma = find(xh_old);
            M_trunc = size(A,1); % M-1;
            if length(gamma) > M_trunc
                [xh_sort ind_sort] = sort(abs(xh_old),'descend');
                xh_old(ind_sort(M_trunc+1:end)) = 0;
                gamma = ind_sort(1:M_trunc);            
            end            
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
    
    tic
    in = out;
    gamma_old = gamma;
    in.gamma = gamma_old;
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
    % in.Atdyu = Atdyu;
    in.W = W;
    in.W_old = W_old;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 0;
    in.err_fun = @(z) (norm(x-z)/norm(x))^2;
    out = l1homotopy(A,y,in);
    xh_dynRWT = out.x_out;
    gamma_dynRWT = out.gamma;
    iter_dynRWT = out.iter;
    time_dynRWT = toc;
    
    xh = xh_dynRWT;
    
    %% Check the solution using BPDN directly with l1homotopy
    %     in = []; x_old = x;
    %     in.W = W;
    %     in.delx_mode = delx_mode;
    %     in.debias = 0;
    %     in.verbose = 0;
    %     in.maxiter = maxiter;
    %     in.plots = 0;
    %     in.record = 1;
    %     in.err_fun = err_fun;
    %     tic
    %     out = l1homotopy(A,y,in);
    %     xh = out.x_out;
    %     iter_bpdn = out.iter;
    %     time_bpdn = toc;
    %     gamma_bpdn = out.gamma;
    %     err_bpdn = out.error_table;
    
    %% SpaRSA 
    x_sparsa = xh_old; W_sparsa = W/tau;
    psi_function = @(x,tau) soft(x,tau*W_sparsa);
    phi_function = @(x) sum(abs(W_sparsa.*x));
    tic;
    [x_sparsa,x_debias_SpaRSA_m,obj_SpaRSA_m_cont,...
        times_SpaRSA_m_cont,debias_start_SpaRSA_m,mse_SpaRSA_m,taus_SpaRSA_m, numA, numAt]= ...
        SpaRSA_adpW(y,A,tau,...
        'Monotone',0,...
        'adp_wt',0,...
        'W_new',W_sparsa,...
        'Debias',0,...
        'Initialization',x_sparsa,...
        'StopCriterion',2,...
        'ToleranceA',1e-4,...
        'psi',psi_function,...
        'phi',phi_function,...
        'Safeguard',1,...
        'MaxiterA',maxiter,...
        'Verbose',0,...
        'True_x',x,...
        'Continuation',1,...
        'Continuationsteps',-1);
    
    time_sparsa = toc;
    iter_sparsa = (numA+numAt)/2;    
    error_sparsa = norm(x-x_sparsa)/norm(x);     
    
    %%
    SIM_stack{sim} = [sim, tau, ...
        norm(x-xh_dynRWT)^2/norm(x)^2, sum(iter_dynRWT,2), sum(time_dynRWT,2),...
        norm(x-x_sparsa)^2/norm(x)^2, sum(iter_sparsa,2), sum(time_sparsa,2)];
    
    fprintf('sim %d. tau = %3.4g, (err,iter,time): dynRWT l1homotopy-%3.4g,%3.4g,%3.4g; SpaRSA-%3.4g,%3.4g,%3.4g. \n', ...
        SIM_stack{sim});
    
    %% plot recovered signals
    figure(1); 
    subplot(211); 
    plot([x xh xh_dynRWT]);
    subplot(212); 
    plot([iW_h(x) iW_h(xh)]);
    title('Comparison betweeen the original and reconstructed signal')
    
    %% Reconstructed signal
    x_vec((sim-1)*N+1:sim*N,1) = x;
    xh_vec((sim-1)*N+1:sim*N,1) = xh_dynRWT; 
end
% mS =  mean(cell2mat(SIM_stack),1);
mS =  mean(cell2mat(SIM_stack(ceil(maxsim*0.2):maxsim,:)),1);
fprintf('Average results for last 80 percent iterations: maxsim %d. tau = %3.4g, (err,iter,time): dynRWT l1homotopy-%3.4g,%3.4g,%3.4g; SpaRSA-%3.4g,%3.4g,%3.4g. \n', maxsim, mS(2:end));
