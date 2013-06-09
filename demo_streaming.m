% demo_streaming
%
% Solves the following streaming problem
% min_x  \sum_i \|x_i-1\|_1 + 1/2*||B_i-1 x_i-1 + A_i x_i - y_i||_2^2
%
% by adding and removing one set of measurements at every time instant 
%
% Applications:
%   	streaming signal recovery from overlapping measurements
%
% We can also add any other regularization operator in the reconstruction
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: August 2012

clear
close all % force

% Limit the number of computational threads (for profiling)
% maxNumCompThreads(1);

%% Setup path
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath utils/
addpath utils/qr/
addpath solvers

disp(['--------------------',datestr(now),'-------------------------'])

% load RandomStates
% 
rseed = 2012;
rand('state',rseed);
randn('state',rseed);

% simulation parameters
mType = 'sign'; % {'randn','orth','rdct'};
sType = 'pcwreg'; % {'randn','sign','highD', 'blocks','pcwPoly'}
SNR = inf;       % additive Gaussian noise

N = 128;   % signal length
M = round(N/4);    % no. of measurements
T = round(M/4);    % sparsity level

% Overlapping measurements  
LM = 2*N; 

% streaming window
P = 4; % size of the working window is P*N       

% remove the top-left edge of the system matrix
LeftEdge_trunc = 1;
% fprintf('Consider oldest set of coefficients to be fully known, and remove their contribution from the measurements... \n');
    
% rank-1 update mode
delx_mode = 'mil'; % mil or qr 
fprintf('Streaming measurements with a time-varying signal ... \n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d, P = %d.', mType, sType, SNR, N, M, T, P);
disp(str0);

% Simulation parameters
maxsim = 128;
SIM_stack = cell(maxsim,1);
SIM_memory = cell(maxsim,1);
 
x_vec = zeros(N*maxsim,1);
xh_vec = zeros(N*maxsim,1);
sig_vec = zeros(N*maxsim,1);
sigh_vec = zeros(N*maxsim,1);

%% Signal and system Setup 
% Generate long sequence of signal and sensing matrices

% Generate a random signal
in = []; in.type = sType; in.T = maxsim*T; in.randgen = 0; 
in.wType = 'daub79'; in.take_fwt = 1; in.J = 5; %log2(N)-3;
genSignal_h = @(N) genSignal(N,in);
[x_long sig wave_struct] = genSignal(maxsim*N,in);

if isempty(wave_struct)
    W_h = @(z) z; iW_h = @(z) z; 
else 
    W_h = wave_struct.W_h;
    iW_h = wave_struct.iW_h;    
end

in = []; in.type = mType;
genAmat_h = @(M,N) genAmat(M,N,in);

x = zeros(N*P+LM-N,1);
A = zeros(M*P,N*P+LM-N);
% A = sparse(M*P,N*(P+1));
y = zeros(M*P,1);
 
% Overlapping measurements (type imagesc(At) to view the structure) 
% LM = 2*N... 
s_ind = 0;
% xt = [genSignal_h(N); genSignal_h(N)];
xt = [W_h(sig(1:N)); W_h(sig(N+1:2*N))];
s_ind = s_ind+2*N;

At = genAmat_h(M,LM);
sigma = sqrt(norm(At*xt)^2/10^(SNR/10)/M);
et = randn(M,1)*sigma;
yt = At*xt+et;

x(1:2*N) = xt;
A(1:M,1:2*N) = At;
y(1:M) = yt;

for ii = 2:P
    % xt = genSignal_h(N);   
    xt = W_h(sig(s_ind+1:s_ind+N));
    s_ind = s_ind+N;
    
    x(ii*N+1:(ii+1)*N) = xt;
    
    At = genAmat_h(M,2*N);
    A((ii-1)*M+1:ii*M,(ii-1)*N+1:(ii+1)*N) = At;
    
    et = randn(M,1)*sigma;
    y((ii-1)*M+1:ii*M) = At*x((ii-1)*N+1:(ii+1)*N)+et;
end

% parameter selection 
% tau = sigma*sqrt(log(N));
tau = max(1e-2*max(abs(A'*y)),sigma*sqrt(log(2*N)));


maxiter = 2*P*N;
err_fun = @(z) (norm(x-z)/norm(x))^2;

%% Initialize by solving a rwt L1 problem
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
    
    W_old = W;
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

for sim = 1:maxsim
    
    %% Update the solution after updating the measurement matrix and/or the
    % sparse signal
    x_old = x;  
    y_old = y; A_old = A;
    
    % Time-varying signal setup
    % Shift the window (apply any prediction?)
    x(1:end-N) = x(N+1:end); 
    xt_old = x(end-N+1:end);

    % xt = genSignal_h(N); 
    if s_ind+N > length(sig)
        break;
    end
    xt = W_h(sig(s_ind+1:s_ind+N));
    s_ind = s_ind+N;
    x(end-N+1:end) = xt; 
    
    % System matrix setup...
    % Shift up and left
    A(1:(P-1)*M,1:P*N) = A(M+1:end,N+1:end);     
    % new measurement matrix
    At = genAmat_h(M,2*N); 
    A(end-M+1:end,end-2*N+1:end) = At;
    
    % shift old measurements and add one new set of measurementts
    y(1:end-M) = y(M+1:end);
    e = randn(M,1)*sigma;
    yt = At*x(end-2*N+1:end)+e;
    y(end-M+1:end) = yt;
    
    % Update the signal estimate (for warm start)
    xh_old = xh;
    xh_old(1:end-N) = xh(N+1:end);
    xh_t = xh(end-N+1:end);
    xh_old(end-N+1:end) = 0*xh_t;   
    xh_old(abs(xh_old)<tau/sqrt(log(P*N))) = 0;
    
    % Remove the top-left edge of the system matrix
    if LeftEdge_trunc
        A_old = A; y_old = y;
        
        alpha0h = xh_old(1:N);
        xh_old = xh_old(N+1:end);
        y = y-A(:,1:N)*alpha0h;
        
        A = A(:,N+1:end);
        alpha0 = x(1:N);
        x = x(N+1:end);
    end
    
    % Update weights 
    alpha = 1; epsilon = 1;
    beta = M*(norm(xh_old,2)/norm(xh_old,1))^1;
    W = tau/alpha./(beta*abs(xh_old)+epsilon);
    % W(end-N+1:1) = tau;
    W_old = W;
    
    homotopy_mode = 'dummy';
    switch homotopy_mode
        case 'dummy'
            % create a dummy variable...
            % use homotopy on the measurements...
            % in principle, we can start with any xh_old with this formulation
            % and any starting value of tau or W...
            gamma = find(xh_old);
            M_trunc = size(A,1); % P*(M-1);
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
    in.record = 1;
    in.err_fun = @(z) (norm(x-z)/norm(x))^2;
    tic
    out = l1homotopy(A,y,in);
    time_streamingRWT = toc;
    xh_streamingRWT = out.x_out;
    gamma_streamingRWT = out.gamma;
    iter_streamingRWT = out.iter;
    
    xh = xh_streamingRWT;
    
    % Check the solution using BPDN directly with l1homotopy
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
        norm(x-xh_streamingRWT)^2/norm(x)^2, sum(iter_streamingRWT,2), sum(time_streamingRWT,2), ... 
        norm(x-x_sparsa)^2/norm(x)^2, sum(iter_sparsa,2), sum(time_sparsa,2)];
    
    % l1homotopy-%3.4g,%3.4g,%3.4g; 
    % norm(x-xh)^2/norm(x)^2, sum(iter_bpdn,2), sum(time_bpdn,2), ...        
    
%     fprintf('sim %d. tau = %3.4g, (err,iter,time): streamingRWT homotopy-%3.4g,%3.4g,%3.4g, SpaRSA-%3.4g,%3.4g,%3.4g. \n', ...
%         SIM_stack{sim});
%     
    %% plot sparse coeffs. 
    fig(1);     
    subplot(211);
    plot([x xh_streamingRWT x_sparsa]); 
    title('Comparison betweeen the original and reconstructed signal')
    drawnow;
    
    %% Reconstructed signal
    if LeftEdge_trunc
        xh = [alpha0h; xh_streamingRWT];
        x = [alpha0; x];
        A = A_old; y = y_old;
    else
        xh = xh_streamingRWT;
    end
    x_vec((sim-1)*N+1:sim*N,1) = x(1:N);
    xh_vec((sim-1)*N+1:sim*N,1) = xh(1:N); 
    sig_vec((sim-1)*N+1:sim*N,1) = iW_h(x(1:N));
    sigh_vec((sim-1)*N+1:sim*N,1) = iW_h(xh(1:N)); 
    
    %% plot recovered signals
    subplot(212);
    plot([sig_vec(1:s_ind(end)) sigh_vec(1:s_ind(end))]); 
        
end
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: maxsim %d. tau = %3.4g, (err,iter,time): streamingRWT l1homotopy-%3.4g,%3.4g,%3.4g; SpaRSA-%3.4g,%3.4g,%3.4g. \n', maxsim, mS(2:end));
% l1homotopy-%3.4g,%3.4g,%3.4g; 

% l1homotopy-%3.4g,%3.4g,%3.4g;
L = N;
err_l1homotopy = norm(sig_vec(L-N+1:N*sim)-sigh_vec(L-N+1:N*sim,1))^2/norm(sig_vec(L-N+1:N*sim))^2; 
fprintf('Signal error: l1homotopy-%3.4g. \n',err_l1homotopy);

