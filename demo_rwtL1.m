% demo_rwtL1
%
% Solves the following reweighted BPDN problem
% min_x  \sum w_i |x_i| + 1/2*||y-Ax||_2^2
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
addpath utils/qr/
addpath utils/utils_Wavelet

disp(['--------------------',datestr(now),'-------------------------'])

% load RandomStates
% 
rseed = 2012;
rand('state',rseed);
randn('state',rseed);

% simulation parameters
mType = 'randn'; % {'randn','orth','rdct'};
sType = 'randn'; % {'randn','sign','highD', 'blocks','pcwPoly'}
SNR = 40;       % additive Gaussian noise

% reweighted setup
rwt = 5;        % number of reweighting iterations

N = 256;   % signal length
M = round(N/2);    % no. of measurements
T = round(M/3);    % sparsity level

% rank-1 update mode 
delx_mode = 'mil'; % mil or qr 
fprintf('Iterative reweighting..\n');
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
        
    %% Iterative reweighting...    
    xh_old = xh; W_old = tau; 
    error_rwt = []; iter_rwt = []; time_rwt = []; 
    for rwt_itr = 1:rwt
        
        % Update weights
        alpha = 1; epsilon = 1;
        beta = M*(norm(xh_old,2)/norm(xh_old,1))^2;
        W = tau/alpha./(beta*abs(xh_old)+epsilon);
        
        homotopy_mode = 'dummy';
        switch homotopy_mode
            case 'dummy'                
                % create a dummy variable...
                % use homotopy on the measurements...
                % in principle, we can start with any xh_old with this formulation
                % and any starting value of W_old (for instance W)...
                W_old = W;
                % yh = A*xh_old;
                % u = -sign(xh_old).*W_old; % out.pk;
                % Atdyu = A'*(yh-y)-u;
                yh = A*xh_old;
                Atr = A'*(A*xh_old-y);
                u =  -W.*sign(xh_old)-Atr;
                pk_old = Atr+u;
            otherwise
                % use homotopy on weights
                u = out.pk;
                Atdyu = 0;
        end
        
        in = out;
        in.xh_old = xh_old;
        in.pk_old = pk_old;
        in.u = u;
        in.W_old = W_old;
        in.W = W;
        % in.Atdyu = Atdyu;
        in.delx_mode = delx_mode;
        in.debias = 0;
        in.verbose = 0;
        in.plots = 0;
        in.record = 1;
        in.err_fun = @(z) (norm(x-z)/norm(x))^2;
        tic
        out = l1homotopy(A,y,in);
        xh_old = out.x_out;
        W_old = W;
        gamma_rwt = out.gamma;
        error_rwt = [error_rwt err_fun(xh_old)];
        iter_rwt = [iter_rwt out.iter];
        time_rwt = [time_rwt toc];
    end
        
%     %% Check the solution using BPDN directly with l1homotopy
%     in = []; x_old = x;
%     in.W = W;
%     in.delx_mode = delx_mode;
%     in.debias = 0;
%     in.verbose = 0;
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
    
    %%
    SIM_stack{sim} = [sim, tau, ...
        norm(x-xh)^2/norm(x)^2, sum(iter_bpdn,2), sum(time_bpdn,2), ...
        norm(x-xh_old)^2/norm(x)^2, sum(iter_rwt,2), sum(time_rwt,2)];  
    
    fprintf('sim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; rwt l1homotopy-%3.4g,%3.4g,%3.4g. \n', ...
        SIM_stack{sim});

    %% plot recovered signals
    figure(1); plot([x xh xh_old]);    
    title('Comparison betweeen the new and old homotopy code')
end
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: maxsim %d. tau = %3.4g, (err,iter,time): l1homotopy from scratch-%3.4g,%3.4g,%3.4g; rwt l1homotopy-%3.4g,%3.4g,%3.4g. \n', maxsim, mS(2:end));
