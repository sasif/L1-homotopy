% BPDN rwt update (initialize with adaptive_wtBPDN)
%
% Solve the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% Initialize with adaptive weighted BPDN
% and dynamically update the weights w_i
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
% Created: April 16, 2011

function out = script_rwtBPDN_adaptive(in)

A = in.A; y = in.y; x = in.x; 
gamma_orig = find(x);
[M N] = size(A);
tau = in.tau; max_rwt = in.max_rwt;
x_init = in.x_init; 

%% parameter selection
rwt_mode = in.rwt_mode;

delx_mode = in.delx_mode;
homotopy_update = 'v1';
debias = in.debias;
verbose = in.verbose;

maxiter = 2*N;
iter_ALL = [];
err_ALL = [];
time_ALL = [];
supp_diff = [];

if norm(x_init) == 0
    %% Adaptive support + weight selection
    % Fast shrinkage on the active set...???
    %
    % use a large value for Tsteps (~10-100)
    %     ewt = 10; shrinkage_mode = 'Tsteps'; shrinkage_flag = 0;
    %
    % use a large value for Trwt
    %     ewt = 100; shrinkage_mode = 'Trwt'; shrinkage_flag = 0; % Gaussian
    % use shrinkage_flag = 0 for Gaussian signal... (amplitude variations)
    % use shrinkage_flag = 2 for sign/ones signal... (flat)
    %
    
    % 'rwt' works better when nonzero components have diverse amplitudes
    % ewt = 2; shrinkage_mode = 'rwt'; shrinkage_flag = 0;
    ewt = 1; shrinkage_mode = 'Trwt'; shrinkage_flag = 0;
    % ewt = 2; shrinkage_mode = 'Tsteps'; shrinkage_flag = 3;
    % ewt = 1; shrinkage_mode = 'OLS'; shrinkage_flag = 0;
    
    if verbose
        fprintf(' ewt = %1.2g, shrinkage_mode: %s, ',ewt, shrinkage_mode);
    end
    
    in = [];
    in.tau = tau;
    in.ewt = ewt; % Setting ewt = 1 with Tsteps solves unweighted LASSO
    in.shrinkage_flag = shrinkage_flag;
    in.shrinkage_mode = shrinkage_mode; % (either fixed or adaptive weighting)
    in.maxiter = maxiter;
    in.debias = debias;
    in.early_terminate = 0;
    in.x_orig = x;
    in.record = 1;
    in.omp = 0;
    in.plots = 0;
    in.plot_wts = 0;
    in.delx_mode = delx_mode;
    % in.Te = T;
    tic;
    out = wtBPDN_adaptive_function(A, y, in);
    time_ALL = [time_ALL toc];
    % time_ALL = [time_ALL out.time];
    x_init = out.x_out;
    gamma_old = out.gamma;
    iter_old = out.iter;
    tau_vec = out.tau_vec;
    iter_ALL = [iter_ALL iter_old];
    err_ALL = [err_ALL norm(x-x_init)/norm(x)];
    supp_diff = [supp_diff length(setxor(gamma_old,gamma_orig))];
    if verbose
        fprintf('err = %3.4g, iter = %3.4g, time = %3.4g \n',err_ALL(end), iter_ALL(end), time_ALL(end));
    end
end
W_new = tau_vec; % These weights are chosen adaptively by the homotopy solver

% THERE SHOULDN'T BE ANY NEED FOR REWEIGHTING HERE
% But just in case...

%% Iterative reweighting
xh_mod = x_init;
gamma_xh = gamma_old;
for rwt_itr = 1:max_rwt
    gamma_old = gamma_xh;
    x_old = xh_mod;
    W_old = W_new;
    
    [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_old,M);    
    W_new = tau/alpha./(beta*abs(x_old)+epsilon);
    
    % W_new = tau/epsilon*ones(N,1);
    % W_new(gamma_old) = min([tau*ones(length(gamma_old),1) tau./(beta*abs(x_old(gamma_old)))],[],2);
    
    switch homotopy_update
        case 'v1'
            % The following Gram matrix and its inverse can be used from the
            % previous homotopy. Too lazy to include that right now...
            % wt BPDN homotopy update
            pk_old = A'*(A*x_old-y);
            pk_old(gamma_old) = sign(pk_old(gamma_old)).*W_old(gamma_old);
            in = [];
            in.x_old = x_old;
            in.gamma = gamma_old;
            in.pk_old = pk_old;
            in.W_old = W_old;
            in.W_new = W_new;
            in.maxiter = maxiter;
            dW = W_new-W_old;
            in.maxiter = maxiter;
            
            % delx = -AtAgx\(-dW(gamma_old).*sign(pk_old(gamma_old)));
            % in.delx = delx;
            in.delx_mode = delx_mode;
            switch in.delx_mode
                case 'mil';
                    % in.delx_mode = 'mil';
                    AtAgx = A(:,gamma_old)'*A(:,gamma_old);
                    iAtAgx = pinv(AtAgx);
                    in.AtA = AtAgx;
                    in.iAtA = iAtAgx;                    
                case 'qr';
                    % in.delx_mode = 'qr';
                    [Q R] = qr(A(:,gamma_old),0);
                    in.Q = Q; in.R = R;
                case 'qrM'
                    % in.delx_mode = 'qrM';
                    [Q0 R0] = qr(A(:,gamma_old));
                    in.Q0 = Q0; in.R0 = R0;
            end
            tic;
            out = wtBPDN_Update_function_v1(A, y, in);
            % time_update = out.time;
            time_update = toc;
            xh_mod = out.x_out;
            gamma_xh = out.gamma;
            iter_update = out.iter;
        case 'v2'
            
            % Homotopy update v2 (weighting embedded in the matrix A)
            AW_old = A*diag(tau./W_old);
            u0_hat = x_old.*(W_old/tau);
            ds = AW_old'*(AW_old*u0_hat-y);
            AW = A*diag(tau./W_new);
            yhat = AW*u0_hat;
            pk_old = ds; pk_old(gamma_old) = sign(pk_old(gamma_old)).*tau;
            in = [];
            in.x_old = u0_hat;
            in.gamma = gamma_old;
            in.pk_old = pk_old;
            in.tau = tau;
            in.maxiter = maxiter;
            in.yhat = yhat;
            
            % ds = pk_old; d = AW'*(yhat-y);
            % delx = AtAgx\(ds(gamma_old)-d(gamma_old));
            % in.delx = delx;
            in.delx_mode = delx_mode;
            switch in.delx_mode
                case 'mil';
                    % in.delx_mode = 'mil';
                    AtAgx = AW(:,gamma_old)'*AW(:,gamma_old);
                    iAtAgx = pinv(AtAgx);                    
                    in.AtA = AtAgx;
                    in.iAtA = iAtAgx;                    
                case 'qr';
                    % in.delx_mode = 'qr';
                    [Q R] = qr(AW(:,gamma_old),0);
                    in.Q = Q; in.R = R;
                case 'qrM'
                    % in.delx_mode = 'qrM';
                    [Q0 R0] = qr(AW(:,gamma_old));
                    in.Q0 = Q0; in.R0 = R0;
                case 'cg'
                    % in.delx_mode = 'cg';
                    in.W_new = tau./W_new;
            end
            tic;
            out = wtBPDN_Update_function_v2(AW, y, in);
            % time_update = out.time;
            time_update = toc;

            xh_mod = out.x_out.*(tau./W_new);
            gamma_xh = out.gamma;
            iter_update = out.iter;
    end
    iter_ALL = [iter_ALL iter_update];
    err_ALL = [err_ALL norm(x-xh_mod)/norm(x)];
    time_ALL = [time_ALL time_update];
    supp_diff = [supp_diff length(setxor(gamma_xh,gamma_orig))];
end

out.x_out = xh_mod;
out.x_init = x_init;
out.iter = iter_ALL;
out.gamma = gamma_xh;
out.W_new = W_new;
out.err = err_ALL;
out.time = time_ALL;
out.supp_diff = supp_diff;