% BPDN rwt update (initialize with standard homotopy)
%
% Solve the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% Initialize with unweighted BPDN
% and dynamically update the weights w_i
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
% Created: March 16, 2011

function out = script_rwtBPDN_iterative(in);
A = in.A; y = in.y; x = in.x;
gamma_orig = find(x);
[M N] = size(A);
tau = in.tau; max_rwt = in.max_rwt;
x_init = in.x_init;

%% parameter selection
% max_rwt = 3;
% lambda = 1e-2;
% tau = min(lambda*max(abs(A'*y)),sigma*sqrt(log(N)));

rwt_mode = in.rwt_mode; 

delx_mode = in.delx_mode;
homotopy_update = 'v1';

maxiter = 2*N;
iter_ALL = [];
err_ALL = [];
time_ALL = [];
supp_diff = [];

if norm(x_init) == 0
    %% Standard BPDN
    in = [];
    in.tau = tau;
    in.maxiter = maxiter;
    in.x_orig = x;
    in.record = 1;
    in.delx_mode = delx_mode;
    tic;
    out = BPDN_homotopy_function(A, y, in); %BPDN
    % time_ALL = [time_ALL out.time];
    time_ALL = [time_ALL toc];
    
    x_init = out.x_out;
    gamma_old = out.gamma;
    iter_old = out.iter;
    iter_ALL = [iter_ALL iter_old];
    err_ALL = [err_ALL norm(x-x_init)/norm(x)];
    
    supp_diff = [supp_diff length(setxor(gamma_old,gamma_orig))];
end
W_new = ones(N,1)*tau; % homotopy solved with same fixed weights. 

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
            pk_old = A'*(A*x_old-y);
            pk_old(gamma_old) = sign(pk_old(gamma_old)).*W_old(gamma_old);
            in = [];
            in.x_old = x_old;
            in.gamma = gamma_old;
            in.pk_old = pk_old;
            in.W_old = W_old;
            in.W_new = W_new;
            dW = W_new-W_old;
            in.maxiter = maxiter;
            
            % delx = -AtAgx\(-dW(gamma_old).*sign(pk_old(gamma_old)));
            % in.delx = delx;
            in.delx_mode = delx_mode;
            switch in.delx_mode
                case 'mil';
                    % in.delx_mode = 'mil';
                    % The following gram matrix and its inverse can be used from the
                    % previous homotopy. Too lazy to include that right now...
                    % wt BPDN homotopy update
                    AtAgx = A(:,gamma_old)'*A(:,gamma_old);
                    iAtAgx = inv(AtAgx);
                    in.AtA = AtAgx;
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
            tic
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
                    iAtAgx = inv(AtAgx);                    
                    in.AtA = AtAgx;
                    in.iAtA = iAtAgx;
                case {'qr','chol'};
                    % in.delx_mode = 'qr';
                    [Q R] = qr(AW(:,gamma_old),0);
                    in.Q = Q; in.R = R;
                case 'qrM'
                    % in.delx_mode = 'qrM';
                    [Q0 R0] = qr(AW(:,gamma_old));
                    in.Q0 = Q0; in.R0 = R0;
            end
            tic
            out = wtBPDN_Update_function_v2(AW, y, in);
            time_update = out.time;
            time_update = toc;
            xh_mod = out.x_out.*(tau./W_new);
            gamma_xh = out.gamma;
            iter_update = out.iter;
    end
    iter_ALL = [iter_ALL iter_update];
    err_ALL = [err_ALL norm(x-xh_mod)/norm(x)];
    time_ALL = [time_ALL time_update];
    supp_diff = [supp_diff length(setxor(gamma_xh,gamma_orig))];    
	
    %     % compare with solution from scratch
    %     in = [];
    %     in.tau = tau;
    %     in.maxiter = maxiter;
    %     in.x_orig = x;
    %     in.record = 1;
    %     AW = A*diag(tau./W_new);
    %     out_new = BPDN_homotopy_function(AW, y, in); %BPDN
    %     xh_mod = out_new.x_out.*(tau./W_new);
    %     gamma_new = out_new.gamma;
    %     new_iter = out_new.iter;
    %     new_time = out_new.time;
    
    %     figure(1);
    %     subplot(221); plot([x xh xh_v2 xh_mod]);
    %     subplot(222); plot([xh-xh_v2 xh-xh_mod]);
    %     subplot(223); plot(xh,x,'.');
    %     subplot(224); plot([W_new W_old]);
    %     err_rec = norm(xh-x)/norm(x);
    %     err_old = norm(x_init-x)/norm(x);
    % disp(sprintf(['RWT iter. = %d, wtBPDN-LASSO homotopy iter. = %d, Update = %d, Update-v2=%d, Rec. error = %3.4g, improve = %3.4g'], rwt_itr, new_iter, iter_update, update_iter_v2, err_rec , (err_old-err_rec)/err_old*100));
end

out.x_out = xh_mod;
out.x_init = x_init;
out.iter = iter_ALL;
out.gamma = gamma_xh;
out.W_new = W_new;
out.err = err_ALL;
out.time = time_ALL;
out.supp_diff = supp_diff;