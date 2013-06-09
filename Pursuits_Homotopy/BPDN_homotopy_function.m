% BPDN_homotopy_function.m
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \tau ||x||_1 + 1/2*||y-Ax||_2^2
%
% Inputs:
% A - m x n measurement matrix
% y - measurement vector
% in - input parameters 
%         tau - final value of regularization parameter
%         maxiter - maximum number of homotopy iterations
%         x_orig - original signal for error computation
%         record - for recording history
%         delx_mode - rank-one update using matrix inversion lemma (mil) or qr factorization (qr)
%         Te - maximum number of non-zero elements allowed
%
% Outputs: 
%   out
%     x_out - output for BPDN
%     gamma - support of the solution
%     iter - number of homotopy iterations taken by the solver
%     time - time taken by the solver
%
% Example usage: 
%     in = [];
%     in.tau = tau; % regularization parameter
%     in.maxiter = maxiter; % maximum allowed iterations
%     in.x_orig = x; % original signal for error computation
%     in.record = 1; % to record history
%     in.delx_mode = 'mil'; % rank-one update using matrix inversion lemma (mil) or qr factorization (qr)
%     in.Te = Te; % maximum number of non-zero elements allowed
%     out = BPDN_homotopy_function(A, y, in); % solve BPDN using matrix A and measurements y
%     xh = out.x_out; % output solution
%     gamma = out.gamma; % support of the solution
%     iter = out.iter; % number of homotopy iterations used
%     time = out.time; % computation time
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
%
% Modified
% May 2012: Added qr update feature
%
%-------------------------------------------+
% Copyright (c) 2007.  Muhammad Salman Asif
%-------------------------------------------+

function out = BPDN_homotopy_function(A, y, in)

% fprintf('Rewrite this function to add qrupdate \n');

t0 = cputime;
N = size(A,2);
M = size(A,1);

tau = in.tau;
Te = inf;
if isfield(in,'Te')
    Te = in.Te;
end
maxiter = 4*N;
if isfield(in,'maxiter');
    maxiter = in.maxiter;
end
early_terminate = 0;
if isfield(in,'early_terminate')
    early_terminate = in.early_terminate;
end
x_orig = in.x_orig;
err_record = in.record;

% Initialization of sign and support
xk_1 = zeros(N,1);
pk_old = -A'*y;
[c idelta] = max(abs(pk_old));

gamma_xh = idelta;
temp_gamma = zeros(N,1);
temp_gamma(gamma_xh) = gamma_xh;
gamma_xc = find([1:N]' ~= temp_gamma);

epsilon = c;
z_x = zeros(N,1);
z_x(gamma_xh) = -sign(pk_old(gamma_xh));
pk_old(gamma_xh) = sign(pk_old(gamma_xh))*epsilon;

% loop parameters
done = 0;
iter = 0;
itr_history = [];
error_table = [];
obj = [];

% initialize delx
in_delx = [];
delx_mode = in.delx_mode;
update_mode = 'init0';
rhs = z_x;
update_delx;

if err_record
    error_table = [epsilon norm(xk_1-x_orig)/norm(x_orig) 1];
    % obj = 1/2*norm(A*xk_1-y)^2+epsilon*norm(xk_1,1);
end

while iter < maxiter
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')
    
    x_k = xk_1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% update support of x %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update direction
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;
    
    if sign(delx_vec(idelta)) == sign(pk_old(idelta)) && iter > 1
        delta = 0; flag = 0;
    else
        pk = pk_old;
        % dk = A'*(A(:,gamma_xh)*delx);
        dk = A'*(A*delx_vec);
        
        %%%--- compute step size
        in = [];
        
        % Setting shrinkage_flag to zero shrinks new active constraint towards the
        % final value instantly if doing so doesn't disturb the active set
        in.delta_flag = 2;
        in.pk = pk; in.dk = dk;
        in.ak = epsilon; in.bk = -1;
        in.gamma = gamma_xh; in.gamma_c = gamma_xc;
        in.delx_vec = delx_vec; in.x = xk_1;
        out = compute_delta(in);
        delta = out.delta; idelta = out.idelta;
        flag = out.flag;
    end
    
    xk_1 = x_k+delta*delx_vec;
    pk_old = pk+delta*dk;
    epsilon_old = epsilon;
    epsilon = epsilon-delta;
    
    if epsilon <= tau;
        xk_1 = x_k + (epsilon_old-tau)*delx_vec;
        total_time= cputime-t0;
        break;
    end
    if length(gamma_xh) >= Te
        total_time = cputime-t0;
        break;
    end
    
    % disp(sprintf(['iter = %d, delta = %3.4g, idelta = %d, flag = %d'], iter, delta, idelta, flag));
    itr_history = [itr_history; idelta delta flag];
    
    %% Check convergence criterion
    if early_terminate
        if length(gamma_x) < M/2
            xhat = zeros(N,1);
            rhs = A(:,gamma_xh)'*y;
            update_mode = 'recompute';
            delx_update;
            xhat(gamma_x) = out;
            if norm(y-A*xhat) < tau
                xk_1 = xhat;
                break;
            end
        end
    end
    
    % update support
    update_supp;
    
    pk_old([gamma_xh; idelta]) = sign(pk_old([gamma_xh; idelta]))*epsilon;
    
    % update delx
    z_x = -sign(pk_old);
    rhs = z_x;
    update_mode = 'update';
    update_delx;
    
    
    constr_violation = nnz((abs(pk_old(gamma_xc))-epsilon)>1e-10);
    sign_violation = nnz((sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh)))>1);
    if constr_violation
        chk = gamma_xc((abs(pk_old(gamma_xc))-epsilon)>1e-10);
        stp = 1;
        fprintf('problem... with constraint violation -- %s \n', mfilename);
        fprintf('Refactorize the matrix... recompute delx \n');
        % some times it comes here due to bad conditioning of AtAgx.
        update_mode = 'init0';
        update_delx;
    end
    if sign_violation>=1
        chk = gamma_xh(sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh))>1);
        stp = 1;
        fprintf('problem... sign mismatch -- %s \n',mfilename);
        fprintf('Refactorize the matrix... recompute delx \n');
        update_mode = 'init0';
        update_delx;
    end
    
    if err_record
        error_table = [error_table; epsilon norm(xk_1-x_orig)/norm(x_orig) length(gamma_xh)];
        % obj = [obj; 1/2*norm(A*xk_1-y)^2+epsilon*norm(xk_1,1)];
    end
end
if err_record
    error_table = [error_table; epsilon norm(xk_1-x_orig)/norm(x_orig) length(gamma_xh)];
    % obj = [obj; 1/2*norm(A*xk_1-y)^2+epsilon*norm(xk_1,1)];
end
total_iter = iter;
total_time = cputime-t0;

out = [];
out.x_out = xk_1;
out.gamma = gamma_xh;
out.iter = total_iter;
out.time = total_time;
out.error_table = error_table;
% out.obj = obj;