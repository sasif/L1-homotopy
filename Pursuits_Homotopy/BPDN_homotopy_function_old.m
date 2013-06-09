% BPDN_homotopy_function.m
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \tau ||x||_1 + 1/2*||y-Ax||_2^2
%
% Inputs:
% A - m x n measurement matrix
% y - measurement vector
% tau - final value of regularization parameter
% maxiter - maximum number of homotopy iterations
%
% Outputs:
% x_out - output for BPDN
% gamma_x - support of the solution
% total_iter - number of homotopy iterations taken by the solver
% total_time - time taken by the solver
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
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
x_orig = in.x_orig;
err_record = in.record;

% Initialization of primal and dual sign and support
z_x = zeros(N,1);
gamma_x = [];       % Primal support

% Initial step
pk_old = -A'*y;
[c i] = max(abs(pk_old));

gamma_xk = i;

epsilon = c;
xk_1 = zeros(N,1);

z_x(gamma_xk) = -sign(pk_old(gamma_xk));
pk_old(gamma_xk) = sign(pk_old(gamma_xk))*epsilon;

z_xk = z_x;

% loop parameters
done = 0;
iter = 0;
data_precision = eps;   % floating point precision

old_delta = 0;
out_x = [];
count_delta_stop = 0;

constraint_plots = 1;

AtAgx = A(:,gamma_xk)'*A(:,gamma_xk);
iAtAgx = inv(A(:,gamma_xk)'*A(:,gamma_xk));

G = @(z) A*z;
Gt = @(z) A'*z;

if err_record
    error_table = [epsilon norm(xk_1-x_orig) 1];
end

while iter < maxiter
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')
    
    gamma_x = gamma_xk;
    z_x = z_xk;
    x_k = xk_1;
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%% update on x %%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Update direction
    % delx = inv(A(:,gamma_x)'*A(:,gamma_x))*z_x(gamma_x);
    %     if iter > 1
    %     [delx_cg, res, iter] = cgsolve2(delx_vec(gamma_x), AtAgx, z_x(gamma_x), 1e-8, 100, 1)
    %     end
    delx = iAtAgx*z_x(gamma_x);
    % delx = (A(:,gamma_x)'*A(:,gamma_x))\z_x(gamma_x);
    delx_vec = zeros(N,1);
    delx_vec(gamma_x) = delx;
    %delx_vec = sparse(delx_vec);
    
    pk = pk_old;
    %dk = A'*(A*delx_vec);
    % Agdelx = A(:,gamma_x)*delx;
    % dk = A'*Agdelx;
    % dk = A'*(A*delx_vec);
    dk_t = G(delx_vec);
    dk = Gt(dk_t);
    
    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW.
    pk_temp = pk_old;
    gammaL_temp = find(abs(abs(pk_old)-epsilon)<min(epsilon,2*eps));
    %     pk_temp(gammaL_temp) = sign(pk_old(gammaL_temp))*epsilon;
    
    xk_temp = x_k;
    gammaX_temp = find(abs(x_k)<1*eps);
    %     xk_temp(gammaX_temp) = 0;
    %%%---
    
    % Compute the step size
    [idelta, out_x, delta, chk_x] = update_primal(gamma_x, gamma_x, z_x,  xk_temp, delx_vec, pk_temp, dk, epsilon, out_x);
    %     disp(sprintf(['iter = %d, delta = %3.4g, idelta = %d, flag = %d'], iter, delta, idelta, chk_x));
    %     pause(1/60)
    if old_delta < 4*eps && delta < 4*eps
        count_delta_stop = count_delta_stop + 1;
    else
        count_delta_stop = 0;
    end
    if count_delta_stop >= 500
        disp('stuck somewhere');
        break;
    end
    old_delta = delta;
    
    xk_1 = x_k+delta*delx_vec;
    pk_old = pk+delta*dk;
    epsilon_old = epsilon;
    epsilon = epsilon-delta;
    
    %% Check convergence criterion
    %     if length(gamma_x) < M/2
    %         xhat = zeros(N,1);
    %         xhat(gamma_x) = AtAgx\(A(:,gamma_x)'*y);
    %         if norm(y-A*xhat) < tau
    %             xk_1 = xhat;
    %             break;
    %         end
    %     end
    
    if epsilon <= tau;
        xk_1 = x_k + (epsilon_old-tau)*delx_vec;
        total_time= cputime-t0;
        break;
    end
    if length(gamma_x)-chk_x >= Te
        total_time = cputime-t0;
        break;
    end
    if chk_x == 1
        % If an element is removed from gamma_x
        gx_old = gamma_x;
        len_gamma = length(gamma_x);
        
        outx_index = find(gamma_x==out_x);
        gamma_x(outx_index) = gamma_x(len_gamma);
        gamma_x(len_gamma) = out_x;
        gamma_xk = gamma_x(1:len_gamma-1);
        
        rowi = outx_index; % ith row of A is swapped with last row (out_x)
        colj = outx_index; % jth column of A is swapped with last column (out_lambda)
        AtAgx_ij = AtAgx;
        temp_row = AtAgx_ij(rowi,:);
        AtAgx_ij(rowi,:) = AtAgx_ij(len_gamma,:);
        AtAgx_ij(len_gamma,:) = temp_row;
        temp_col = AtAgx_ij(:,colj);
        AtAgx_ij(:,colj) = AtAgx_ij(:,len_gamma);
        AtAgx_ij(:,len_gamma) = temp_col;
        iAtAgx_ij = iAtAgx;
        temp_row = iAtAgx_ij(colj,:);
        iAtAgx_ij(colj,:) = iAtAgx_ij(len_gamma,:);
        iAtAgx_ij(len_gamma,:) = temp_row;
        temp_col = iAtAgx_ij(:,rowi);
        iAtAgx_ij(:,rowi) = iAtAgx_ij(:,len_gamma);
        iAtAgx_ij(:,len_gamma) = temp_col;
        
        AtAgx = AtAgx_ij(1:len_gamma-1,1:len_gamma-1);
        iAtAgx = update_inverse(AtAgx_ij, iAtAgx_ij,2);
        xk_1(out_x) = 0;
    else
        % If an element is added to gamma_x
        gamma_xk = [gamma_x; idelta];
        new_x = idelta;
        
        AtgxAnx = A(:,gamma_x)'*A(:,new_x);
        AtAgx_mod = [AtAgx AtgxAnx; AtgxAnx' A(:,new_x)'*A(:,idelta)];
        
        AtAgx = AtAgx_mod;
        iAtAgx = update_inverse(AtAgx, iAtAgx,1);
        xk_1(idelta) = 0;
        gamma_x = gamma_xk;
    end
    
    z_xk = zeros(N,1);
    z_xk(gamma_xk) = -sign(pk_old(gamma_xk));
    pk_old([gamma_x]) = sign(pk_old([gamma_x]))*epsilon;
    %     figure(1); plot(pk_old)
    
    if err_record
        error_table = [error_table; epsilon norm(xk_1-x_orig) length(gamma_x)];
    end
end
if err_record
    error_table = [error_table; epsilon norm(xk_1-x_orig) length(gamma_x)];
end
total_iter = iter;
total_time = cputime-t0;

out = [];
out.x_out = xk_1;
out.gamma = gamma_x;
out.iter = total_iter;
out.time = total_time;
out.error_table = error_table;