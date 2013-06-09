% Iterative reweighting updates for BPDN
%
% v1 (each weight treated as a separate regularization parameter...)
% This version is reported in the paper.
% 
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% and dynamically update the weights w_i
%
% Inputs:
% A - m x n measurement matrix
% y - measurement vector
% in - input structure
%   x_old - old estimate
%   gamma - old signal support
%   AtA - A(:,gamma_xh)'*A(:,gamma_xh);
%   iAtA - inv(AtA)
%	or QR/Cholesky factors
%   pk_old - old constrtaints set
%   W_old - old values of regularization parameter
%   W_new - new values of regularization parameter
%   maxiter - maximum number of homotopy iterations
%
% Outputs:
% out - output structure
%   x_out - output for BPDN
%   gamma - support of the solution
%   iter - number of homotopy iterations taken by the solver
%   time - time taken by the solver
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: March 16, 2011

function out = wtBPDN_Update_function_v1(A, y, in);

% x_old, gamma_xh, AtAgx, iAtAgx, pk_old, W_old, W_new, maxiter
x_old = in.x_old;
gamma_xh = in.gamma;
pk_old = in.pk_old;
W_old = in.W_old;
W_new = in.W_new;
maxiter = in.maxiter;

t0 = cputime;
N = size(A,2);
M = size(A,1);

% pk_old = A'*(A*x_k-y);
% pk_old(gamma_xh) = sign(pk_old(gamma_xh)).*W_old(gamma_xh);
xk_1 = x_old;

epsilon = 0;
Wk = (1-epsilon)*W_old + epsilon*W_new;
dW = W_new-W_old;
z_x = -sign(pk_old);

done = 0;
iter = 0;
idelta = gamma_xh(1); flag = 1;
itr_history = [];

% initialize delx
in_delx = []; 
delx_mode = in.delx_mode;
rhs = -dW.*z_x; 
opts = in; % update_delx is now using opts as input parameter
update_mode = 'init1';
update_delx;


while ~done
    iter = iter+1;
    x_k = xk_1;
    
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;
    
    if sign(delx_vec(idelta)) == sign(pk_old(idelta)) && iter > 1
        delta = 0; flag = 0;
    else
        pk = pk_old;
        dk = A'*(A(:,gamma_xh)*delx);
        
        if flag == 0 && dW(idelta) < sign(pk(idelta))*dk(idelta)
            delta = 0; flag = 1;
        else
            temp_gamma = zeros(N,1);
            temp_gamma(gamma_xh) = gamma_xh;
            gamma_xc = find([1:N]' ~= temp_gamma);
            
            %% Calculate homotopy step size
            in.delta_flag = 2;
            in.pk = pk; in.dk = dk;
            in.ak = Wk; in.bk = dW;
            in.gamma = gamma_xh; in.gamma_c = gamma_xc;
            in.delx_vec = delx_vec; in.x = xk_1;
            out = compute_delta(in);
            delta = out.delta; idelta = out.idelta;
            flag = out.flag;
        end
    end
    %% Update support
    e0 = epsilon;
    epsilon = e0 + delta;
    
    if epsilon > 1
        delta_end = 1-e0;
        xk_1 = x_k + delta_end*delx_vec;
        break;
    end
    
    xk_1 = x_k + delta*delx_vec;
    gamma_x_old = gamma_xh;
    
    % disp(sprintf(['iter = %d, delta = %3.4g, idelta = %d, flag = %d'], iter, delta, idelta, flag));
    itr_history = [itr_history; idelta delta flag];
    
    % update support
    update_supp; 

    pk_old = pk+delta*dk;
    Wk_old = Wk;
    Wk = (1-epsilon)*W_old + epsilon*W_new;
    pk_old([gamma_xh; idelta]) = sign(pk_old([gamma_xh; idelta])).*Wk([gamma_xh; idelta]);
    
    % update delx
    z_x = -sign(pk_old);
    rhs = -dW.*z_x;
    update_mode = 'update';
    update_delx;
    %     AtAgx = (A(:,gamma_xh)'*A(:,gamma_xh));
    %     delx2 = AtAgx\rhs(gamma_xh); % -AtAgx\(dW(gamma_xh).*z_x);
    %     figure(112); plot([delx delx2]);
    %     if norm(delx-delx2) > 1e-5
    %         stp = 1;
    %     end
    
    constr_violation = nnz((abs(pk_old(gamma_xc))-Wk(gamma_xc))>1e-10);
    sign_violation = nnz((sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh)))>1);
    if constr_violation
        chk = gamma_xc((abs(pk_old(gamma_xc))-Wk(gamma_xc))>1e-10);
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
end

total_iter = iter;
total_time = cputime-t0;

out = [];
out.x_out = xk_1;
out.gamma = gamma_xh;
out.iter = total_iter;
out.time = total_time;
