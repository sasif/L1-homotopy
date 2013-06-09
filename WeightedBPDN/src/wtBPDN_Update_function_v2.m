% BPDN_homotopy_ReWeighted.m
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \|x\|_1 + 1/2*||y-AWx||_2^2
%
% and dynamically update the weights w_i
%
% Inputs:
% A - m x n measurement matrix
% y - measurement vector
% in - input structure
%   tau - value of regularization parameter
%   x_old - old estimate
%   gamma - old signal support
%   AtA - A(:,gamma_xh)'*A(:,gamma_xh);
%   iAtA  - inv(AtA)
%   or QR/Cholesky factors
%   pk_old - old constrtaints set
%   yhat - dummy vector for old measuremnents
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
% Email: sasif@ece.gatech.edu
% Created: March 16, 2011

function out = wtBPDN_Update_function_v2(A, y, in);

% x_old, gamma_xh, AtAgx, iAtAgx, pk_old, W_old, W_new, maxiter
x_old = in.x_old;
gamma_xh = in.gamma;
pk_old = in.pk_old;
yhat = in.yhat; % dummy y (for warm-start)
% ds = in.ds; % dummy d (for warm-start)
tau = in.tau;
maxiter = in.maxiter;

t0 = cputime;
N = size(A,2);
M = size(A,1);

% pk_old = A'*(A*x_k-y);
% pk_old(gamma_xh) = sign(pk_old(gamma_xh)).*W_old(gamma_xh);
xk_1 = x_old;
ds = pk_old;
d = A'*(yhat-y);

epsilon = 0;

done = 0;
iter = 0;
itr_history = [];
idelta = gamma_xh(1); flag = 1;

% initialize delx
in_delx = [];
delx_mode = in.delx_mode;
rhs = ds-d;
update_mode = 'init1';
update_delx;

while ~done
    iter = iter+1;
    x_k = xk_1;
    z_x = -sign(pk_old);
    
    % AtAgx_inv = inv(A(:,gamma_xh)'*A(:,gamma_xh));
    % AtAgx = (A(:,gamma_xh)'*A(:,gamma_xh));
    % norm(iAtAgx-AtAgx_inv)
    % delx = AtAgx\(ds(gamma_xh)-d(gamma_xh));
    % delx computation moved below...
    
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;
    
    if sign(delx_vec(idelta)) == sign(pk_old(idelta)) && iter > 1
        delta = 0; flag = 0;
    else
        pk = pk_old;
        dk = A'*(A(:,gamma_xh)*delx)+d-ds;
        
        temp_gamma = zeros(N,1);
        temp_gamma(gamma_xh) = gamma_xh;
        gamma_xc = find([1:N]' ~= temp_gamma);
        
        %% Calculate homotopy step size
        delta1_constr = (tau-pk)./(dk);
        delta1_constr = delta1_constr(gamma_xc);
        delta2_constr = (-tau-pk)./(dk);
        delta2_constr = delta2_constr(gamma_xc);
        delta3_constr = (-x_k(gamma_xh)./delx);
        idelta_1 = (delta1_constr>0);
        idelta_2 = (delta2_constr>0);
        idelta_3 = (delta3_constr>0);
        delta1 = min(delta1_constr(idelta_1));
        delta2 = min(delta2_constr(idelta_2));
        delta3 = min(delta3_constr(idelta_3));
        if isempty(delta1); delta1 = inf; end
        if isempty(delta2); delta2 = inf; end
        if isempty(delta3); delta3 = inf; end
        if delta1>delta2
            delta = delta2;
            idelta = gamma_xc(delta2_constr==delta2);
            flag = 1;
        else
            delta = delta1;
            idelta = gamma_xc(delta1_constr==delta1);
            flag = 1;
        end
        if delta3 <= delta
            delta = delta3;
            idelta = gamma_xh(delta3_constr == delta3);
            flag = 0;
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
    update_supp
    
    pk_old = pk+delta*dk;
    pk_old([gamma_xh; idelta]) = sign(pk_old([gamma_xh; idelta])).*tau;
    
    % update delx
    rhs = ds-d;
    update_mode = 'update';
    update_delx;
    
    %     figure(1); plot([Wk(gamma_xh)-abs(A(:,gamma_xh)'*(A*xk_1-y))])
    %     figure(2); plot([Wk(gamma_xh) abs(A(:,gamma_xh)'*(A*xk_1-y))])
    %     figure(3); plot(delx_vec);
    
    % disp(sprintf(['idelta = %d, delta = %3.4g, flag = %d'],idelta, delta, flag));
    
    constr_violation = nnz((abs(pk_old(gamma_xc))-tau)>1e-10);
    sign_violation = nnz((sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh)))>1);
    if constr_violation
        chk = gamma_xc((abs(pk_old(gamma_xc))-tau)>1e-10);
        stp = 1;
        fprintf('problem... with constraint violation -- %s', mfilename);
        % some times it comes here due to bad conditioning of AtAgx.
        fprintf('Refactorize the matrix... recompute delx \n');
        update_mode = 'init0';
        update_delx;
    end
    if sign_violation>=1
        chk = gamma_xh(sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh))>1);
        stp = 1;
        fprintf('problem... sign mismatch -- %s',mfilename);
        fprintf('Refactorize the matrix... recompute delx \n');
        update_mode = 'init0';
        update_delx;
    end
    % figure(1);
    % plot([abs(A(:,gamma_xh)'*(A*x_k-(1-e0)*yhat-e0*y)+(1-e0)*ds(gamma_xh))-tau abs(A(:,gamma_xh)'*(A*xk_1-(1-epsilon)*yhat-epsilon*y)+(1-epsilon)*ds(gamma_xh))-tau])
    % plot(abs(A(:,gamma_xh)'*(A*xk_1-(1-epsilon)*yhat-epsilon*y)+(1-epsilon)*ds(gamma_xh))-tau)
end

total_iter = iter;
total_time = cputime-t0;

out = [];
out.x_out = xk_1;
out.gamma = gamma_xh;
out.iter = total_iter;
out.time = total_time;
