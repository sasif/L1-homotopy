% wtBPDN_homotopy_function
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum w_i |x_i| + 1/2*||y-Ax||_2^2
%
% using homotopy from scratch. 
%
% (A simpler way to solve weighted BPDN is to modulate columns of A using
% the weights and solve simple BPDN, i.e., min \|x_i\|_1 + 1/2\|AWx-y\|_2^2
% Inputs:
% A - m x n measurement matrix
% y - measurement vector
% in - input structure
%   W - final values of regularization parameter
%   maxiter - maximum number of homotopy iterations
%   Te -
%   record - record iteration history
%   x_orig - origianl signal for error history
%
% Outputs:
% out - output structure
%   x_out - output for BPDN
%   gamma - support of the solution
%   iter - number of homotopy iterations taken by the solver
%   time - time taken by the solver
%   error_table - error table with iteration record
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
%
%-------------------------------------------+
% Copyright (c) 2011.  Muhammad Salman Asif
%-------------------------------------------+

function out = wtBPDN_homotopy_function(A, y, in)

N = size(A,2);
M = size(A,1);

W_vec = in.W;
maxiter = in.maxiter;
Te = inf;
if isfield(in,'Te')
    Te = in.Te;
end
err_record = 0;
if isfield(in,'record');
    err_record = in.record;
    if err_record
        x_orig = in.x_orig;
    end
end
t0 = cputime;

% Regularization parameters
unique_eps = sort(unique(W_vec),'descend');

% Initialization of primal sign and support
z_x = zeros(N,1);
gamma_x = [];       % Primal support

% Initial step
pk_old = -A'*y;
constr_mask = abs(pk_old)>W_vec;
[c idelta] = max(abs(pk_old.*constr_mask));
eps_iter = sum(unique_eps>c)+1;

gamma_xh = idelta;
temp_gamma = zeros(N,1);
temp_gamma(gamma_xh) = gamma_xh;
gamma_xc = find([1:N]' ~= temp_gamma);

z_x(gamma_xh) = -sign(pk_old(gamma_xh));
epsilon = c;
pk_old(gamma_xh) = sign(pk_old(gamma_xh))*epsilon;
xk_1 = zeros(N,1);

% loop parameters
done = 0;
iter = 0;
itr_history = [];

error_table = [];
if err_record
    error_table = [epsilon norm(xk_1-x_orig) 1];
end

% initialize delx
in_delx = [];
delx_mode = in.delx_mode;
indicator_temp = epsilon>W_vec;
rhs = indicator_temp.*z_x;
update_mode = 'init0';
update_delx;
flag = 1;

while iter < maxiter
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')
    
    x_k = xk_1;
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%% update on x %%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Update direction
    % update direction
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;
    
    if sign(delx_vec(idelta)) == sign(pk_old(idelta)) && iter > 1 && flag == 0
        delta = 0; flag = 0;
    else
        pk = pk_old;
        % dk = AtA*del_x_vec;
        dk_temp = A*delx_vec;
        dk = A'*dk_temp;
        
        %%%--- compute step size
        in = [];
        
        % Setting shrinkage_flag to zero shrinks new active constraint towards the
        % final value instantly if doing so doesn't disturb the active set
        
        epsilon_temp = epsilon.*(epsilon>W_vec)+W_vec.*(epsilon<=W_vec);
        one_temp = epsilon>W_vec;
        
        in.delta_flag = 2;
        in.pk = pk; in.dk = dk;
        in.ak = epsilon_temp; in.bk = -one_temp;
        in.gamma = gamma_xh; in.gamma_c = gamma_xc;
        in.delx_vec = delx_vec; in.x = xk_1;
        out = compute_delta(in);
        delta = out.delta; idelta = out.idelta;
        flag = out.flag;
        
        xk_1 = x_k+delta*delx_vec;
        pk_old = pk+delta*dk;
        epsilon_old = epsilon;
        epsilon = epsilon-delta;
        
        if epsilon <= unique_eps(eps_iter)
            epsilon = unique_eps(eps_iter);
            delta_end = epsilon_old-epsilon;
            pk_old = pk+delta_end*dk;
            epsilon_temp = epsilon.*(epsilon>W_vec)+W_vec.*(epsilon<=W_vec);
            pk_old([gamma_xh]) = sign(pk_old([gamma_xh])).*epsilon_temp([gamma_xh]);
            
            xk_1 = x_k + delta_end*delx_vec;
            eps_iter = eps_iter+1;
            if eps_iter > length(unique_eps)
                %disp('done!');
                break;
            else
                %disp('switch epsilon!');
                flag = 1;
                z_x = -sign(pk_old);
                indicator_temp = epsilon>W_vec;
                rhs = indicator_temp.*z_x;
                update_mode = 'recompute';
                update_delx;
                continue;
            end
        end
        if epsilon <= min(W_vec); %sqrt(2*log(N))*sigma; %1e-7 %|| iter > 5*T || (length(gamma_lambda) == K)
            delta_end = epsilon_old-thresh;
            pk_old = pk+delta_end*dk;
            xk_1 = x_k + delta_end*delx_vec;
            % disp('done!');
            break;
        end
        if length(gamma_x)-flag >= Te
            total_time = cputime-t0;
            break;
        end
    end
    % disp(sprintf(['iter = %d, delta = %3.4g, idelta = %d, flag = %d'], iter, delta, idelta, flag));
    itr_history = [itr_history; idelta delta flag];
    
    % update support
    update_supp;
    
    epsilon_temp = epsilon.*(epsilon>W_vec)+W_vec.*(epsilon<=W_vec);
    pk_old([gamma_xh; idelta]) = sign(pk_old([gamma_xh; idelta])).*epsilon_temp([gamma_xh; idelta]);
    
    % update delx
    z_x = -sign(pk_old);
    indicator_temp = epsilon>W_vec;
    rhs = indicator_temp.*z_x;
    update_mode = 'update';
    update_delx;
        
    constr_violation = nnz((abs(pk_old(gamma_xc))-epsilon_temp(gamma_xc))>1e-10);
    sign_violation = nnz(abs(sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh)))>1);
    if constr_violation
        chk = gamma_xc((abs(pk_old(gamma_xc))-epsilon_temp(gamma_xc))>1e-10);
        stp = 1;
        fprintf('problem... with constraint violation -- %s\n', mfilename);
        fprintf('Refactorize the matrix... recompute delx \n');
        % some times it comes here due to bad conditioning of AtAgx.
        update_mode = 'init0';
        update_delx;
    end
    if sign_violation>1
        chk = gamma_xh(abs(sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh)))>1);
        stp = 1;
        fprintf('problem... sign mismatch -- %s\n',mfilename);
        fprintf('Refactorize the matrix... recompute delx \n');
        update_mode = 'init0';
        update_delx;
    end
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
out.gamma = gamma_xh;
out.iter = total_iter;
out.time = total_time;
out.error_table = error_table;