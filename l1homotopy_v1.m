% l1homotopy_v1.m
%
% A general program that solves homotopy for
%   a weighted LASSO/BPDN problem with or without a warm-start vector
%   
%   Some examples for dynamic updating include
%   sequential measurements
%   time-varying signal
%   iterative reweighting
%   measurement replacement
%   dictionary learning
%   Kalman-type filtering
%   Streaming signal recovery 
%
%   or any other problem that can be written in the following form:
%
% minimize_x  \Sum [1-epsilon]w_old + epsilon w_new]|x_i| + 1/2*||Ax-(1-epsilon)yh-epsilon y||_2^2 + (1-epsilon)u'x,
%
% where solution is updating by changing epsilon from 0 toward 1.
%   yh = A*xh_old, 
%   u = -w.*sign(xh_old) (where w can be the old or the new weights)
%
% Optimality conditions:
%
% A'(Ax-yh)+u + epsilon(A'(yh-y)-u) = -(w_old + epsilon(w_new-w_old))z    on Gamma
% |A'(Ax-yh)+u + epsilon(A'(yh-y)-u)| < w_old + epsilon(w_new-w_old)    off Gamma
%
%   
% Inputs:
%  A -- M x N measurement matrix
%  y -- measurements
%
%  opts - input structure
%
%   xh_old  -- old signal estimate
%   gamma   -- support of xh_old
%
%   (LHS terms)
%   pk_old  -- A'(A*xh_old-yh)+u
%   Atdyu   -- A'(yh-y)-u
%
%   (RHS terms)
%   W_old  -- w_old
%   W_new  -- w_new
%
%   AtAgx and iAtAgx (i.e., A(:,gamma)'*A(:,gamma) and its inverse)
%   or QR/Cholesky factors
%
%   delx_mode -- mode for rank-1 update ('mil', 'chol', or 'qr')
%   nonneg  -- add nonneg constraint on the solution? (default = 0)
%   maxiter -- maximum number of homotopy iterations
%   Te      -- maximum support size allowed
%   record  -- record iteration history
%   x_orig  -- origianl signal for error history
%   debias  -- debias the solution at the end
%   early_terminate -- terminate early if the support is identified
%                   (useful only in high SNR settings)
%   verbose -- print output after every verbose steps
%   plots   -- plot the solution at every iteration after verbose steps
%
% Outputs:
% out - output structure
%   x_out - output for BPDN
%   gamma - support of the solution
%   iter - number of homotopy iterations taken by the solver
%   time - time taken by the solver
%   error_table - error table with iteration record
%   iAtA on Gamma, or QR/Cholesky factors
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
%
%-------------------------------------------+
% Copyright (c) 2012.  M. Salman Asif
%-------------------------------------------+

function out = l1homotopy_v1(A, y, opts)

N = size(A,2);
M = size(A,1);

% Use prior, related information if provided 
% old solution
if isfield(opts,'xh_old'); xh_old = opts.xh_old; else xh_old = zeros(N,1); end
% old solution constraints 
if isfield(opts,'pk_old'); pk_old = opts.pk_old; else pk_old = -A'*y; end
% old solution support 
if isfield(opts,'gamma'); gamma_xh = opts.gamma; else gamma_xh = []; end
if isempty(gamma_xh); [tau_old gamma_xh] = max(abs(pk_old)); end
% dummy variable (..)
if isfield(opts,'Atdyu'); Atdyu = opts.Atdyu; else Atdyu = 0; end
% new regularization parameters
if isfield(opts,'W_new'); W_new = opts.W_new; else W_new = opts.tau; end
% old regularization parameters 
if isfield(opts,'W_old'); W_old = opts.W_old; else W_old = tau_old; end
% no initial solution
if norm(xh_old) == 0; pk_old = -A'*y; Atdyu = 0; [tau_old gamma_xh] = max(abs(pk_old)); W_old = tau_old; end
% input is a zero vector
if norm(y) == 0; 
    out = opts;
    out.x_out = zeros(N,1);
    out.gamma = []; % find(abs(xk_1)>0);
    out.iter = 0;
    out.time = 0;
    disp('input is a zero vector');
    return;
end
% output is a zero vector
if nnz(abs(A'*y) < W_new) == N
    out = opts;
    out.x_out = zeros(N,1);
    out.gamma = []; % find(abs(xk_1)>0);
    out.iter = 0;
    out.time = 0;
    disp('output is a zero vector');
    return;
end

% Make vectors out of scalar regularization parameters..
W_old = ones(N,1).*W_old; W_new = ones(N,1).*W_new;

% maximum iterations 
if isfield(opts,'maxiter'); maxiter = opts.maxiter; else maxiter = 2*N; end
% non-negativity constraint
if isfield(opts,'nonneg'); nonneg = opts.nonneg; else nonneg = 0; end
% maximum support size
if isfield(opts,'Te'); Te = opts.Te; else Te = inf; end
% record error/history
if isfield(opts,'record'); err_record = opts.record; else err_record = 0; end
if err_record; err_fun = opts.err_fun; end 
% debiasing step at the end (solve LS on the support)
if isfield(opts,'debias'); debias = opts.debias; else debias = 0; end
% early terminate if residual of restricted LS falls below certain
% threshold
if isfield(opts,'early_terminate'); early_terminate = opts.early_terminate; else early_terminate = 0; end
% print output
if isfield(opts,'verbose'); verbose = opts.verbose; else verbose = 0; end
% debug plots
if isfield(opts,'plots'); plots = opts.plots; else plots = 0; end

%% GO
t0 = cputime;

%% Initial step
epsilon = 0;
Wk = (1-epsilon)*W_old + epsilon*W_new;
dW = W_new-W_old;

temp_gamma = zeros(N,1);
temp_gamma(gamma_xh) = gamma_xh;
gamma_xc = find([1:N]' ~= temp_gamma);

z_x = zeros(N,1);
z_x(gamma_xh) = -sign(pk_old(gamma_xh));
pk_old(gamma_xh) = sign(pk_old(gamma_xh)).*W_old(gamma_xh);
pk = pk_old;
dk = 0*pk;
xk_1 = xh_old;

% Initial step setup
idelta = gamma_xh(1); flag = 1;    

% initialize delx
in_delx = [];
delx_mode = opts.delx_mode;
rhs = -dW.*z_x-Atdyu;

if norm(xh_old)==0 && length(gamma_xh) == 1    
    update_mode = 'init0';
else    
    update_mode = 'init1';
end
update_delx;

%% loop parameters
done = 0;
iter = 0;
itr_history = [];
error_table = [];
if err_record
    error_table = [epsilon err_fun(xk_1) length(gamma_xh)];
end

while iter < maxiter
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')
    
    %% Homotopy
    x_k = xk_1;
    
    % Update direction
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;
    
    if ~isempty(idelta) && (sign(delx_vec(idelta)) == sign(pk_old(idelta)) && abs(x_k(idelta)) == 0)
        delta = 0; flag = 0;
    else
        pk = pk_old;
        % dk = AtA*delx_vec;
        % dk_temp = A*delx_vec;   
        dk_temp = mvprod(A,delx_vec,gamma_xh,0);
        dk = A'*dk_temp+Atdyu;
        
        %%%--- compute step size
        in = [];
        
        % Setting shrinkage_flag to zero shrinks new active constraint towards the
        % final value instantly if doing so doesn't disturb the active set        
        in.shrinkage_flag = 2; % determines how to select the stepsize/support
        in.nonneg = nonneg; % imposes non-negativity constraint on the solution
        in.pk = pk; in.dk = dk;
        in.ak = Wk; in.bk = dW;
        in.gamma = gamma_xh; in.gamma_c = gamma_xc;
        in.delx_vec = delx_vec; in.x = xk_1; 
        out = compute_delta(in);
        delta = out.delta; idelta = out.idelta;
        flag = out.flag;
    end
    e0 = epsilon;
    epsilon = e0 + delta;
    
    if epsilon > 1
        delta_end = 1-e0;
        xk_1 = x_k + delta_end*delx_vec;
        pk_old = pk + delta_end*dk;
        Wk = W_new;
        pk_old(gamma_xh) = sign(pk_old(gamma_xh)).*Wk(gamma_xh);
        break;
    end
    
    xk_1 = x_k + delta*delx_vec;
    gamma_old = gamma_xh;
    
    itr_history = [itr_history; idelta delta flag];
    
    % update support
    update_supp;

    pk_old = pk+delta*dk;
    Wk_old = Wk;
    Wk = (1-epsilon)*W_old + epsilon*W_new;
    pk_old([gamma_xh; idelta]) = sign(pk_old([gamma_xh; idelta])).*Wk([gamma_xh; idelta]);
    
    % update delx
    z_x = -sign(pk_old);
    rhs = -dW.*z_x-Atdyu;
    update_mode = 'update';
    update_delx;
    
    %     AtAgx = (A(:,gamma_xh)'*A(:,gamma_xh));
    %     delx2 = AtAgx\rhs(gamma_xh); % -AtAgx\(dW(gamma_xh).*z_x);
    %     figure(112); plot([delx delx2]);
    %     if norm(delx-delx2) > 1e-5
    %         stp = 1;
    %     end
    
    
    % Check convergence criterion (this can be useful)...
    if early_terminate
        if length(gamma_xh) < M/2
            xhat = zeros(N,1);
            % xhat(gamma_xh) = AtAgx\(A(:,gamma_xh)'*y);
            switch delx_mode
                case 'mil'
                    xhat(gamma_xh) = iAtA*(A(:,gamma_xh)'*y);
                case {'qr','chol'}
                    xhat(gamma_xh) = R\(R'\(A(:,gamma_xh)'*y));
            end
            if norm(y-A*xhat) < tau
                xk_1 = xhat;
                break;
            end
        end
    end
    
    
    if err_record        
        error_table = [error_table; epsilon err_fun(xk_1) length(gamma_xh)];
    end
    
    %% debug
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
    
    %% print and plot
    if mod(iter-1,verbose) == 0 && verbose
        fprintf(['iter = %d, delta = %3.4g, idelta = %d, flag = %d.\n'], iter, delta, idelta, flag);
    end
    if mod(iter-1,plots) == 0 && plots
        fig1 = figure(1234);
        subplot(2,1,1)
        hold off
        plot(pk,'.r', 'MarkerSize',14);
        hold on;
        plot(pk_old, 'LineWidth',1);
        
        if flag == 1
            plot(idelta, pk_old(idelta),'or','MarkerSize',18,'LineWidth',2);
            text(idelta, pk_old(idelta)*1.1, ['Incoming \gamma = ',num2str(idelta)],'FontSize',14);
        else
            plot(idelta, pk_old(idelta),'ok','MarkerSize',18,'LineWidth',2);
            text(idelta, pk_old(idelta)*1.1, ['Outgoing \gamma = ',num2str(idelta)],'FontSize',14);
        end
        set(gca,'FontSize',16, 'XLim',[1 N] );
        title(sprintf('BPDN shrinkage constraints: N = %d, M = %d', N, M));
        plot(1:N, Wk,'--k','MarkerSize',12);
        plot(1:N, -Wk, '--k','MarkerSize',12);
        plot(1:N, Wk_old,'--m','MarkerSize',12);
        plot(1:N, -Wk_old, '--m','MarkerSize',12);
        
        figure(fig1);
        subplot(2,1,2)
        hold off
        plot(x_k,'.r','MarkerSize',14); hold on;
        plot(xk_1,'LineWidth',1);
        if flag == 0
            plot(idelta, 0,'ok', 'MarkerSize',18,'LineWidth',2);
        else
            plot(idelta, 0,'or', 'MarkerSize',18,'LineWidth',2);
        end
        set(gca,'FontSize',16,'XLim',[1 N]);
        title(['Solution estimate at \epsilon = ',num2str(epsilon), ', iter. = ', num2str(iter)]);
        
        if iter == 1 && verbose
            disp('  ');
            disp('Every frame in the figure corresponds to a critical point on the homotopy path.')
            disp('Circle represents an incoming element, star represents an outgoing element.');
            disp(' ');
            disp('Put pause somewhere in the code to see this. ');
            disp('For now press some key to continue...');
            pause
        end
    end
    
end

%% debiasing step? 
if debias
    x_out = zeros(N,1);
    switch delx_mode
        case 'mil'
            x_out(gamma_xh) = iAtA*(A(:,gamma_xh)'*y);
        case {'qr','chol'}
            x_out(gamma_xh) = R\(R'\(A(:,gamma_xh)'*y));
    end
else
    x_out = xk_1;
end

%
if err_record
    error_table = [error_table; epsilon err_fun(x_out) length(gamma_xh)];
end
total_iter = iter;
total_time = cputime-t0;

%% Output the results
out = opts;
out.x_out = x_out;
out.gamma = gamma_xh; % find(abs(xk_1)>0);
out.iter = total_iter;
out.time = total_time;
out.error_table = error_table;
out.pk = pk_old;
switch delx_mode
    case 'mil'
        out.iAtA = iAtA;
    case 'qr'
        out.Q = Q;
        out.R = R;
    case 'chol'
        out.R = R;
    case 'qrM'
        out.Q0 = Q0;
        out.R0 = R0;
end

