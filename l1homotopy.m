% l1homotopy.m
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
% minimize_x  \|W x\|_1 + 1/2*\|Ax-y\|_2^2 + (1-epsilon)u'x,
%
%   u is defined as u = -W*sign(xh_old)-A'*(A*xh_old-y) 
%   and xh_old is an arbitrary warm-start vector 
%   (zero vector if no warm-start is available). 
%
%   The homotopy is solved by changing epsilon from 0 to 1. 
%   
%
% Optimality conditions:
%
%  A'(Ax-y)+u - epsilon(u)  = - W z    on Gamma
% |A'(Ax-y)+u - epsilon(u)| <   W      off Gamma
%
%   
% Inputs:
%  A -- M x N measurement matrix
%  y -- measurements
%
%  opts - input structure
%
%   xh_old  -- old signal estimate (if warm-start is not provided, then xh_old is set to zero)
%   gamma   -- support of xh_old
%
%   pk_old  -- A'(A*xh_old-y)+u
%   u       -- -W*sign(xh_old)-A'(A*xh_old-y)
%   W       -- weights for the L1 term... 
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
% out -- output structure
%   x_out -- output for BPDN
%   gamma -- support of the solution
%   iter -- number of homotopy iterations taken by the solver
%   time -- time taken by the solver
%   error_table -- error table with iteration record
%   iAtA on Gamma, or QR/Cholesky factors
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Web: http://users.ece.gatech.edu/~sasif/
%
% References: 
% 1. Sparse recovery of streaming signals using L1-homotopy, 
%       by M. Salman Asif and Justin Romberg 
% 2. Dynamic compressive sensing: Sparse recovery of streaming signals and video, 
%       by M. Salman Asif (Ph.D. thesis) 
%
%-------------------------------------------+
% Copyright (c) 2013.  M. Salman Asif
%-------------------------------------------+

% Change history
% 
% 05-29-13 -- name changed from l1homotopy_v2 to l1homotopy

function out = l1homotopy(A, y, opts)

N = size(A,2);
M = size(A,1);

% Use prior, related information if provided 
% Weights or regularization parameter... 
if isfield(opts,'W'); W = opts.W; else W = opts.tau; end
if isfield(opts,'xh_old'); 
    % old solution
    xh_old = opts.xh_old; 
    % old solution constraints
    pk_old = opts.pk_old;
    % old solution support
    gamma_xh = opts.gamma;
    % dummy variable (..)
    u = opts.u;
else
    xh_old = zeros(N,1); 
end

% no initial solution
if norm(xh_old) == 0;     
    pk_old = -A'*y;
    [tau_old gamma_xh] = max(abs(pk_old));
    z_x = zeros(N,1);
    z_x(gamma_xh) = -sign(pk_old(gamma_xh));
    u = -W.*z_x-pk_old;
    pk_old = pk_old+u;
end

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
if nnz(abs(A'*y) < W) == N
    out = opts;
    out.x_out = zeros(N,1);
    out.gamma = []; % find(abs(xk_1)>0);
    out.iter = 0;
    out.time = 0;
    disp('output is a zero vector');
    return;
end

% Make vectors out of scalar regularization parameters..
W = ones(N,1).*W;

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

temp_gamma = zeros(N,1);
temp_gamma(gamma_xh) = gamma_xh;
gamma_xc = find([1:N]' ~= temp_gamma);

z_x = zeros(N,1);
z_x(gamma_xh) = -sign(pk_old(gamma_xh));
pk_old(gamma_xh) = sign(pk_old(gamma_xh)).*W(gamma_xh);
pk = pk_old;
dk = 0*pk;
xk_1 = xh_old;

% Initial step setup
idelta = gamma_xh(1); flag = 1;    

% initialize delx
in_delx = [];
delx_mode = opts.delx_mode;
rhs = u;

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
        dk = A'*dk_temp-u;
        
        %%%--- compute step size
        in = [];
        
        % Setting shrinkage_flag to zero shrinks new active constraint towards the
        % final value instantly if doing so doesn't disturb the active set        
        in.shrinkage_flag = 2; % determines how to select the stepsize/support
        in.nonneg = nonneg; % imposes non-negativity constraint on the solution
        in.pk = pk; in.dk = dk;
        in.ak = W; 
        in.gamma = gamma_xh; in.gamma_c = gamma_xc;
        in.delx_vec = delx_vec; in.x = xk_1; 
        out = compute_delta_v2(in);
        delta = out.delta; idelta = out.idelta;
        flag = out.flag;
    end
    e0 = epsilon;
    epsilon = e0 + delta;
    
    if epsilon > 1
        delta_end = 1-e0;
        xk_1 = x_k + delta_end*delx_vec;
        pk_old = pk + delta_end*dk;
        pk_old(gamma_xh) = sign(pk_old(gamma_xh)).*W(gamma_xh);
        break;
    end
    
    xk_1 = x_k + delta*delx_vec;
    gamma_old = gamma_xh;
    
    itr_history = [itr_history; idelta delta flag];
    
    % update support
    update_supp;

    pk_old = pk+delta*dk;  
    pk_old([gamma_xh; idelta]) = sign(pk_old([gamma_xh; idelta])).*W([gamma_xh; idelta]);
    
    % update delx
    z_x = -sign(pk_old);
    rhs = u;
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
    constr_violation = nnz((abs(pk_old(gamma_xc))-W(gamma_xc))>1e-10);
    sign_violation = nnz((sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh)))>1);
    if constr_violation && nonneg == 0
        chk = gamma_xc((abs(pk_old(gamma_xc))-W(gamma_xc))>1e-10);
        stp = 1;
        fprintf('problem... with constraint violation -- %s \n', mfilename);
        fprintf('Refactorize the matrix... recompute delx (consider using qr delx_update instead of mil) \n');
        % some times it comes here due to bad conditioning of AtAgx.
        update_mode = 'init0';
        update_delx;
    end
    if sign_violation>=1 && nonneg == 0
        chk = gamma_xh(sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh))>1);
        stp = 1;
        fprintf('problem... sign mismatch -- %s \n',mfilename);
        fprintf('Refactorize the matrix... recompute delx (consider using qr delx_update instead of mil)\n');
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
        plot(1:N, W,'--k','MarkerSize',12);
        plot(1:N, -W, '--k','MarkerSize',12); 
        
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

