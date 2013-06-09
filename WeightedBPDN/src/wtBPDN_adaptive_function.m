% wtBPDN_adaptive_function
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \tau \e_i |x_i| + 1/2*||y-Ax||_2^2
% where at every homotopy step the e_i corresponding to the incoming element
% is reduced to a very small value (e.g., 1e-6). This way active and
% inactive indices have different values of the regularization parameter.
% The active elements have smaller weight, so they can stay active as long
% as there is no change in sign. The inactive elements have higher weight,
% so they are pushed towards zero. The hope is that, in this adaptive
% procedure, only the true elements will become active and stay active
% along the homotopy path.
%
% In homotopy, such an adaptive weight selection strategy can be included
% at every step without any additional cost.
%
% Inputs:
% A - m x n measurement matrix
% y - measurement vector
% in - input structure
%   tau - final value of regularization parameter
%
%   shrinkage_mode - (fixed or adaptive weight selection)
%       (for more details, see shrinkage_update.m)
%           'Tsteps':   active weights are reduced to tau/ewt
%           'frac':     active weights are divided by ewt
%           'Trwt':     active weights updated as w_i = tau/(beta*x_i)
%           'rwt':      ...
%           'OLS': set according to LS solution on the active support
%               such as w_gamma = 1./abs((A_gamma'*A_gamma)^-1*A_gamma'*y);
%               
%           "Tsteps" signifies the observation that using Tsteps along with
%           a large value of ewt often yields solution in T steps.
%
%   ewt - Selects weight factor for the active elements (Use a large value)
%           ewt controls tradeoff b/w speed and accuracy
%           higher ewt --> quicker but potentially unstable
%           (ewt=1, shrinkage_mode = Tsteps && shrinkage_flag=2) solves 
%           standard LASSO homotopy path
%
%   shrinkage_flag - 0, 1, or 2 (default is 0)
%           0 - Instantaneously set active constraints to the "desired"
%               value as long as it does not interfere with the active set.
%               And if any constraint is violated, take care of that by
%               resetting the running value of epsilon
%               FAST BUT POTENTIALLY UNSTABLE
%           1 - Step size that causes a constraint violation by a member
%               of the inactive set
%           2 - Gradually change the active constraints, step size takes
%               into consideration both sign of elements in the active set
%               and constraint violations by the members of inactive set
%               (as is the case in standard homotopy)
%           3 - FLASH variation???
%
%   maxiter - maximum number of homotopy iterations
%   Te - maximum support size allowed
%   omp - comparison with OMP
%   record - record iteration history
%   x_orig - origianl signal for error history
%   debias - debias the solution at the end
%   early_terminate - terminate early if the support is identified
%                   (useful only in high SNR settings)
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
% Email: sasif@gatech.edu
%
%-------------------------------------------+
% Copyright (c) 2012.  M. Salman Asif
%-------------------------------------------+

function out = wtBPDN_adaptive_function(A, y, in)

N = size(A,2);
M = size(A,1);

% Regularization parameters
tau = in.tau;
ewt = in.ewt;
shrinkage_mode = in.shrinkage_mode;
shrinkage_flag = 0;
if isfield(in,'shrinkage_flag')
    shrinkage_flag = in.shrinkage_flag;
end

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
omp = 0; % compare results with OMP
if isfield(in,'omp');
    omp = in.omp;
end
plots = 0; % debug plots
if isfield(in,'plots');
    plots = in.plots;
end
plot_wts = 0; % plot evolution of weights
if isfield(in,'plot_wts');
    plot_wts = in.plot_wts;
end
debias = 0;
if isfield(in,'debias')
    debias = in.debias;
end
early_terminate = 0;
if isfield(in,'early_terminate')
    early_terminate = in.early_terminate;
end

t0 = cputime;

%% Phase I (support selection)
% Initial step
z_x = zeros(N,1);
pk_old = -A'*y;
[c idelta] = max(abs(pk_old));

gamma_xh = idelta;
temp_gamma = zeros(N,1);
temp_gamma(gamma_xh) = gamma_xh;
gamma_xc = find([1:N]' ~= temp_gamma);

z_x(gamma_xh) = -sign(pk_old(gamma_xh));
epsilon = c;
pk_old(gamma_xh) = sign(pk_old(gamma_xh))*epsilon;
xk_1 = zeros(N,1);

%% loop parameters
done = 0;
iter = 0;
rwt_step2 = 0;

gamma_omp = gamma_xh;

error_table = [];
if err_record
    error_table = [epsilon norm(xk_1-x_orig) 1];
end


%% (selective support shrinkage)
epsilon_old = epsilon;
Supp_ledger = zeros(N,1);
Supp_ledger(idelta) = 1;

tau_vec = ones(N,1)*tau; % Final value of the regularization parameters
epsilon_vec = ones(N,1)*epsilon; % Running values of the regularization parameters

% Update target weights
shrinkage_update

% initialize delx
in_delx = [];
delx_mode = in.delx_mode;
rhs = (epsilon_vec-tau_vec).*z_x;
update_mode = 'init0';
update_delx;

while iter < maxiter
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')
    %% OMP comparison
    if omp && iter < M
        x_omp = zeros(N,1);
        x_omp(gamma_omp) = A(:,gamma_omp)\y;
        p_omp = A'*(y-A*x_omp);
        gamma_ompC = setdiff([1:N],gamma_omp);
        [val_omp, ind_omp] = max(abs(p_omp(gamma_ompC)));
        gamma_omp = [gamma_omp; gamma_ompC(ind_omp)];
    end
    
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
        dk_temp = A*delx_vec;
        dk = A'*dk_temp;
        
        %%%--- compute step size
        in = [];
        
        % Setting shrinkage_flag to zero shrinks new active constraint towards the
        % final value instantly if doing so doesn't disturb the active set
        in.shrinkage_flag = shrinkage_flag;
        in.pk = pk; in.dk = dk;
        in.ak = epsilon_vec; in.bk = tau_vec-epsilon_vec;
        in.gamma = gamma_xh; in.gamma_c = gamma_xc;
        in.delx_vec = delx_vec; in.x = xk_1;
        out = compute_delta(in);
        delta = out.delta; idelta = out.idelta;
        flag = out.flag;
        
        %% FLASH stepsize selection???        
        if shrinkage_flag == 3 && flag == 1
            delta_l = 0.5;
            delta_avg = out.delta_in*(1-delta_l)+delta_l; % Select stpe size as a convex combination of forward selection (FS) and LASSO...
            if delta_avg < out.delta_out
                delta = delta_avg;
                idelta = out.idelta_in;
            else
                delta = out.delta_out;
                idelta = out.idelta_out;
                flag = 0;
            end
        end
        
        if delta > 1
            delta = 1;
            flag = 1;
        end
        
        xk_1(gamma_xh) = x_k(gamma_xh)+delta*delx_vec(gamma_xh);
        pk_old = pk+delta*dk;
        
        epsilon_vec_old = epsilon_vec;
        % epsilon_vec(gamma_xh) = (1-delta)*epsilon_vec(gamma_xh)+delta*tau_vec(gamma_xh);
        epsilon_vec = (1-delta)*epsilon_vec+delta*tau_vec;
                
        pk_old(gamma_xh) = sign(pk_old(gamma_xh)).*epsilon_vec(gamma_xh);
        
        
%         fig(333); plot(abs([A'*(A(:,gamma_xh)*x_orig(gamma_xh)-y) A'*(A(:,gamma_xh)*xk_1(gamma_xh)-y) tau_vec epsilon_vec]));
%         pause;
        
        % Check convergence criterion (this can be useful)...
        if early_terminate
                if length(gamma_xh) < M/2
                    xhat = zeros(N,1);
                    % xhat(gamma_xh) = AtAgx\(A(:,gamma_xh)'*y);
                    switch delx_mode
                        case 'mil'
                            xhat(gamma_xh) = iAtA*(A(:,gamma_xh)'*y);
                        case 'qr'
                            xhat(gamma_xh) = R\(R'\(A(:,gamma_xh)'*y));
                    end
                    if norm(y-A*xhat) < tau
                        xk_1 = xhat;
                        break;
                    end
                end
        end
        
        if max(abs(pk_old)) <= tau
            % if you want to solve exactly according to tau, uncomment the
            % following lines:
            %
            % delta_end = epsilon_old-tau;
            % xk_1(gamma_xh) = x_k(gamma_xh)+delta_end*delx_vec(gamma_xh);
            % pk_old = pk+delta_end*dk;
            % epsilon_vec = epsilon_vec_old;
            % epsilon_vec(gamma_xh) = (1-delta_end)*epsilon_vec(gamma_xh)+delta_end*tau_vec(gamma_xh);
            
            % disp('epsilon reduce below threshold');
            % fig(303); plot([xk_1(gamma_xh)-(A(:,gamma_xh)'*A(:,gamma_xh))\(A(:,gamma_xh)'*y-epsilon_vec(gamma_xh).*z_x(gamma_xh))])
            if flag == 0
                outx_index = find(gamma_xh==idelta);
                gamma_xh = [gamma_xh(1:outx_index-1); gamma_xh(outx_index+1:end)];
                
                xk_1(idelta) = 0;
                epsilon_vec(idelta) = epsilon;
            end
            break;
        end
        if length(gamma_xh) >= Te
            total_time = cputime-t0;
            % disp('support size exceeds limit');
            % fig(303); plot([xk_1(gamma_xh) (A(:,gamma_xh)'*A(:,gamma_xh))\(A(:,gamma_xh)'*y-tau*z_x(gamma_xh)) x_orig(gamma_xh)])
            % setxor(gamma_xh,find(abs(x_orig)>0))
            break;
        end
        
        %% Search for new element
        % The one that violates the constraint
        epsilon_old = epsilon;
        [epsilon index] = max(abs(pk_old));
        if flag == 1 && delta == 1
            if nnz(index == gamma_xh)
                % iter = iter-1;
                shrinkage_update;
                z_x = -sign(pk_old);
                rhs = (epsilon_vec-tau_vec).*z_x;
                
                switch delx_mode
                    case 'mil'
                        delx = iAtA*rhs(gamma_xh);
                    case 'qr'
                        delx = R\(R'\rhs(gamma_xh));
                end
                continue;
                % because the index already exists in the active set
            end
            idelta = index;
            epsilon_vec(idelta) = epsilon;
        elseif flag == 1
            epsilon_vec(idelta) = abs(pk_old(idelta)); % (1-delta)*epsilon_vec(idelta)+delta*tau_vec(idelta);
        end
    end

    if err_record
        error_table = [error_table; epsilon norm(xk_1-x_orig) length(gamma_xh)];
    end
    
    
    % update support
    update_supp;
    
    temp_gamma = zeros(N,1);
    temp_gamma(gamma_xh) = gamma_xh;
    gamma_xc = find([1:N]' ~= temp_gamma);
    epsilon_vec(gamma_xc) = epsilon;
    % epsilon_vec(gamma_xc) = max(abs(A'*y));
    
    if flag == 0
        Supp_ledge(idelta) = 0;
    else
        Supp_ledger(gamma_xh) = Supp_ledger(gamma_xh)+1;
    end
    
    %% Shrinkage parameters selection
    shrinkage_update;
    
    % update delx
    z_x = -sign(pk_old);
    rhs = (epsilon_vec-tau_vec).*z_x;
    in_delx.max_rec = 1;
    update_mode = 'update';
    update_delx;
    %     AtAgx = A(:,gamma_xh)'*A(:,gamma_xh);
    %     delx2 =  AtAgx\rhs(gamma_xh);% AtAgx\((epsilon_vec(gamma_xh)-tau_vec(gamma_xh)).*z_x(gamma_xh));
    %     fig(111); plot([delx delx2]);
    %     if norm(delx-delx2) > 1e-5
    %         stp = 1;
    %     end
    
    %% debug...
    if plots
        fig(101); plot([A'*(A*xk_1-y) pk_old epsilon_vec -epsilon_vec]);
        fig(102); clf; hold on;
        subplot(311); plot([abs(pk_old) epsilon_vec epsilon_vec_old]);
        title(sprintf('iter %d',iter));
        subplot(312); plot(delx_vec);
        subplot(313); hold on; stem(x_orig,'Marker','.'); plot([x_k xk_1]);
        pause(1/60);
        
        % fig(303); plot([xk_1(gamma_xh)-(A(:,gamma_xh)'*A(:,gamma_xh))\(A(:,gamma_xh)'*y-epsilon_vec(gamma_xh).*z_x(gamma_xh))])
        % [pk(gamma_xh) x_k(gamma_xh) pk_old(gamma_xh) xk_1(gamma_xh) x_orig(gamma_xh) gamma_xh]
        if (max(abs(A'*(A*xk_1-y))-abs(pk_old)) > 1e-8)
            disp('constraints mismatch...')
        end
    end
    constr_violation = nnz((abs(pk_old(gamma_xc))-epsilon_vec(gamma_xc))>1e-10);
    sign_violation = nnz(abs(sign(pk_old(gamma_xh))+sign(xk_1(gamma_xh)))>1);
    if constr_violation
        chk = gamma_xc((abs(pk_old(gamma_xc))-epsilon_vec(gamma_xc))>1e-10);
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
    
    %% Figure to view evolution of weights
    if plot_wts
        if mod(iter,10)==1
            if ~exist('marker_iter','var')
                weight_marker = {'b','r','k','m','g'};
                fig(121); clf;
                set(gca,'FontSize',16);
                semilogy(1:N,epsilon_vec([gamma_xh; gamma_xc]),'Color',weight_marker{1},'LineWidth',2);
                marker_iter = 1;
                hold on;
                EPS_VEC = [];
            else
                marker_iter = mod(marker_iter,length(weight_marker))+1;
                fig(121); semilogy(1:N,epsilon_vec([gamma_xh; gamma_xc]),'Color',weight_marker{marker_iter},'LineWidth',2);
            end
            axis tight;
            YLim = get(gca,'YLim');
            YLim1 = YLim;
            YLim1(1) = YLim(1)-(YLim(2)-YLim(1))*0.1;
            YLim1(2) = YLim(2)+(YLim(2)-YLim(1))*0.3;
            set(gca,'YLim',YLim1);
            str1 = sprintf('step %d',iter);
            % text(iter+5,epsilon*1.2, str1,'FontSize',16);
            EPS_VEC = [EPS_VEC epsilon_vec([gamma_xh; gamma_xc])];
        end
    end
end

if debias
    x_out = zeros(N,1);
    switch delx_mode
        case 'mil'            
            x_out(gamma_xh) = iAtA*(A(:,gamma_xh)'*y);
        case 'qr'
            x_out(gamma_xh) = R\(R'\(A(:,gamma_xh)'*y));
    end
else
    x_out = xk_1;
end

if err_record
    error_table = [error_table; epsilon norm(x_out-x_orig) length(gamma_xh)];
end
total_iter = iter;
total_time = cputime-t0;

out = [];
out.x_out = x_out;
out.gamma = gamma_xh; % find(abs(xk_1)>0);
out.iter = total_iter;
out.time = total_time;
out.error_table = error_table;
out.tau_vec = epsilon_vec;

% if ~isempty(setxor(gamma_xh,find(abs(x_orig)>0))) || iter > nnz(x_orig)
%     stp = 1;
% end
% fig(303); plot([xk_1(gamma_xh)-(A(:,gamma_xh)'*A(:,gamma_xh))\(A(:,gamma_xh)'*y-tau_vec(gamma_xh).*z_x(gamma_xh))])
% fig(303); plot([xk_1(gamma_xh) (A(:,gamma_xh)'*A(:,gamma_xh))\(A(:,gamma_xh)'*y-tau_vec(gamma_xh).*z_x(gamma_xh)) x_orig(gamma_xh)])
