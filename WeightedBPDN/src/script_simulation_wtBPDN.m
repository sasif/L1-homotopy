% partial script for running different solvers

%% Main loop
SIM_stack = cell(maxsim,1);
SIM_memory = cell(maxsim,1);

for sim = 1:maxsim
    %% Generate a random signal
    in = []; in.type = sType; in.T = T; in.randgen = 1; in.take_fwt = 0;
    switch sType
        case 'blocks'
            in.take_fwt = 1;
            in.wType = 'haar';
        case 'heavisine'
            in.take_fwt = 1;
            in.wType = 'daub4';
        case 'pcwpoly';
            in.take_fwt = 1;
            in.wType = 'daub8';
    end
    x = genSignal(N,in);
    % gamma_orig = find(abs(x)>0);
    [val ind] = sort(abs(x),'descend');
    ind_pos = ind(val>0);
    gamma_orig = ind_pos(1:min(length(ind_pos),M-1));
    
    % measurement matrix
    in = []; in.type = mType;
    A = genAmat(M,N,in);
    
    % measurements
    sigma = sqrt(norm(A*x)^2/10^(SNR/10)/M);
    %sigma = .05;
    e = randn(M,1)*sigma;
    y = A*x+e;
    
    % % orthogonalize rows of A even if A itself is not orthogonal
    %     [Q, R] = qr(A',0);
    %     A = Q'; y = R' \ y;
    
    % parameter selection
    if lambda > 0
        tau = lambda*max(abs(A'*y));
    else
        % tau = sigma*sqrt(log(N));
        tau = max(1e-4*max(abs(A'*y)),sigma*sqrt(log(N)));        
    end
    % tau = max(lambda*max(abs(A'*y)),sigma*sqrt(log(N)));
    
    maxiter = 2*N;
    
    %% Adaptive support & weight selection
    in = [];
    in.A = A; in.y = y; in.x = x;
    in.x_init = zeros(N,1); in.max_rwt = rwt_adp;
    in.tau = tau;
    in.rwt_mode = rwt_mode;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    if in.verbose
        fprintf('sType-%s, mType-%s, (N,M,T)-%d,%d,%d ',sType,mType,N,M,T);    
    end
    out = script_rwtBPDN_adaptive(in);
    xh_adp = out.x_out;
    iter_adp = out.iter;
    gamma_adp = out.gamma;
    W_new = out.W_new;
    xh_adp_init = out.x_init;
    err_adp = out.err;
    time_adp = out.time;
    supp_diff_adp = out.supp_diff;
    
    %% Iterative reweighted BPDN (IRW-H)
    in = [];
    in.A = A; in.y = y; in.x = x;
    in.x_init = zeros(N,1); in.max_rwt = rwt;
    in.tau = tau;
    in.rwt_mode = rwt_mode;
    in.delx_mode = delx_mode;
    out = script_rwtBPDN_iterative(in);
    xh_rwt = out.x_out;
    xh_rwt_init = out.x_init;
    iter_rwt = out.iter;
    gamma_rwt = out.gamma;
    err_rwt = out.err;
    time_rwt = out.time;
    supp_diff_rwt = out.supp_diff;
    
    %% Oracle rwt BPDN
    alpha = 5; beta = 10;
    W_new = ones(N,1); W_new(gamma_orig) = 1/alpha./(beta*abs(x(gamma_orig)));
    % epsilon = 0.1; W_new = 1/alpha./(beta*(abs(x))+epsilon);
    
    % To check the accuracy of rwtBPDN_adaptive set
    % W_new = out.W_new/tau;
    
    in = [];
    in.tau = tau;
    in.maxiter = maxiter;
    in.x_orig = x;
    in.record = 1;
    in.delx_mode = 'mil';
    %     in.Te = rwt_itr;
    AW = A*diag(1./W_new);
    tic
    out_new = BPDN_homotopy_function(AW, y, in); %BPDN
    % time_orac = out_new.time;
    time_orac = toc;
    xh_orac = out_new.x_out.*(1./W_new);
    gamma_orac = out_new.gamma;
    iter_orac = out_new.iter;
    
    %% oracle LS
    x_LS = zeros(N,1);
    x_LS(gamma_orig) = A(:,gamma_orig)\y;
    
    %% OMP
    in = [];
    in.Te = round(1.2*T);
    out = OMP_function(y,A,in);
    x_omp = out.x_out;
    iter_omp = out.iter;
    gamma_omp = out.gamma;
    
    %% FPC_AS
    %     opts = [];
    %     opts.gtol = 1e-6;
    %     opts.xtol = 1e-6;
    %     opts.fullMu = false;
    %     opts.x0 = zeros(N,1);
    %     opts.eta = 4;
    %     opts.mxitr = N;
    %     % opts.PrintOutput = 0;
    %     % opts.f_value_tol = tolA;
    %     % opts.PrintOptions = 0;
    %     % opts.record = -1;
    %     % opts.scale_A = 1;
    %     W_new = ones(N,1);
    %     AW = A;
    %     x_fpc = [];
    %     iter_fpc = [];
    %     err_fpc = [];
    %     supp_diff_fpc = [];
    %     for rwt_itr = 1:rwt+1
    %         % [x_FPC, out_FPC] = FPC_AS(N,AW,y,tau,[],opts); % call FPC_AS
    %         out_FPC = fpc(N,AW,y,1/tau,[],opts);
    %         out_FPC.nProdA = out_FPC.itr;
    %         out_FPC.nProdAt = out_FPC.itr;
    %         x_FPC = out_FPC.x;
    %
    %         x_fpc = x_FPC./W_new;
    %
    %         [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_fpc,M);
    %         W_new = 1/alpha./(beta*abs(x_fpc)+epsilon);
    %         opts.x0 = x_fpc.*W_new;
    %         AW = A*diag(1./W_new);
    %         iter_fpc = [iter_fpc (out_FPC.nProdA+out_FPC.nProdAt)/2];
    %         err_fpc = [err_fpc norm(x-x_fpc)/norm(x)];
    %         supp_diff_fpc = [supp_diff_fpc length(setxor(gamma_orig,find(x_fpc)))];
    %     end
    
    %% SPGL1
    % call SPGL1
    iter_spgl1 = []; err_spgl1 = []; supp_diff_spgl1 = []; time_spgl1 = [];
    if ~exist('spgSetParms','file'); error('Solver SPGL1 is not found.'); end
    
    W_new = ones(N,1);
    % delta = norm(e);
    delta = sqrt(M)*sigma;
    x_spgl1 = [];
    
    for rwt_itr = 1:rwt+1
        spg_opts = spgSetParms('verbosity',0,'weights',W_new,'iterations',maxiter);
        [x_spgl1,r,g,info] = spgl1(A,y,0,delta,x_spgl1,spg_opts);
        
        rerr = norm(x_spgl1-x)/norm(x);
        iter_spgl1 = [iter_spgl1  (info.nProdA+info.nProdAt)/2];
        err_spgl1 = [err_spgl1 rerr];
        supp_diff_spgl1 = [supp_diff_spgl1 length(setxor(gamma_orig,find(x_spgl1)))];
        time_spgl1 = [time_spgl1 info.timeTotal];
        [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_spgl1,M);
        W_new = 1/alpha./(beta*abs(x_spgl1)+epsilon);
    end
    
    %% YALL1
    % set options
    digit = 6; if sigma > 0; digit = 4; end
    opts = [];
    opts.tol = 10^(-digit);
    W_new = ones(N,1);
    opts.weights = W_new;
    opts.print = 0;
    opts.maxit = maxiter/2;
    x_yall1 = [];
    
    iter_yall1 = [];
    err_yall1 = [];
    time_yall1 = [];
    supp_diff_yall1 = [];
    for rwt_itr = 1:rwt+1
        opts.nu = 0; opts.rho = tau; 
        tic;
        [x_yall1,Out] = yall1(A,y,opts);
        % time_yall1 = [time_yall1 Out.cputime];
        time_yall1 = [time_yall1 toc];
        
        iter_yall1 = [iter_yall1 (Out.cntA+Out.cntAt)/2];
        rerr = norm(x-x_yall1)/norm(x);
        err_yall1 = [err_yall1 rerr];
        
        supp_diff_yall1 = [supp_diff_yall1 length(setxor(gamma_orig,find(x_yall1)))];
        
        [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_yall1,M);
        W_new = 1/alpha./(beta*abs(x_yall1)+epsilon);
        opts.x0 = x_yall1;
        opts.weights = W_new;
    end
    
    %% SpaRSA
    W_new = 1;
    % AW = A;
    x_sparsa = 0;
    iter_sparsa = [];
    err_sparsa = [];
    time_sparsa = [];
    supp_diff_sparsa = [];
    for rwt_itr = 1:rwt+1
        psi_function = @(x,tau) soft(x,tau*W_new);
        phi_function = @(x) sum(abs(W_new.*x));
        tic;
        [x_SpaRSA,x_debias_SpaRSA_m,obj_SpaRSA_m_cont,...
            times_SpaRSA_m_cont,debias_start_SpaRSA_m,mse_SpaRSA_m,taus_SpaRSA_m, numA, numAt]= ...
            SpaRSA_adpW(y,A,tau,...
            'Monotone',0,...
            'adp_wt',0,...
            'W_new',W_new,...
            'Debias',0,...
            'Initialization',x_sparsa,...
            'StopCriterion',2,...            
            'ToleranceA',1e-4,...
            'psi',psi_function,...
            'phi',phi_function,...
            'Safeguard',1,...
            'MaxiterA',maxiter,...
            'Verbose',0,...
            'True_x',x,...
            'Continuation',1,...
            'Continuationsteps',-1);
        x_sparsa = x_SpaRSA;
        
        time_sparsa = [time_sparsa toc];
        % time_sparsa = [time_sparsa times_SpaRSA_m_cont(end)];

        iter_sparsa = [iter_sparsa (numA+numAt)/2];
        rerr = norm(x-x_sparsa)/norm(x);
        err_sparsa = [err_sparsa rerr];
        
        
        if rerr > 100
            fprintf('rerr big... iter # %d in %s\n',sim,str0);
        end
        
        supp_diff_sparsa = [supp_diff_sparsa length(setxor(gamma_orig,find(x_sparsa)))];
        
        [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_sparsa,M);
        W_new = 1/alpha./(beta*abs(x_sparsa)+epsilon);
    end
    
    %% NESTA
    %     opts = [];
    %     opts.Verbose = 0;
    %     opts.tolvar = 1e-8;
    %     delta = 0;
    %     muf = 1e-8;
    %     W_new = ones(N,1);
    %
    %     iter_nesta = [];
    %     err_nesta = [];
    %     supp_diff_nesta = [];
    %     x_nesta = [];
    %
    %     for rwt_itr = 1:rwt+1
    %         opts.U = spdiags(W_new,0,N,N);
    %         opts.Ut = opts.U;
    %         opts.normU = max(W_new);
    %         opts.xplug = x_nesta;  % use old solution as starting value
    %
    %         % constrained
    %         [x_nesta,niter,resid,outData,optsOut] = NESTA(A,[],y,muf,delta,opts);
    %         % unconstrained
    %         La = norm( A*A' );
    %         % [x_nesta,niter,resid,outData,optsOut] = NESTA_UP(A,[],y,tau,La,muf,opts);
    %
    %         iter_nesta = [iter_nesta niter*3];
    %         rerr = norm(x-x_nesta)/norm(x);
    %         err_nesta = [err_nesta rerr];
    %         supp_diff_nesta = [supp_diff_nesta length(setxor(gamma_orig,find(x_nesta)))];
    %
    %         [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_nesta,M);
    %         W_new = 1/alpha./(beta*abs(x_nesta)+epsilon);
    %     end
    
    %% SALSA
    %     W_new = ones(N,1);
    %     x_salsa = 0;
    %     mu = 1/tau;
    %     invLS = @(x) inv(A'*A+mu*eye(N))*x;
    %
    %     iter_salsa = [];
    %     err_salsa = [];
    %     supp_diff_salsa = [];
    %     for rwt_itr = 1:rwt+1
    %         psi_function = @(x,tau) soft(x,tau*W_new);
    %         phi_function = @(x) sum(abs(W_new.*x));
    %         [x_salsa, numA, numAt, objective, distance,  times, mses] = ...
    %             SALSA_v2(y, A, tau,...
    %             'Mu', 1e-3, ...
    %             'True_x', x, ...
    %             'psi',psi_function,...
    %             'phi',phi_function,...
    %             'INITIALIZATION', x_salsa, ...
    %             'ToleranceA', 1e-6,...
    %             'MAXITERA', maxiter, ...
    %             'VERBOSE', 0);
    %
    %         iter_salsa = [iter_salsa (numA+numAt)/2];
    %         rerr = norm(x-x_salsa)/norm(x)
    %         err_salsa = [err_salsa rerr];
    %         supp_diff_salsa = [supp_diff_salsa length(setxor(gamma_orig,find(x_salsa)))];
    %
    %         [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_salsa,M);
    %         W_new = 1/alpha./(beta*abs(x_salsa)+epsilon);
    %     end
    
    %% Save results...
    
    % fprintf('sim iter. %d: \n',sim);
    % SIM_stack = [SIM_stack; sim, rwt_itr, iter_adp, iter_rwt_ALL, norm(x-xh_adp), norm(x-xh_rwt),length(setxor(gamma_adp,gamma_rwt)),length(setdiff(gamma_adp,gamma_orig)),length(setdiff(gamma_rwt,gamma_orig))];
    SIM_stack{sim} = [sim, tau, norm(x-x_LS)/norm(x), norm(x-xh_orac)/norm(x), iter_orac, ...
        norm(x-xh_adp)/norm(x), sum(iter_adp,2), sum(time_adp,2), ...
        norm(x-xh_rwt)/norm(x), sum(iter_rwt,2), sum(time_rwt,2), ...
        norm(x-x_yall1)/norm(x), sum(iter_yall1,2), sum(time_yall1,2), ...
        norm(x-x_sparsa)/norm(x), sum(iter_sparsa,2), sum(time_sparsa,2), ...
        norm(x-x_spgl1)/norm(x), sum(iter_spgl1,2), sum(time_spgl1,2), ...
        norm(x-x_omp)/norm(x), iter_omp];
    % norm(x-x_fpc)/norm(x), sum(iter_fpc,2), ...
    % sum(iter_nesta,2), norm(x-x_nesta)/norm(x), ...
    % sum(iter_salsa,2), norm(x-x_salsa)/norm(x), ...
    
    % adaptive rwt
    exp = 1;
    SIM_memory{sim}{exp,1} = 'adp';
    SIM_memory{sim}{exp,2} = iter_adp;
    SIM_memory{sim}{exp,3} = err_adp;
    SIM_memory{sim}{exp,4} = supp_diff_adp;
    SIM_memory{sim}{exp,5} = time_adp;
    % iterative rwt
    exp = exp+1;
    SIM_memory{sim}{exp,1} = 'rwt';
    SIM_memory{sim}{exp,2} = iter_rwt;
    SIM_memory{sim}{exp,3} = err_rwt;
    SIM_memory{sim}{exp,4} = supp_diff_rwt;
    SIM_memory{sim}{exp,5} = time_rwt;
    % yall1
    exp = exp+1;
    SIM_memory{sim}{exp,1} = 'yall1';
    SIM_memory{sim}{exp,2} = iter_yall1;
    SIM_memory{sim}{exp,3} = err_yall1;
    SIM_memory{sim}{exp,4} = supp_diff_yall1;
    SIM_memory{sim}{exp,5} = time_yall1;
    % sparsa
    exp = exp+1;
    SIM_memory{sim}{exp,1} = 'sparsa';
    SIM_memory{sim}{exp,2} = iter_sparsa;
    SIM_memory{sim}{exp,3} = err_sparsa;
    SIM_memory{sim}{exp,4} = supp_diff_sparsa;
    SIM_memory{sim}{exp,5} = time_sparsa;
    % spgl1
    exp = exp+1;
    SIM_memory{sim}{exp,1} = 'spgl1';
    SIM_memory{sim}{exp,2} = iter_spgl1;
    SIM_memory{sim}{exp,3} = err_spgl1;
    SIM_memory{sim}{exp,4} = supp_diff_spgl1;
    SIM_memory{sim}{exp,5} = time_spgl1;
    %     % fpc
    %     exp = exp+1;
    %     SIM_memory{sim}{exp,1} = 'fpc';
    %     SIM_memory{sim}{exp,2} = iter_fpc;
    %     SIM_memory{sim}{exp,3} = err_fpc;
    %     SIM_memory{sim}{exp,4} = supp_diff_fpc;
    %     % nesta
    %     exp = exp+1;
    %     SIM_memory{sim}{6,1} = 'nesta';
    %     SIM_memory{sim}{6,2} = iter_nesta;
    %     SIM_memory{sim}{6,3} = err_nesta;
    %     SIM_memory{sim}{6,4} = supp_diff_nesta;
    %     % salsa
    %     exp = exp+1;
    %     SIM_memory{sim}{6,1} = 'salsa';
    %     SIM_memory{sim}{6,2} = iter_salsa;
    %     SIM_memory{sim}{6,3} = err_salsa;
    %     SIM_memory{sim}{6,4} = supp_diff_salsa;
    
    % mS = mean(cell2mat(SIM_stack),1);
    % fprintf('maxsim %d. oracLS=%3.4g. (err,iter,time): oracWT-%3.4g,%3.4g; adp-%3.4g,%3.4g,%3.4g; rwt-%3.4g,%3.4g,%3.4g; yall1-%3.4g,%3.4g,%3.4g; sparsa-%3.4g,%3.4g,%3.4g; omp-%3.4g,%3.4g; fpc-%3.4g,%3.4g; spgl1-%3.4g,%3.4g. \n', maxsim, mS(2:end));
end

mS = mean(cell2mat(SIM_stack),1);
str1 = 'maxsim %d. tau = %3.4g, oracLS=%3.4g. (err,iter,time): oracWT-%3.4g,%3.4g; adp-%3.4g,%3.4g,%3.4g; rwt-%3.4g,%3.4g,%3.4g; yall1-%3.4g,%3.4g,%3.4g; sparsa-%3.4g,%3.4g,%3.4g; spgl1-%3.4g,%3.4g,%3.4g; omp-%3.4g,%3.4g.';
str2 = sprintf([str1,' \n'], maxsim, mS(2:end));
disp(str2);

