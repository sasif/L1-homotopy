% Simulation to test the performance of REC

clear
close all

% % load fixed random states
% load RandomStates
% rand('state', s_rand);
% randn('state', s_randn);

sim_runs = 10;
additional_elements = 100;

% No. of new measurements added at a time
m_u = 10;
inner_iterations = additional_elements/m_u;

stack_time = zeros(sim_runs*inner_iterations,2);
stack_iter = zeros(sim_runs*inner_iterations,5);


for outer_iter = 1:sim_runs
    % % source length
    N = 150;
    % codeword length
    M = 2*N;
    M_st = M;
    
    % number of perturbations
    max_sparse_errors = .2; % e.g., 20% of M
    T = round(max_sparse_errors*M_st);
    
    T = round(N/2.5);
    
    % coding matrix
    Orth_mat = randn(M,M);
    A = randn(M,N)/sqrt(N);
    A = orth(A);

    % Annihilating projection matrix
    AtA = A'*A;
    iAtA = inv(AtA);
    AiAtA = A*iAtA;
    AiAtAAt = AiAtA*A';
    Q = eye(M)-AiAtAAt;

    % source word
    x = (randn(N,1));

    % channel: perturb T randomly chosen entries
    q = randperm(M);
    
    % Introduce sparse errors
    e = zeros(M,1);
    % e(q(1:T)) = randsrc(T,1); % Arbitrary sparse errors
    e(q(1:T)) = -A(q(1:T),:)*x; % Erasures!

    % Small noise
    x0 = randn(N,1);
    Ax0 = A*x0;
    sigma = median(abs(Ax0))/20; % control the power in small noise
    q_y = randn(M,1)*sigma;

    % Received codeword
    y = A*x+e+q_y;

    % Regularization parameter
    tau = 0.01;%*max(abs(Q*y)); % l1_ls
%     if sigma>0
%         tau = sigma * sqrt(log(N)*2); % BPDN
%         % tau = max(abs(Q'*q_y)); % ideal ???
%     end

    % Data recovery
    Qy = Q*y;
    [ep, gamma_e, ep_iter, t1] =  BPDN_homotopy_function(Q, Qy, tau, 4*M);
    % xp = inv(A'*A)*A'*(y-ep);
    xp = AiAtA'*(y-ep);

    tolA_h = tau*sum(abs(ep))+1/2*(norm(Q*(ep-y)))^2;

    for inn_iter = 1:inner_iterations
        [outer_iter inn_iter]
        if inn_iter > 1
%             % Update the parameters from the previous run of inner iterations
%             AtA = AtA_n;
%             iAtA = iAtA_n;
%             AiAtA = AiAtA_n;
%             AiAtAAt = AiAtAAt_n;
%             Q = Q_n;
%             y = y_n;
%             A = A_n;
%             q_y = q_y_n;
%             e = e_n;
%             xp = xp_n;
%             ep = ep_n;
%             M = M_n;
%             gamma_e = gamma_e_n;
%             e_BB_mono = e_BB_mono_n;
            % Update the parameters for next run of inner iterations
            AtA = FtF;
            iAtA = iFtF;
            AiAtA = FiFtF;
            AiAtAAt = FiFtFFt;
            Q = P;
            y = s;
            A = F;
            q_y = q_yw;
            e = c;
            xp = xp_h;
            ep = cp_h;
            M = M+m_u;
            gamma_e = gamma_h;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup for adding m_u new observations %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % No. of new observations
        % m_u =10;

        b = randn(m_u,N)/sqrt(N); % new rows in coding matrix
        max_sparse_errors_new = max_sparse_errors;
        max_sparse_errors_new = .2;
        T_new = round((m_u*max_sparse_errors_new)*(rand>=.5)+rand*.5); % new gross/sparse errors
        d = zeros(m_u,1);
        q_new = randperm(m_u);
        
        % Introduce gross errors in the new elements of the codeowrd
        % d(q_new(1:T_new)) = randsrc(T_new,1); % Arbitrary sparse errors
        d(q_new(1:T_new)) = -b(q_new(1:T_new),:)*x; % Erasures!

        % small noise in the new measurements
        q_w = randn(m_u,1)*sigma;

        w = b*x+d+q_w;
        F = [A; b];
        s = [y;w];
        c = [e; d];
        q_yw = [q_y; q_w];


        iAtAbt = iAtA*b';
        biAtAbt = b*iAtAbt;
        S_biAtAbt = inv(eye(m_u)+biAtAbt);
        AiAtAbt = AiAtA*b';
        iAtAbt_S = (iAtAbt*S_biAtAbt);
        AiAtAbt_S = A*iAtAbt_S;

        FtF = AtA+b'*b; %F'*F;
        iFtF = iAtA - iAtAbt_S*iAtAbt'; % inv(FtF);
        FiFtF = [AiAtA - AiAtAbt_S*iAtAbt'; iAtAbt'-biAtAbt*S_biAtAbt*iAtAbt']; % F*iFtF;
        FiFtFFt = [AiAtAAt - AiAtAbt_S*AiAtAbt' AiAtAbt-AiAtAbt_S*biAtAbt; AiAtAbt'-biAtAbt*AiAtAbt_S' biAtAbt-biAtAbt*S_biAtAbt*biAtAbt]; %FiFtF*F';
        P = eye(M+m_u)-FiFtFFt;

        dp = w-b*xp;
        z_d = sign(dp);
        cp_h = [ep; dp];
        gamma_n = M+find(abs(dp)>2*eps);
        gamma_n_old = gamma_n;
        gamma_h = [gamma_e; gamma_n];
        epsilon = 0;
        e0 = 0;

        QgQ = Q(gamma_e,gamma_e);
        PgP = P(gamma_h,gamma_h);
        uQ = AiAtAbt_S(gamma_e,:);
        vQ = AiAtAbt(gamma_e,:);
        QgQ_update = QgQ + uQ*vQ';
        PgP_update = [[QgQ_update; P(gamma_n,gamma_e)] P(gamma_h,gamma_n)];
        PgP = PgP_update;

        iQgQ = inv(QgQ);
        iQgQ_update = iQgQ - (iQgQ*uQ)*(inv(eye(m_u)+vQ'*iQgQ*uQ))*(vQ'*iQgQ);

        P11 = QgQ_update;
        P12 = P(gamma_e,gamma_n);
        P21 = P12';
        P22 = P(gamma_n,gamma_n);
        S_P = inv(P22-P21*iQgQ_update*P12);
        iP11_P12 = iQgQ_update*P12;
        iPgP_update = [iQgQ_update+(iP11_P12*S_P)*iP11_P12' -iP11_P12*S_P; -S_P*iP11_P12' S_P];
        iPgP = iPgP_update;

        pk_old = P*(cp_h-s);
        [cp_h gamma_h cp_h_iter th] = DynamicSeq_REC_function(P, s, iPgP, cp_h, gamma_h, gamma_n, pk_old, tau, M, m_u, 4*M);
        xp_h = FiFtF'*(s-cp_h);

        % % Check the solution using homotopy from scratch.
        [cp2, gamma_c2, cp2_iter, t2] = BPDN_homotopy_function(P, P*s, tau, 4*M);
        
        stack_time((outer_iter-1)*inner_iterations+inn_iter,:) = [t2 th];
        stack_iter((outer_iter-1)*inner_iterations+inn_iter,:) = [cp2_iter cp_h_iter T_new norm(x-xp_h)/norm(x) norm(cp_h-cp2)];
%         iter_count = (outer_iter-1)*inner_iterations+inn_iter
    end
end

['BPDN      REC']
Average_time = mean(stack_time,1)
Average_iterations = mean(stack_iter(:,1:2))
