% Simulation for dynamic update with sequential measurements.
% Compare warm start of GPSR and FPC with the homotopy update of BPDN.
% Need GPSR and FPC packages in your matlab path to run this simulation. 

clear
close all
clear classes;

sim_runs = 20; % Use a large number of simulation run to compare the performance
tau_table = [.5 .1 .05 .01];
outer_iterations = length(tau_table);
stack_time_solve = zeros(sim_runs*outer_iterations,4);
stack_time_update = zeros(sim_runs*outer_iterations,3);
stack_iter = zeros(sim_runs*outer_iterations,7);

for out_iter = 1:outer_iterations
    for sim_iter = 1:sim_runs
        % n is the original signal length
        n = 1024;
        % m is number of observations to make
        m = 512;

        % number of spikes to put down
        % n_spikes = floor(.01*n);
        ratio_sp = 5;
        n_spikes = round(m/ratio_sp);
        T = n_spikes;

        % random +/- 1 signal
        x = zeros(n,1);
        q = randperm(n);
        x(q(1:n_spikes)) = sign(randn(n_spikes,1));
        % x(q(1:n_spikes)) = randn(n_spikes,1);

        % measurement matrix
        G = randn(m,n)/sqrt(m);
        % orthonormalize rows
        % G = orth(G')';
        m_u =1;
        A = G(1:end-m_u,:);
        b = G(end-m_u+1:end,:);

        hR = @(x) A*x;
        hRt = @(x) A'*x;

        % noisy observations
        sigma = 0.01; % Add noise depending on the amplitude of y
        e = sigma*randn(m-m_u,1);
        y = A*x + e;

        % regularization parameter
        tau = tau_table(out_iter)*max(abs(A'*y));

        % Homotopy method
        tic
        [xp, gamma_x, iter_xp, t_homotopy1] = BPDN_homotopy_function(A, y, tau,4*n); %Lasso
        t_homotopy1 = toc;
        tolA_h = tau*sum(abs(xp))+1/2*(norm(A*xp-y))^2;

        % GPSR
        pdg_scale = 1e-12;
        tolA = tolA_h*(1+pdg_scale);
        stopCri=4;
        debias = 0;

        first_tau_factor = 0.8*(max(abs(A'*y))/tau);
        steps = 5;
        tic;
        x_BB_mono = xp;
        [x_BB_mono,x_debias_BB_mono,obj_BB_mono,...
            times_BB_mono,debias_start_BB_mono,mse]= ...
            GPSR_BB(y,hR,tau,...
            'Debias',debias,...
            'Monotone',1,...
            'AT',hRt,...
            'Initialization',0,...
            'Continuation',1,...
            'ContinuationSteps',steps,...
            'FirstTauFactor',first_tau_factor,...
            'StopCriterion',stopCri,...
            'ToleranceA',tolA,...
            'Verbose', 0, ...
            'ToleranceD',0.001);
        t_BB_mono = toc;
        
        % FPC
        opts.gtol = 1e-6;
        %opts.f_value_tol = tolA;
        opts.PrintOptions = 0;
        M = [];
        tic;
        [x_FPC, Out_FPC] = FPC_AS(n,A,y,tau,M,opts); % call FPC_AS
        t_FPC = toc;
        
        % New measurement
        chk_e0 = 1;  % 0 or 1: This selects if we want to take U = (A'A+b'b) (1) or (A'A+e0*b'b) (0)...!
        % In simulation U = (A'A+b'b) works better

        % Gram matrix update
        AtAgx = A(:,gamma_x)'*A(:,gamma_x)+chk_e0*b(:,gamma_x)'*b(:,gamma_x);
        iAtAgx = inv(AtAgx);
        pk_old = A'*(A*xp-y);
        w = b*x+randn(m_u,1)*sigma;

        % Update using homotopy
        tic
        [xp_h, gamma_xh, iter_xp_update, t_homotopy_update] =  DynamicSeq_BPDN_function(A, b, AtAgx, iAtAgx, y, w, xp, gamma_x, pk_old, tau, chk_e0, 4*n);
        t_homotopy_update = toc;
        
        
        
        % Compute using homotopy from scratch
        y2 = [y; w];
        G = [A; b];
        tic
        [xp2, gamma_x2, iter_xp2, t_homotopy2] = BPDN_homotopy_function(G, y2, tau,4*n); %Lasso
        t_homotopy2 = toc;
        
        tolA_h = tau*sum(abs(xp2))+1/2*(norm(G*xp2-y2))^2;

        % GPSR
        tolA = tolA_h*(1+pdg_scale);
        stopCri=4;
        debias = 0;

        % PUTTIN CONTINUATION IN HERE MAKES THINGS WORSE
        tic;
        [x_BB_mono_update,x_debias_BB_mono,obj_BB_mono,...
            times_BB_mono,debias_start_BB_mono,mse]= ...
            GPSR_BB(y2,G,tau,...
            'Debias',debias,...
            'Monotone',1,...
            'Initialization',x_BB_mono,...
            'StopCriterion',stopCri,...
            'Verbose', 0, ...
            'ToleranceA',tolA,...
            'ToleranceD',0.001);
        t_BB_mono_update = toc;

        % FPC
        opts.x0 = x_FPC;
%       opts.f_value_tol = tolA;
        M = [];
        tic
        [x_FPC_update, Out_FPC_update] = FPC_AS(n,G,y2,tau,M,opts); % call FPC_AS
        t_FPC_update = toc;
%         stack_time_solve((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy1 t_BB_mono];
%         stack_time_update((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy_update t_BB_mono_update];
%         stack_iter((out_iter-1)*sim_runs+sim_iter,:) = [iter_xp iter_xp_update norm(x-xp) norm(x-xp_h) norm(x-x_BB_mono_update) length(times_BB_mono)];

        stack_time_solve((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy2 t_homotopy1 t_BB_mono t_FPC];
        stack_time_update((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy_update t_BB_mono_update t_FPC_update];
        stack_iter((out_iter-1)*sim_runs+sim_iter,:) = [iter_xp2 iter_xp_update , length(times_BB_mono), (Out_FPC_update.nProdA+Out_FPC_update.nProdAt)/2 norm(x-xp_h)/norm(x) norm(xp_h-x_BB_mono_update)/norm(xp_h) norm(xp_h-x_FPC_update)/norm(xp_h)];
        [out_iter sim_iter]
    end
end

% Generate table for performance
disp(['      DynamicSeq         BPDN homotopy           GPSR            FPC_AS']);

table_stats = zeros(4,8);
table_stats(:,3) = mean(reshape(stack_iter(:,1),sim_runs,outer_iterations));
table_stats(:,1) = mean(reshape(stack_iter(:,2),sim_runs,outer_iterations));
table_stats(:,5) = mean(reshape(stack_iter(:,3),sim_runs,outer_iterations));
table_stats(:,7) = mean(reshape(stack_iter(:,4),sim_runs,outer_iterations));
table_stats(:,4) = mean(reshape(stack_time_solve(:,1),sim_runs,outer_iterations));
table_stats(:,2) = mean(reshape(stack_time_update(:,1),sim_runs,outer_iterations));
table_stats(:,6) = mean(reshape(stack_time_update(:,2),sim_runs,outer_iterations));
table_stats(:,8) = mean(reshape(stack_time_update(:,3),sim_runs,outer_iterations));
table_stats