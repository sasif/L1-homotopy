% Simulation for BPDN homotopy update of DynamicX.
% Compare warm start of GPSR and FPC with the homotopy update of BPDN.
% Need GPSR and FPC packages in your matlab path to run this simulation. 

clear 
clear classes; close all; 

% % load fixed random states
% load RandomStates
% rand('state', s_rand);
% randn('state', s_randn);

sim_runs = 10; % Run for large number of simulations to compare the performance.
tau_table = [.5, .1, .05, .01];
outer_iterations = length(tau_table);

stack_time_solve = zeros(sim_runs*outer_iterations,4);
stack_time_update = zeros(sim_runs*outer_iterations,3);
stack_iter = zeros(sim_runs*outer_iterations,9);

for out_iter = 1:outer_iterations
    for sim_iter = 1:sim_runs
        % n is the original signal length
        n = 1024;

        % k is number of observations to make
        k = 512;

        % number of spikes to put down
        % n_spikes = floor(.01*n);
        ratio_sp = 5;
        n_spikes = round(k/ratio_sp);
        T = n_spikes;

        % random +/- 1 signal
        x = zeros(n,1);
        q = randperm(n);
        x(q(1:n_spikes)) = sign(randn(n_spikes,1));
        % x(q(1:n_spikes)) = randn(n_spikes,1);

        % measurement matrix
        A = randn(k,n)/sqrt(k);
        % orthonormalize rows
        % A = orth(A')';

        % noisy observations
        sigma = 0.01; % how much noise to add should depend on power of Ax
        e = randn(k,1)*sigma;
        y = A*x + e;

        % regularization parameter
        tau = tau_table(out_iter)*max(abs(A'*y)); % l1_ls
        % if sigma>0
        %     tau = sigma * sqrt(log(n)*2); % BPDN
        %     %tau = max(abs(A'*(A*x-y))); % ideal ???
        % end

        % Homotopy method
        tic
        [xp, gamma_x, iter_xp, t_homotopy1] = BPDN_homotopy_function(A, y, tau,4*n); %BPDN
        t_homotopy1 = toc;
        tolA_h = tau*sum(abs(xp))+1/2*(norm(A*xp-y))^2;

        % Adjust tolerance for the objective value, relative to the
        % homotopy solution (which we assume to be accurate!).
        pdg_scale = 1e-6;

        tolA = tolA_h+pdg_scale;
        stopCri=4;
        debias = 0;
        x_BB_mono = xp;
        t_BB_mono = 0;
        first_tau_factor = 0.8*(max(abs(A'*y))/tau);
        steps = 5;
        
        tic;
        [x_BB_mono,x_debias_BB_mono,obj_BB_mono,...
            times_BB_mono,debias_start_BB_mono,mse]= ...
            GPSR_BB(y,A,tau,...
            'Debias',debias,...
            'Monotone',1,...
            'Initialization',0,...
            'Continuation',1,...
            'ContinuationSteps',steps,...
            'FirstTauFactor',first_tau_factor,...
            'StopCriterion',stopCri,...
            'Verbose', 0,...
            'True_x',xp,...
            'MaxiterA', 1000,...
            'ToleranceA',tolA,...
            'ToleranceD',0.001);
        t_BB_mono = toc;

        % FPC
        opts.mxitr = 4*n;
        opts.f_value_tol = tolA;
        M = [];
        tic
        [x_FPC, Out_FPC] = FPC_AS(n,A,y,tau,M,opts); % call FPC_AS
        t_FPC = toc;
        
        % Model the change in signal
        dx_model = 1;
        switch dx_model
            case  1
                % Add some elements (T_in) or remove some existing elements (T_out)
                dx  = zeros(n,1);
                T_in = randint(1,1,round(T/20));
                T_out = randint(1,1,round(T/20));
                qin = randperm(n);
                dx(qin(1:T_in)) = (randn(T_in,1));
                qout = randperm(T);
                oldx_indices = q(1:T);
                dx(oldx_indices(qout(1:T_out))) = -x(oldx_indices(qout(1:T_out))).*rand(T_out,1);
                dx(oldx_indices(qout(T_out+1:end))) = randn(T-T_out,1)*.1;
            case  2;
                % Perturb the current locations by Gaussian noise
                oldx_indices = q(1:T);
                dx = zeros(n,1);
                dx(oldx_indices) = .1*randn(T,1);
                T_in = 0;
                T_out = 0;
            otherwise
                disp('Nooooo');
        end
        xt = x+dx;
        yt = y+A*dx;

        pk_old = A'*(A*xp-y);
        AtAgx = A(:,gamma_x)'*A(:,gamma_x);
        iAtAgx = inv(AtAgx);

        % Update using homotopy
        tic
        [xp_h, gamma_xh, iter_xp_update, t_homotopy_update2] = DynamicX_BPDN_function(A, AtAgx, iAtAgx, y, yt, xp, gamma_x, pk_old, tau,4*n);
        t_homotopy_update = toc;
        tolA_h = tau*sum(abs(xp_h))+1/2*(norm(A*xp_h-yt))^2;

        tic
        [xp2, gamma_x2, iter_xp2, t_homotopy2] = BPDN_homotopy_function(A, y, tau,4*n); %BPDN
        t_homotopy2 = toc;

        % % l1_ls
        % [x_l1_ls_update,status,history] = l1_ls(A,yt,2*tau,1e-4);
        % t_l1_ls_update = history(7,end);
        % tolA_ls = history(2,end)/2;

        % GPSR
        tolA = tolA_h+pdg_scale;
        stopCri=4;
        debias = 0;
        % x_BB_mono = xp;
        % t_BB_mono = 0;
        % PUTTIN CONTINUATION IN HERE MAKES THINGS WORSE
        tic;
        [x_BB_mono_update,x_debias_BB_mono,obj_BB_mono,...
            times_BB_mono,debias_start_BB_mono,mse]= ...
            GPSR_BB(yt,A,tau,...
            'Debias',debias,...
            'Monotone',1,...
            'Initialization',x_BB_mono,...
            'StopCriterion',stopCri,...
            'ToleranceA',tolA,...
            'True_x',xp_h,...
            'MaxiterA', 1000,...
            'Verbose', 0,...
            'ToleranceD',0.001);
        t_BB_mono_update = toc;
        
        % FPC
        opts.x0 = x_FPC;
        opts.f_value_tol = tolA;
        M = [];
        tic
        [x_FPC_update, Out_FPC_update] = FPC_AS(n,A,yt,tau,M,opts); % call FPC_AS
        t_FPC_update = toc;
        
        stack_time_solve((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy2 t_homotopy1 t_BB_mono t_FPC];
        stack_time_update((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy_update t_BB_mono_update t_FPC_update];
        % New mat files rangTau2 last entry is iteration count for GPSR
        stack_iter((out_iter-1)*sim_runs+sim_iter,:) = [iter_xp2 iter_xp_update , length(times_BB_mono), (Out_FPC_update.nProdA+Out_FPC_update.nProdAt)/2 norm(xt-xp_h)/norm(xt) norm(xp_h-x_BB_mono_update)/norm(xp_h) norm(xp_h-x_FPC_update)/norm(xp_h) norm(xt-x_FPC_update)/norm(xt) Out_FPC_update.itr];
        [out_iter sim_iter]        
    end
end

% Generate table for performance
disp(['       DynamicX          BPDN homotopy           GPSR            FPC_AS']);
table_stats = zeros(outer_iterations,8);
table_stats(:,3) = mean(reshape(stack_iter(:,1),sim_runs,outer_iterations));
table_stats(:,1) = mean(reshape(stack_iter(:,2),sim_runs,outer_iterations));
table_stats(:,5) = mean(reshape(stack_iter(:,3),sim_runs,outer_iterations));
table_stats(:,7) = mean(reshape(stack_iter(:,4),sim_runs,outer_iterations));
table_stats(:,4) = mean(reshape(stack_time_solve(:,1),sim_runs,outer_iterations));
table_stats(:,2) = mean(reshape(stack_time_update(:,1),sim_runs,outer_iterations));
table_stats(:,6) = mean(reshape(stack_time_update(:,2),sim_runs,outer_iterations));
table_stats(:,8) = mean(reshape(stack_time_update(:,3),sim_runs,outer_iterations))
