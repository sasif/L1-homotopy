% Simulation of homotopy update for signals with large dynamic range.
% Comparison with FPC_AS and simple BPDN homotopy
% Need FPC package and accompanying problems to test parts of this file.

clear
% close all
% clear classes;

% % load fixed random states
% load RandomStates
% rand('state', s_rand);
% randn('state', s_randn);

sim_runs = 100;
% tau_table = [1e-1 1e-3 1e-5];
outer_iterations = 1;

stack_time_solve = zeros(sim_runs*outer_iterations,2);
stack_time_update = zeros(sim_runs*outer_iterations,2);
stack_iter = zeros(sim_runs*outer_iterations,5);

% Test difficult problems, adopted from FPC_AS package.
FPC_problems = 2;

if FPC_problems == 1
    % Following list of problems can be used from FPC_AS package
    Problist = {'Ameth6Xmeth2n1024m512k154seed200', 'Ameth6Xmeth6n1024m512k154seed200', ...
        'CaltechTest1', 'CaltechTest2', 'CaltechTest3', 'CaltechTest4'};
    di = 6; % choose the problem type
    load( Problist{di}, 'b','Omega','n','xs');

    m = length(b);
    x = xs;
    q = (find(abs(xs)>0));
    T = length(q);
    large_mag = max(abs(x(q)));
    small_mag = min(abs(x(q)));
    DCT_mat = dct(eye(n));
    A = DCT_mat(Omega,:);
else
    n = 512;
    m = 228;
    T = 30;
    q = randperm(n);
    T1 = 15;
    T2 = T-T1;
    x = zeros(n,1);
    large_mag = 1e6;
    small_mag = 1e2;
    x(q(1:T1)) = randsrc(T1,1)*large_mag;
    x(q(T1+1:T1+T2)) = randsrc(T2,1)*small_mag;

    % measurement matrix
    A = randn(m,n)/sqrt(m);
    % orthonormalize rows
    % A = orth(A')';

end

% noisy observations
sigma = 0.0; % how much noise to add should depend on power of Ax
e = randn(m,1)*sigma;
y = A*x + e;

% regularization parameter
tau = .01*small_mag/large_mag*max(abs(A'*y)); 
%tau = tau_table; 
% if sigma>0
%     tau = sigma * sqrt(log(n)*2); % BPDN
%     %tau = max(abs(A'*(A*x-y))); % ideal ???
% end

% Homotopy method
tic
[xp, gamma_x, iter_xp, t_homotopy1] = BPDN_homotopy_function(A, y, tau,4*n); %BPDN
t_homotopy1 = toc;
tolA_h = tau*sum(abs(xp))+1/2*(norm(A*xp-y))^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The relaxation in the tolerance of objective value for FPC
% The performance of FPC is best for the small values of these parameters 
% and that increases the iteration count.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdg_scale = 1e-3; 
opts.gtol = 1e-8;

tolA = tolA_h*(1+pdg_scale);

M = [];
opts.f_value_tol = tolA;

tic
[x_FPC, Out_FPC] = FPC_AS(n,A,y,tau,M,opts); % call FPC_AS
t_FPC = toc
[t_homotopy1 t_FPC]

for out_iter = 1:outer_iterations
    for sim_iter = 1:sim_runs

        % Model the change in signal
        dx_model = 1;
        switch dx_model
            case  1
                % Add some elements (T_in) or remove some existing elements (T_out)
                dx  = zeros(n,1);
                T_in = randint(1,1,round(T/10))+1;
                T_out = randint(1,1,round(T/10));
                new_large_small = round(rand(T_in,1));
                qin = randperm(n);
                dx(qin(1:T_in)) = randn(T_in,1).*(new_large_small.*large_mag + (1-new_large_small).*small_mag);
                
                qout = randperm(T);
                oldx_indices = q(1:T);
                dx(oldx_indices(qout(1:T_out))) = -x(oldx_indices(qout(1:T_out))).*rand(T_out,1);
                %dx(oldx_indices(qout(T_out+1:end))) = randn(T-T_out,1)*.01;
            case  2;
                % Perturb the current locations by Gaussian noise
                oldx_indices = q(1:T);
                dx = zeros(n,1);
                dx(oldx_indices) = .1*randn(T,1);
                T_in = 0;
                T_out = 0;
            otherwise
                disp('DIY');
        end
        xt = x+dx;
        yt = y+A*dx;

        pk_old = A'*(A*xp-y);
        pk_old(gamma_x) = sign(pk_old(gamma_x))*tau;
        AtAgx = A(:,gamma_x)'*A(:,gamma_x);
        iAtAgx = inv(AtAgx);

        % Update using homotopy
        tic;
        [xp_h, gamma_xh, iter_xp_update, t_homotopy_update] = DynamicX_BPDN_function(A, AtAgx, iAtAgx, y, yt, xp, gamma_x, pk_old, tau,2*n);
        t_homotopy_update = toc;
        
        tolA_h = tau*sum(abs(xp_h))+1/2*(norm(A*xp_h-yt))^2;

        % Homotopy method from scratch
        tic
        [xp2, gamma_x2, iter_xp2, t_homotopy2] = BPDN_homotopy_function(A, yt, tau,2*n); %BPDN
        t_homotopy2 = toc;
        tolA_h = tau*sum(abs(xp))+1/2*(norm(A*xp-y))^2;
        
        % % l1_ls
        % [x_l1_ls_update,status,history] = l1_ls(A,yt,2*tau,1e-4);
        % t_l1_ls_update = history(7,end);
        % tolA_ls = history(2,end)/2;

        tolA = tolA_h*(1+pdg_scale);
        
        opts.x0 = x_FPC;
        opts.f_value_tol = tolA;
        M = [];

        tic;
        [x_FPC_update, Out_FPC_update] = FPC_AS(n,A,yt,tau,M,opts); % call FPC_AS
        t_FPC_update = toc;
        
        stack_time_solve((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy2 t_FPC];
        stack_time_update((out_iter-1)*sim_runs+sim_iter,:) = [t_homotopy_update t_FPC_update];
        stack_iter((out_iter-1)*sim_runs+sim_iter,:) = [iter_xp2 iter_xp_update (Out_FPC_update.nProdA+Out_FPC_update.nProdAt)/2 norm(xt-xp_h) norm(xt-x_FPC_update)];
        [out_iter sim_iter]
    end
end

% Generate table for performance
disp(['       DynamicX          BPDN homotopy          FPC_AS']);
disp(['   nProdAtA    CPU     nProdAtA     CPU    nProdAtA     CPU'])
table_stats = zeros(length(outer_iterations),6);
table_stats(:,3) = mean(reshape(stack_iter(:,1),sim_runs,outer_iterations));
table_stats(:,1) = mean(reshape(stack_iter(:,2),sim_runs,outer_iterations));
table_stats(:,5) = mean(reshape(stack_iter(:,3),sim_runs,outer_iterations));
table_stats(:,4) = mean(reshape(stack_time_solve(:,1),sim_runs,outer_iterations));
table_stats(:,2) = mean(reshape(stack_time_update(:,1),sim_runs,outer_iterations));
table_stats(:,6) = mean(reshape(stack_time_update(:,2),sim_runs,outer_iterations));
table_stats
