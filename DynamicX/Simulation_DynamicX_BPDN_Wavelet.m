% BPDN update with dynamic change in x
% Author: Salman Asif
% Created: November 2008
% Modified: February 2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the submitted paper
% 03/02/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DynamicX_..._Blocks_final_new
% DynamicX_..._PcwPoly_final_new
% DynamicX_..._Houste_final_new_Haar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need WaveLab environment in your MATLAB path to run this simulation.

clear
close all
% Model of the signal
model_sig = 1; % Test only 1 and 2

sim_runs = 10; % Use large number of simulations for evaluations
if model_sig == 2
    sim_runs = 255; % sim_runs should be the size of image for slices of image
end

stack_time=zeros(sim_runs,6);
stack_iter=zeros(sim_runs,9);

for sim_iter = 1:sim_runs

    J = 1;
    if model_sig == 1
        % Use wavelets on randomly perturbed standard signals from WaveLab
        %  Name string: 'Blocks', 'Piece-Polynomial' (Piece-Wise 3rd degree polynomial)

        % signal length
        N = 2048;
        % number of observations to make
        m = 1024;

        Name = 'Blocks';
        sig = MakeSignal(Name,N)';
        qmf = MakeONFilter('Haar');
        
        Name = 'Piece-Polynomial';
        sig = MakeSignal(Name,N)';
        qmf = MakeONFilter('Daubechies',8);

        wc_sig  = FWT_PO(sig,J,qmf);
        x = wc_sig;
    elseif model_sig == 2
        % Use strips from some image
        
        % I = imread('cameraman.tif');
        % I = imread('pirate.tif')';
        % I = imread('barb.gif');
        % I = (imread('Lena.jpg'));
        % load boats;
        % I = boats;
        % load blocks
        % I = blocks;
        I = rgb2gray(imread('house.tiff'));
        
        N = size(I,1);
        m = N/2;
        sig = double(I(:,sim_iter))/256;
        qmf = MakeONFilter('Haar');
        % qmf = MakeONFilter('Daubechies',4);
        wc_sig  = FWT_PO(sig,J,qmf);
        x = wc_sig;
    else
        disp('Do it yourself');
    end


    % measurement matrix
    G = randn(m,N)/sqrt(m);
    % A = (orth(G'))';
    A = G;

    % observations
    sigma = .00;
    e = randn(m,1)*sigma;
    y = A*x + e;
    %figure(101); clf; plot(y); hold on; plot(e,'r');

    % regularization parameter
    tau = 0.005*max(abs(A'*y));    
    % tau = sigma * sqrt(log(N)*2); % BPDN

    % Solve for measurements y
    tic
    [xp, gamma_x, iter_xp, t_homotopy1] = BPDN_homotopy_function(A, y, tau,4*N); %BPDN
    t_homotopy1 = toc;
    tolA_h = tau*sum(abs(xp))+1/2*(norm(A*xp-y))^2;

    % GPSR
    pdg_scale = 1e-8;
    t_BB_mono = 0;
    tolA = tolA_h+pdg_scale;
    stopCri=4;
    debias = 0;
    tic;
    [x_BB_mono,x_debias_BB_mono,obj_BB_mono,...
        times_BB_mono,debias_start_BB_mono,mse]= ...
        GPSR_BB(y,A,tau,...
        'Debias',debias,...
        'Monotone',1,...
        'Initialization',0,...
        'Continuation',0,...
        'StopCriterion',stopCri,...
        'Verbose', 0,...
        'MaxiterA', 1000,...
        'ToleranceA',tolA,...
        'ToleranceD',0.001);
    t_BB_mono = toc;

    opts.f_value_tol = tolA;
    M = [];
    tic
    n = N;
    tic;
    [x_FPC, Out_FPC] = FPC_AS(n,A,y,tau,M,opts); % call FPC_AS
    t_FPC = toc;
    
    % Model the change in signal
    if model_sig == 1
        % Perturb the signal slightly at random
        sig_mod = DynamicX_varying_model(Name,N)';
        wc_mod  = FWT_PO(sig_mod,J,qmf);
        xt = wc_mod;
        
    elseif model_sig == 2
        % Use strips from some image
        pm1 = randsrc;
        sig_mod = double(I(:, sim_iter+1))/256;
        wc_mod  = FWT_PO(sig_mod,J,qmf);
        xt = wc_mod;
    else
        disp('hehehe...');
    end
    dx = xt-x;
    yt = A*xt;

    pk_old = A'*(A*xp-y);
    AtAgx = A(:,gamma_x)'*A(:,gamma_x);
    iAtAgx = inv(AtAgx);

    if norm(y-yt)>2*eps
        tic;
        [xp_h, gamma_xh, iter_xp_update, t_homotopy_update] = DynamicX_BPDN_function(A, AtAgx, iAtAgx, y, yt, xp, gamma_x, pk_old, tau, 4*N);
        t_homotopy_update = toc;
    else
        xp_h = xp;
        gamma_xh = gamma_x;
        iter_xp_update = 0;
        t_homotopy_update = 0;
    end
    tolA_h = tau*sum(abs(xp_h))+1/2*(norm(A*xp_h-yt))^2;

    % GPSR
    tolA = tolA_h+pdg_scale;
    stopCri=4;
    debias = 0;
    tic;
    % PUTTING CONTINUATION IN HERE MAKES THINGS WORSE
    [x_BB_mono_update,x_debias_BB_mono,obj_BB_mono,...
        times_BB_mono,debias_start_BB_mono,mse]= ...
        GPSR_BB(yt,A,tau,...
        'Debias',debias,...
        'Monotone',1,...
        'Initialization',x_BB_mono,...
        'StopCriterion',stopCri,...
        'ToleranceA',tolA,...
        'MaxiterA', 1000,...
        'Verbose', 0,...
        'ToleranceD',0.001);
    t_BB_mono_update = toc;

    opts.x0 = x_FPC;
    opts.f_value_tol = tolA;
    M = [];
    tic
    [x_FPC_update, Out_FPC_update] = FPC_AS(n,A,yt,tau,M,opts); % call FPC_AS
    t_FPC_update = toc;
    
    stack_time(sim_iter,:) = [t_homotopy1 t_BB_mono t_FPC t_homotopy_update t_BB_mono_update t_FPC_update];
    % For old DynamicX_GPSR_Blocks_final or similar
    %stack_iter(sim_iter,:) = [iter_xp iter_xp_update norm(dx) norm(xt-xp_h)/norm(xt) norm(xt-x_BB_mono_update)/norm(xt) tau length(gamma_x)];

    % For new mat files
    stack_iter(sim_iter,:) = [iter_xp iter_xp_update length(times_BB_mono) (Out_FPC_update.nProdA+Out_FPC_update.nProdAt)/2 norm(xt-xp_h)/norm(xt) norm(xp_h-x_BB_mono_update)/norm(xt) norm(xt-x_FPC_update)/norm(xt) tau length(gamma_x)];
    % I_rec(:,sim_iter) = IWT_PO(xp,J,qmf);
end

% Generate table for performance
disp(['       DynamicX          BPDN homotopy           GPSR            FPC_AS']);
disp(['   nProdAtA    CPU     nProdAtA    CPU    nProdAtA   CPU    nProdAtA   CPU']);
outer_iterations = 1;
table_stats = zeros(1,8);
table_stats(:,3) = mean(reshape(stack_iter(:,1),sim_runs,outer_iterations));
table_stats(:,1) = mean(reshape(stack_iter(:,2),sim_runs,outer_iterations));
table_stats(:,5) = mean(reshape(stack_iter(:,3),sim_runs,outer_iterations));
table_stats(:,7) = mean(reshape(stack_iter(:,4),sim_runs,outer_iterations));
table_stats(:,4) = mean(reshape(stack_time(:,1),sim_runs,outer_iterations));
table_stats(:,2) = mean(reshape(stack_time(:,4),sim_runs,outer_iterations));
table_stats(:,6) = mean(reshape(stack_time(:,5),sim_runs,outer_iterations));
table_stats(:,8) = mean(reshape(stack_time(:,6),sim_runs,outer_iterations));
table_stats
