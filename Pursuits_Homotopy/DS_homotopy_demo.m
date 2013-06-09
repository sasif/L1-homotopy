% DS_homotopy_demo_test.m
%
% Solves the following Dantzig selector problem
% min_x  ||x||_1  subject to  ||A'(Ax-y)||_\infty <= epsilon
%
% using Primal Dual Pursuit Homotopy method.
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
%
%-------------------------------------------+
% Copyright (c) 2008.  Muhammad Salman Asif 
%-------------------------------------------+

clear; close all; 

rseed = 0;
rand('state',rseed);
randn('state',rseed);

N = 512;   % signal length
T = 20;    % sparsity level
M = 250;    % no. of measurements

% Generate a random signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1));

% measurement matrix
A = randn(M,N)/sqrt(M); % Gaussian

% Hadamard
% H = hadamard(N);
% A = H(q(1:M),:)/sqrt(M);

% Bernoulli
% A = randsrc(M,N)/sqrt(M);

% Random Projection
% A = orth(A')';

% AtA = A'*A;

% measurements
SNR = inf;
s = A*x;
sigma = sqrt((norm(s)^2/(10^(SNR/10)))/M);
% sigma = 0.01;
e = randn(M,1)*sigma;
y = A*x+e;
thresh = 1e-3; % sqrt(2*log(N))*sigma;

% Initialization of primal and dual sign and support
z_x = zeros(N,1);
z_lambda = zeros(N,1);
gamma_lambda = [];  % Dual support
gamma_x = [];       % Primal support
gamma_lambda_app = [];

% Initial step
Primal_constrk = -A'*y;
[c i] = max(abs(Primal_constrk));

gamma_lambdak = i;
gamma_xk = gamma_lambdak;

z_lambda(gamma_lambdak) = sign(Primal_constrk(gamma_lambdak));
epsilon = c;
Primal_constrk(gamma_lambdak) = sign(Primal_constrk(gamma_lambdak))*epsilon;
xk_1 = zeros(N,1);

lambdak_1 = zeros(N,1);
lambdak_1(gamma_lambdak) = inv(A(:,gamma_lambdak)'*A(:,gamma_lambdak))*z_lambda(gamma_lambdak);

Dual_constrk = A'*(A*lambdak_1);

z_x(gamma_xk) = -sign(Dual_constrk(gamma_xk));
Dual_constrk(gamma_xk) = sign(Dual_constrk(gamma_xk));

z_xk = z_x;
z_lambdak = z_lambda;

% loop parameters
done = 0;
iter = 0;
data_precision = eps;   % floating point precision

old_delta = 0;
out_lambda = [];
out_x = [];
count_delta_stop = 0;

constraint_plots = 0;

AtglAgx = A(:,gamma_lambdak)'*A(:,gamma_xk);
iAtglAgx = inv(A(:,gamma_lambdak)'*A(:,gamma_xk));
AtgxAgl = AtglAgx';
iAtgxAgl = iAtglAgx';
time_s = clock;

itr_history = [iter epsilon 0 0 length(gamma_xk)]

while ~done
    iter = iter+1;
    warning('off','MATLAB:divideByZero')

    %%%% IDEA implementation %%%%
    gamma_x = gamma_xk;
    gamma_lambda = gamma_lambdak;
    z_lambda = z_lambdak;
    z_x = z_xk;
    x_k = xk_1;
    lambda_k = lambdak_1;

    %%%%%%%%%%%%%%%%%%%%%
    %%%% update on x %%%%
    %%%%%%%%%%%%%%%%%%%%%

    % UPDATE INVERSE MATRIX

    % Update direction
    %del_x = -inv(A(:,gamma_lambda)'*A(:,gamma_x))*z_lambda(gamma_lambda);
    del_x = -iAtglAgx*z_lambda(gamma_lambda);
    del_x_vec = zeros(N,1);
    del_x_vec(gamma_x) = del_x;

    pk = Primal_constrk;
    % dk = AtA*del_x_vec;
    dk = A'*(A(:,gamma_x)*del_x);
    
    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW. 
    pk_temp = Primal_constrk;
    % gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,1e-12));
    % pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
    
    xk_temp = x_k;
    % gammaX_temp = find(abs(x_k)<2*eps);
    % xk_temp(gammaX_temp) = 0;
    %%%---
    
    % Step size computation
    [i_delta, out_x, delta, chk_x] = update_primal(gamma_x, gamma_lambda, z_x,  xk_temp, del_x_vec, pk_temp, dk, epsilon, out_lambda);

    if old_delta < 4*eps & delta < 4*eps
        count_delta_stop = count_delta_stop + 1;
    else
        count_delta_stop = 0;
    end
    if count_delta_stop >= 50
        disp('stuck in some corner');
        break;
    end
    old_delta = delta;

    xk_1 = x_k+delta*del_x_vec;
    Primal_constrk = pk+delta*dk;
    epsilon_old = epsilon;
    epsilon = epsilon-delta;

    if chk_x == 1
        % If an element is removed from gamma_x
        out_x = out_x(1);
        gl_old = gamma_lambda;
        gx_old = gamma_x;

        outx_index = find(gamma_x==out_x);
        gamma_x(outx_index) = gamma_x(end);
        gamma_x(end) = out_x;
        gamma_x = gamma_x(1:end-1);

        %CHECK THIS SINGULARITY CONDITION USING SCHUR COMPLEMENT IDEA !!!
        % X = [A B; C D];
        % detX = detA detS
        % S = D-C A^{-1} B

        % So the idea here is slightly different, we don't want the last
        % entry (Q22) in the inverse of modified Gram matrix to be zero, because
        % that will make inversion impossible.

        % This is the best thing, pick the index with maximum absolute
        % value in the column corresponding to out_x in iAtgxgl. This will
        % ensure that new matrix will be invertible.
        [max_val outl_index] = max(abs(iAtgxAgl(:,outx_index)));
        new_lambda = gamma_lambda(outl_index);
        some_lambda = new_lambda;

        outl_index = find(gamma_lambda==some_lambda);
        gamma_lambda(outl_index) = gamma_lambda(end);
        gamma_lambda(end) = some_lambda;
        gamma_lambdak = gamma_lambda;
        gamma_lambda = gamma_lambda(1:end-1);

        rowi = outx_index; % ith row of A is swapped with last row (out_x)
        colj = outl_index; % jth column of A is swapped with last column (out_lambda)
        AtgxAgl_ij = AtgxAgl;
        temp_row = AtgxAgl_ij(rowi,:);
        AtgxAgl_ij(rowi,:) = AtgxAgl_ij(end,:);
        AtgxAgl_ij(end,:) = temp_row;
        temp_col = AtgxAgl_ij(:,colj);
        AtgxAgl_ij(:,colj) = AtgxAgl_ij(:,end);
        AtgxAgl_ij(:,end) = temp_col;
        iAtgxAgl_ij = iAtgxAgl;
        temp_row = iAtgxAgl_ij(colj,:);
        iAtgxAgl_ij(colj,:) = iAtgxAgl_ij(end,:);
        iAtgxAgl_ij(end,:) = temp_row;
        temp_col = iAtgxAgl_ij(:,rowi);
        iAtgxAgl_ij(:,rowi) = iAtgxAgl_ij(:,end);
        iAtgxAgl_ij(:,end) = temp_col;

        AtgxAgl = AtgxAgl_ij(1:end-1,1:end-1);
        AtglAgx = AtgxAgl';
        iAtgxAgl = update_inverse(AtgxAgl_ij, iAtgxAgl_ij,2);
        iAtglAgx = iAtgxAgl';
        xk_1(out_x) = 0;
    else
        % If a new element is added to gamma_lambda
        gamma_lambdak = [gamma_lambda; i_delta];
        new_lambda = i_delta;
        lambda_k(new_lambda) = 0;
    end

    out_lambda = [];
    z_lambdak = zeros(N,1);
    z_lambdak(gamma_lambdak) = sign(Primal_constrk(gamma_lambdak));
    sgn_new_lambda = sign(Primal_constrk(new_lambda));
    Primal_constrk(gamma_lambdak) = sign(Primal_constrk(gamma_lambdak))*epsilon;

    if epsilon < thresh; %sqrt(2*log(N))*sigma; %1e-7 %|| iter > 5*T || (length(gamma_lambda) == K)
        delta_end = epsilon_old-thresh;
        Primal_constrk = pk+delta_end*dk;
        xk_1 = x_k + delta_end*del_x_vec;
        disp('done!');
        break;
    end

    if iter > 500*T
        disp('too many iterations ooooooooops');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% update on lambda %%%
    %%%%%%%%%%%%%%%%%%%%%%%%

    sgn_new_lambda = sign(Primal_constrk(new_lambda));

    AtgxAnl = A(:,gamma_x)'*A(:,new_lambda);
    %AtgxAnl = AtA(gamma_x,new_lambda);

    % Update direction
    %del_lambda = -inv(A(:,gamma_x)'*A(:,gamma_lambda))*(A(:,gamma_x)'*A(:,new_lambda));
    del_lambda = -iAtgxAgl*AtgxAnl;
    del_lambda_p = zeros(N,1);
    del_lambda_p(gamma_lambda) = del_lambda*sgn_new_lambda;
    del_lambda_p(new_lambda) = 1*sgn_new_lambda;

    ak = Dual_constrk;
    % bk = AtA*del_lambda_p;
    bk = A'*(A(:,[gamma_lambda; new_lambda])*del_lambda_p([gamma_lambda; new_lambda]));

    % check if the sign of update direction is correct
    if sign(bk(out_x)) == sign(ak(out_x)) & abs(bk(out_x)) >=eps
        bk = -bk;
        del_lambda_p = -del_lambda_p;
    end

    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW. 
    ak_temp = Dual_constrk;
    % gammaX_temp = find(abs(abs(Dual_constrk)-1)<2*eps);
    % ak_temp(gammaX_temp) = sign(Dual_constrk(gammaX_temp));
    
    lambdak_temp = lambda_k;
    % gammaL_temp = find(abs(lambda_k)<1*eps);
    % lambdak_temp(gammaL_temp) = 0;
    %%%---
    
    % Step size computation
    [i_theta, out_lambda, theta, chk_lambda] = update_dual(gamma_x, gamma_lambda, z_lambda, lambdak_temp, del_lambda_p, ak_temp, bk, new_lambda, out_x);

    out_x = [];

    lambdak_1 = lambda_k + theta*del_lambda_p;
    Dual_constrk = ak+theta*bk;

    if chk_lambda == 1
        % If an element is removed from gamma_lambda
        outl_index = find(gamma_lambda==out_lambda);
        i_theta = gamma_x(end);

        if ~isempty(outl_index) % because otherwise new_lambda is removed.
            % Use Sherman-Woodburry-Morrison formula
            % (A+BC)^{-1} = A^{-1} - A^{-1}B(I+CA^{-1}B)^{-1}CA^{-1};
            iA = iAtgxAgl;
            B = AtgxAnl-A(:,gamma_x)'*A(:,out_lambda);
            C = zeros(1,length(gamma_lambda));
            C(outl_index) = 1;
            AtgxAgl(:,outl_index) = AtgxAnl;
            AtglAgx = AtgxAgl';

            iAB = iA*B;
            CiA = iA(outl_index,:);
            iAtgxAgl = iA-iAB*((CiA)/(1+iAB(outl_index)));
            iAtglAgx = iAtgxAgl';

            gamma_lambda(outl_index) = new_lambda;
        end

        z_xk = zeros(N,1);
        z_xk(gamma_x) = -sign(Dual_constrk(gamma_x));
        lambdak_1(out_lambda) = 0;
    else
        % If an element is added to gamma_x

        %%% adding columns
        % AtglAgx_mod = zeros(length(gamma_x));
        AtglAgx_mod = [AtglAgx A(:,gamma_lambda)'*A(:,i_theta); AtgxAnl' A(:,new_lambda)'*A(:,i_theta)];
        %AtglAgx_mod = [AtglAgx AtA(gamma_lambda,i_theta); AtgxAnl' AtA(new_lambda,i_theta)];

        %CHECK THIS SINGULARITY CONDITION USING SCHUR COMPLEMENT IDEA !!!
        % X = [A B; C D];
        % detX = detA detS
        % S = D-C A^{-1} B
        A11 = AtglAgx;
        A12 = AtglAgx_mod(1:end-1,end);
        A21 = AtglAgx_mod(end,1:end-1);
        A22 = AtglAgx_mod(end,end);
        S = A22 - A21*(iAtglAgx*A12);
        if abs(S) < 2*eps
            disp('matrix has become singular');
            % done = 1;
            % break;
        end
        AtglAgx = AtglAgx_mod;
        AtgxAgl = AtglAgx';
        iAtglAgx = update_inverse(AtglAgx, iAtglAgx,1);
        iAtgxAgl = iAtglAgx';

        out_lambda = [];
        gamma_lambda = [gamma_lambda; new_lambda];

        gamma_x = [gamma_x; i_theta];
        new_x = i_theta;
        z_xk = zeros(N,1);
        z_xk(gamma_x) = -sign(Dual_constrk(gamma_x));
        xk_1(new_x) = 0;
    end
    Dual_constrk(gamma_x) = sign(Dual_constrk(gamma_x));
    gamma_lambdak = gamma_lambda;
    gamma_xk = gamma_x;
   
    if constraint_plots
        fig = figure(1);
        subplot(211)
        hold off;
        plot(pk,'.r', 'MarkerSize',14);        
        hold on;
        plot(Primal_constrk, 'LineWidth',1);
        set(gca,'FontSize',16, 'XLim',[1 N] );
        title({'BPDN shrinkage constraints,'; ['n = ',num2str(N), ', m = ', num2str(M), ', T = ',num2str(T)]});
        plot(1:N, epsilon*ones(1,N),'--k','MarkerSize',12);
        plot(1:N, -epsilon*ones(1,N), '--k','MarkerSize',12);
        plot(1:N, epsilon_old*ones(1,N),'--m','MarkerSize',12);
        plot(1:N, -epsilon_old*ones(1,N), '--m','MarkerSize',12);
        
        
        subplot(212);
        hold off;
        plot(ak,'.r');
        hold on;
        plot(Dual_constrk);
        %legend('Old constraint', 'New constraint', 'Element change','Location','NorthEast');
        title('Dual constraints');
        plot(1:N, 1, '--k');
        plot(1:N, -1, '--k');
        axis([1 N  -1.5 1.5])
        
        if chk_x == 0
            subplot(211);
            plot(new_lambda, Primal_constrk(new_lambda),'or','MarkerSize',18,'LineWidth',2);
            text(new_lambda, Primal_constrk(new_lambda)*1.1, ['Incoming \gamma = ',num2str(new_lambda)],'FontSize',14);
        else
            subplot(212)
            plot(out_x, (ak(out_x)),'*k','MarkerSize',18,'LineWidth',2);
            text(out_x,ak(out_x)*1.1, ['Outgoing \gamma = ',num2str(out_x)],'FontSize',14);
        end
        
        if chk_lambda == 0
            subplot(212)
            plot(new_x, Dual_constrk(new_x),'or','MarkerSize',18,'LineWidth',2);
            text(new_x, Dual_constrk(new_x)*1.1, ['Incoming \gamma = ',num2str(new_x)],'FontSize',14);
        else
            subplot(211);
            plot(out_lambda, (Primal_constrk(out_lambda)),'*k','MarkerSize',18,'LineWidth',2);
            text(out_lambda,Primal_constrk(out_lambda)*1.1, ['Outgoing \gamma = ',num2str(out_x)],'FontSize',14);
        end
    end
    itr_history_temp = [iter epsilon delta theta length(gamma_xk)]
    itr_history = [itr_history; itr_history_temp];
end
time_e = clock;
etime(time_e,time_s);

figure(1); clf; 
subplot(211); plot(pk,'.r'); hold on;
plot(Primal_constrk); 
plot(1:N, epsilon,'k'); plot(1:N, -epsilon, 'k'); plot(1:N, epsilon_old,'g'); plot(1:N, -epsilon_old, 'g');
title('Primal constraints');

subplot(212); plot(ak,'.r'); hold on; 
plot(Dual_constrk); 
plot(1:N, 1, '-k'); plot(1:N, -1, '-k'); axis([1 N  -1.5 1.5])
title('Dual constraints');

figure(2); clf; 
plot(xk_1); hold on; plot(x,'.r');
legend('estimation', 'original');
title('Estimated signal')

% %%
% x0 = A'*inv(A*A')*y;
% % large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1eq_pd(x0, Afun, Atfun, y, 1e-8, 1000, 1e-8, 1000);
% toc