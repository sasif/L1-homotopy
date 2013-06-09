% BPDN_homotopy_demo.m
% 
% Solves the following basis pursuit denoising problem
% min_x  \epsilon ||x||_1 + 1/2*||y-Ax||_2^2
%
% This script is a modified version of PD_pursuit method which
% uses only Primal update in PD_pursuit method.
%
% Here support of primal and dual vectors remain same at every step.
% The update direction here can be considered as equivalent to dual in Dantzig
% selector.

% In S-step solution property (for noiseless case) this becomes more clear when,
% if primal update direction and dual vector match at every step
% (maybe opposite in sign, depending on the formulation) then
% homotopy path taken by BPDN and Dantzig selector is exactly same.
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
%
%-------------------------------------------+
% Copyright (c) 2008.  Muhammad Salman Asif 
%-------------------------------------------+

clear; close all

rseed = 0;
rand('state',rseed);
randn('state',rseed);

N = 512;   % signal length
T = 50;    % sparsity level
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
thresh = 1e-3; %sqrt(2*log(N))*sigma;

% Initialization of primal and dual sign and support
z_x = zeros(N,1);
z_lambda = zeros(N,1);
gamma_lambda = [];  % Dual support
gamma_x = [];       % Primal support

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

% Turn the plots on and off
constraint_plots = 1; % 1 - draw the plots after every homotopy iteration.

AtglAgx = A(:,gamma_lambdak)'*A(:,gamma_xk);
iAtglAgx = inv(A(:,gamma_lambdak)'*A(:,gamma_xk));
AtgxAgl = AtglAgx';
iAtgxAgl = iAtglAgx';

while ~done
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')

    gamma_x = gamma_xk;
    gamma_lambda = gamma_lambdak;
    z_lambda = z_lambdak;
    z_x = z_xk;
    x_k = xk_1;
    lambda_k = lambdak_1;

    %%%%%%%%%%%%%%%%%%%%%
    %%%% update on x %%%%
    %%%%%%%%%%%%%%%%%%%%%

    % Update direction
%     del_x = -inv(A(:,gamma_lambda)'*A(:,gamma_x))*z_lambda(gamma_lambda);
%     diff_inv = [ max(max(abs(inv(A(:,gamma_lambda)'*A(:,gamma_x))-iAtglAgx)))]
%     if diff_inv > 1e-8
%         stp =1;
%     end
    del_x = -iAtglAgx*z_lambda(gamma_lambda);
    del_x_vec = zeros(N,1);
    del_x_vec(gamma_x) = del_x;

    pk = Primal_constrk;
    % dk = AtA*del_x_vec;
    dk = A'*(A(:,gamma_x)*del_x);
    
    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW. 
    pk_temp = Primal_constrk;
    gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,1e-12));
    pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
    
    xk_temp = x_k;
    gammaX_temp = find(abs(x_k)<1*eps);
    xk_temp(gammaX_temp) = 0;
    %%%---
    
    % Compute the step size
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
    
    if epsilon < thresh; %sqrt(2*log(N))*sigma; %1e-7 %|| iter > 5*T || (length(gamma_lambda) == K)
        delta_end = epsilon_old-thresh;
        Primal_constrk = pk+delta_end*dk;
        xk_1 = x_k + delta_end*del_x_vec;
        disp('done!');
        break;
    end
    
    if length(gamma_x)==M & chk_x ==0
        stp =1;
        disp('Cannot do it Sire'); % Commondos: are you out of your mind Sire!
        break;
    end

    if chk_x == 1
        % If an element is removed from gamma_x
        gl_old = gamma_lambda;
        gx_old = gamma_x;

        outx_index = find(gamma_x==out_x);
        gamma_x(outx_index) = gamma_x(end);
        gamma_x(end) = out_x;
        gamma_x = gamma_x(1:end-1);

        outl_index = outx_index;
        gamma_lambda = gamma_x;
        gamma_lambdak = gamma_lambda;
        gamma_xk = gamma_x;

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
        AtglAgx = AtgxAgl;
        iAtgxAgl = update_inverse(AtgxAgl_ij, iAtgxAgl_ij,2);
        iAtglAgx = iAtgxAgl;
        xk_1(out_x) = 0;
    else
        % If an element is added to gamma_x
        gamma_lambdak = [gamma_lambda; i_delta];
        gamma_xk = gamma_lambdak;
        i_theta = i_delta;
        new_lambda = i_delta;
        AtgxAnl = A(:,gamma_x)'*A(:,new_lambda);
        AtglAgx_mod = [AtglAgx A(:,gamma_lambda)'*A(:,i_theta); AtgxAnl' A(:,new_lambda)'*A(:,i_theta)];

        AtglAgx = AtglAgx_mod;
        AtgxAgl = AtglAgx;
        iAtglAgx = update_inverse(AtglAgx, iAtglAgx,1);
        iAtgxAgl = iAtglAgx;
        xk_1(i_theta) = 0;
    end

    z_lambdak = zeros(N,1);
    z_lambdak(gamma_lambdak) = sign(Primal_constrk(gamma_lambdak));
    z_xk = -z_lambdak;
    Primal_constrk(gamma_lambdak) = sign(Primal_constrk(gamma_lambdak))*epsilon;
    
    if iter > 500*T
        disp('too many iterations ooooooooops');
        break;
    end

    if constraint_plots
        fig = figure(1);
        subplot(2,1,1)
        hold off
        plot(pk,'.r', 'MarkerSize',14);
        hold on;
        plot(Primal_constrk, 'LineWidth',1);

        if chk_x == 0
            plot(new_lambda, Primal_constrk(new_lambda),'or','MarkerSize',18,'LineWidth',2);
            text(new_lambda, Primal_constrk(new_lambda)*1.1, ['Incoming \gamma = ',num2str(new_lambda)],'FontSize',14);
        else
            plot(out_x, Primal_constrk(out_x),'*k','MarkerSize',18,'LineWidth',2);
            text(out_x,Primal_constrk(out_x)*1.1, ['Outgoing \gamma = ',num2str(out_x)],'FontSize',14);
        end
        set(gca,'FontSize',16, 'XLim',[1 N] );
        title({'BPDN shrinkage constraints,'; ['n = ',num2str(N), ', m = ', num2str(M), ', T = ',num2str(T)]});
        plot(1:N, epsilon*ones(1,N),'--k','MarkerSize',12);
        plot(1:N, -epsilon*ones(1,N), '--k','MarkerSize',12);
        plot(1:N, epsilon_old*ones(1,N),'--m','MarkerSize',12);
        plot(1:N, -epsilon_old*ones(1,N), '--m','MarkerSize',12);

        subplot(2,1,2)
        hold off
        plot(x_k,'.r','MarkerSize',14); hold on;
        plot(xk_1,'LineWidth',1);
        if chk_x == 1
            plot(out_x, 0,'ok', 'MarkerSize',18,'LineWidth',2);
        else
            plot(new_lambda, 0,'or', 'MarkerSize',18,'LineWidth',2);
        end
        set(gca,'FontSize',16,'XLim',[1 N]);
        title(['BPDN estimate at \tau = ',num2str(epsilon), ', iter = ', num2str(iter)]);
        
        if iter == 1
            disp('  ');
            disp('Every frame in the figure corresponds to a critical point on the homotopy path.')
            disp('Circle represents an incoming element, star represents an outgoing element.');
            disp(' ');
            disp('Put pause somewhere in the code to see this. ');
            disp('For now press some key to continue...');
            pause
        end

        %drawnow;
        %print(gcf,'-dbmp','-r200',['BPDNPath_',num2str(iter)])
    end
    % pause
    [iter epsilon delta]
end

fig = figure(1);
subplot(2,1,1)
hold off
plot(pk,'.r', 'MarkerSize',14);
hold on;
plot(Primal_constrk, 'LineWidth',1);


set(gca,'FontSize',16, 'XLim',[1 N] );
title({'BPDN shrinkage constraints,'; ['n = ',num2str(N), ', m = ', num2str(M), ', T = ',num2str(T)]});
plot(1:N, epsilon*ones(1,N),'--k','MarkerSize',12);
plot(1:N, -epsilon*ones(1,N), '--k','MarkerSize',12);
plot(1:N, epsilon_old*ones(1,N),'--m','MarkerSize',12);
plot(1:N, -epsilon_old*ones(1,N), '--m','MarkerSize',12);

figure(1);
subplot(2,1,2)
hold off
plot(x,'.r','MarkerSize',14); hold on;
plot(xk_1,'LineWidth',1);
set(gca,'FontSize',16,'XLim',[1 N]);
title(['BPDN estimate at \tau = ',num2str(thresh), ', iter = ', num2str(iter)]);
