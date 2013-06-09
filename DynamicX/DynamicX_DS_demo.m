% Test file for Dantzig selector (DS) update with dynamic change in x

% Author: Muhammad Salman Asif @ Georgia Tech
% E-mail: sasif@ece.gatech.edu
% Created: February 2009
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

clear; clc;

addpath('../Pursuits_Homotopy')

% % load fixed random states
% load RandomStates
% rand('state', s_rand);
% randn('state', s_randn);

% signal length
N = 512;
% number of spikes to put down
T = 20;
% number of observations to make
M = 256;

% random signal supported on q(1:T)
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1));


% measurement matrix
A = randn(M,N)/sqrt(M);

% observations
sigma = .01;
e = randn(M,1)*sigma;
y = A*x + e;

% regularization parameter
tau = 0.05*max(abs(A'*y)); %l1_ls
if sigma>0
    tau = sigma * sqrt(log(N)*2); % BPDN
    %tau = max(abs(A'*(A*x-y))); % ideal ???
end

% Solve for measurements y
[xp, lame, gamma_x, gamma_lambda, xp_iter, t1] = DS_homotopy_function(A, y, tau,4*N); %Lasso

% Model the change in signal
dx_model = 1; % 0 - Perturb only non-zero locations
              % 1 - Add some new elements (T_in) and remove some existing
              % elements (T_out)
switch dx_model
    case  1
        % Add some elements (T_in) or remove some existing elements (T_out)
        dx  = zeros(N,1);
        dx(q(1:T)) = randn*.1;
        T_in = round(T/20);
        T_out = round(T/20);
        qin = randperm(N);

        dx(qin(1:T_in)) = sign(randn(T_in,1));

        qout = randperm(T);
        oldx_indices = q(1:T);
        dx(oldx_indices(qout(1:T_out))) = -x(oldx_indices(qout(1:T_out)));

    case  2;
        % Perturb the current locations by Gaussian noise
        oldx_indices = q(1:T);
        dx = zeros(N,1);
        dx(oldx_indices) = randn(T,1)*.2;
    otherwise
        disp('Nooooo');
end

xt = x+dx;
yt = y+A*dx;

pk_old = A'*(A*xp-y);
ak_old = A'*(A*lame);

Aglgx = A(:,gamma_lambda)'*A(:,gamma_x); 
[Q_glgx R_glgx] = qr(Aglgx);

[xp_h, lambda_h, gamma_xh, gamma_lh, xp_h_iter, th] = DynamicX_DS_function(A, Q_glgx, R_glgx, y, yt, xp, lame, gamma_x, gamma_lambda, pk_old, ak_old, tau, 4*N);

[xp2, lame2, gamma_x2, gamma_lambda2, xp2_iter, t2] = DS_homotopy_function(A, yt, tau, 4*N);

disp(' ');
disp('old_homotopy , new_homotopy , update_homotopy')
cputime_comparison = [t1 t2 th]
iter_comparison = [xp_iter xp2_iter xp_h_iter]

% figure(2); plot(xp_h-xp2); shg
