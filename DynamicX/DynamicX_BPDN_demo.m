% BPDN/LASSO update homotopy comparison
% Author: Muhammad Salman Asif
% Created: February 2009:
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

clear; clc;

% % load fixed random states
% load RandomStates
% randn('state',s_randn);
% rand('state',s_rand);

% signal length
N = 512;
% number of spikes to put down
T = 40;
% number of observations to make
M = 350;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = randn(T,1);%sign(randn(T,1)).*randint(T,1,256);

% measurement matrix
A = randn(M,N)/sqrt(M);

% observations
sigma = .01;
e = randn(M,1)*sigma;
y = A*x + e;

% regularization parameter
tau = 0.005*max(abs(A'*y)); %l1_ls
if sigma>0
    tau = sigma * sqrt(log(N)*2); % BPDN
    %tau = max(abs(A'*(A*x-y))); % ideal ???
end

% Initial update using BPDN
[xp, gamma_x, xp_iter, t1] = BPDN_homotopy_function(A, y, tau,4*N); %Lasso

% Model the change in signal
dx_model = 1; % 0 - Perturb only non-zero locations
              % 1 - Add some new elements (T_in) and remove some existing
              % elements (T_out)
switch dx_model
    case  1
        % Add some elements (T_in) or remove some existing elements (T_out)
        dx  = zeros(N,1);
        dx(q(1:T)) = randn*.1;
        T_in = round(T/10);
        T_out = round(T/10);
        qin = randperm(N);

        dx(qin(1:T_in)) = (randn(T_in,1));

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

% Gram matrix update
AtAgx = A(:,gamma_x)'*A(:,gamma_x);
iAtAgx = inv(AtAgx);
pk = A'*(A*xp-y);

[xp_h, gamma_xh, xp_h_iter, th] =  DynamicX_BPDN_function(A, AtAgx, iAtAgx, y, yt, xp, gamma_x, pk, tau, 4*N);

[xp2, gamma_x2, xp2_iter, t2] = BPDN_homotopy_function(A, yt, tau, 4*N); %Lasso

% [old_homotopy new_homotopy update_homotopy]
Cputime_comparison = [t1 t2 th]
Iteration_comparison = [xp_iter xp2_iter xp_h_iter]

% figure(2); plot(xp_h-xp2); shg
