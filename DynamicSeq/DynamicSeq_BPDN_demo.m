% BPDN update with sequential measurements homotopy comparison
% Author: Salman Asif
% Created: February 2009:
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

clear; clc

% load fixed random states
load RandomStates
randn('state',s_randn);
rand('state',s_rand);

% signal length
N = 256;
% number of spikes to put down
T = 30;
% number of observations to make
K = 100;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = randn(T,1);%sign(randn(T,1)).*randint(T,1,256);

% measurement matrix
G = randn(K,N)/sqrt(K);
% G = orth(G')';
A = G; %(1:end-1,:);
% b = G(end,:);

% observations
sigma = .01;
e = randn(K,1)*sigma;
y = A*x + e;

% regularization parameter
tau = 0.005*max(abs(A'*y)); %l1_ls
if sigma>0
    tau = sigma * sqrt(log(N)*2); % BPDN
    %tau = max(abs(A'*(A*x-y))); % ideal ???
end

% Initial update using Lasso
[xp, gamma_x, xp_iter, t1] = BPDN_homotopy_function(A, y, tau,4*N); %Lasso

chk_e0 = 1;  % 0 or 1: This selects if we want to take U = (A'A+b'b) (1) or (A'A+e0*b'b) (0)...!
% In simulation U = (A'A+b'b) works better

b = randn(1,N);
ew = randn(1,1)*sigma;
w = b*x+ew;

% Gram matrix update
AtAgx = A(:,gamma_x)'*A(:,gamma_x)+chk_e0*b(:,gamma_x)'*b(:,gamma_x);
iAtAgx = inv(AtAgx);
pk = A'*(A*xp-y);

[xp_h, gamma_xh, xp_h_iter, th] =  DynamicSeq_BPDN_function(A, b, AtAgx, iAtAgx, y, w, xp, gamma_x, pk, tau, chk_e0, 4*N);


y2 = [y; w];
G = [A; b];
[xp2, gamma_x2, xp2_iter, t2] = BPDN_homotopy_function(G, y2, tau, 4*N); %Lasso

% [old_homotopy new_homotopy update_homotopy]
disp(' ');
disp('Results for ');
disp('old_homotopy , new_homotopy , update_homotopy')
cputime_comparison = [t1 t2 th]
Iter_comparison = [xp_iter xp2_iter xp_h_iter]

% figure(2); plot(xp_h-xp2); shg
