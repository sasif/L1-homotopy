% Dantzig selector sequential measurement update homotopy comparison
% Author: Salman Asif
% Created: February 2009:
% Modified: June 2009
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

clear; clc

% % load fixed random states
% load RandomStates
% randn('state',s_randn);
% rand('state',s_rand);

% signal length
N = 256;
% number of spikes to put down
T = 20;
% number of observations to make
K = 100;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = randn(T,1);%sign(randn(T,1)).*randint(T,1,256);

% measurement matrix
G = randn(K,N)/sqrt(K);
% G = orth(G')';
A = G(1:end-1,:);
b = G(end,:);

x(q(1:T)) = randsrc(T,1);

% observations
sigma = .01;
e = randn(K-1,1)*sigma;
y = A*x + e;

% regularization parameter
tau = 0.005*max(abs(A'*y)); %l1_ls
if sigma>0
    tau = sigma * sqrt(log(N)*2); % BPDN or DS
    %tau = max(abs(A'*(A*x-y))); % ideal ???
end

% Initial update for the DS
[xp, lame, gamma_x, gamma_lambda, xp_iter, t1] = DS_homotopy_function(A, y, tau, 4*N);


Aglgx = A(:,gamma_lambda)'*A(:,gamma_x); 
[Q_glgx R_glgx] = qr(Aglgx);

pk = A'*(A*xp-y);
ak = A'*(A*lame);
w = b*x+randn*sigma;

[xp_h, lambda_h, gamma_xh, gamma_lh, xp_h_iter, th] = DynamicSeq_DS_function(A, b, Q_glgx, R_glgx, y, w, xp, lame, gamma_x, gamma_lambda, pk, ak, tau, 4*N);

y2 = [y; w];
G = [A; b];
[xp2, lame2, gamma_x2, gamma_lambda2, xp2_iter, t2] = DS_homotopy_function(G, y2, tau, 4*N);

disp(' ');
disp('Results for ');
disp('old_homotopy , new_homotopy , update_homotopy')
cputime_comparison = [t1 t2 th]
Iter_comparison = [xp_iter xp2_iter xp_h_iter]
