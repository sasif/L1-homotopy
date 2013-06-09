% DynamicSeq_REC_demo.m
% Robust ell_1 decoding with homotopy update scheme for new measurements.
%
% Created by Salman Asif @ Georgia Tech
% February 2009
%
% Model changes in the error patterns with new measurements and observe the
% effect on the solution.
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

clear
close all

% % load fixed random states
% load RandomStates
% rand('state', s_rand);
% randn('state', s_randn);
  
% % source length
N = 256;
% codeword length
M = 3*N;

% number of perturbations
T = round(.2*M);

% coding matrix
Orth_mat = randn(M,M);
G = Orth_mat(:,1:N);
A = randn(M,N)/sqrt(M);

% Annihilating projection matrix 
AtA = A'*A;
iAtA = inv(AtA);
AiAtA = A*iAtA;
AiAtAAt = AiAtA*A';
Q = eye(M)-AiAtAAt;

% source word
x = randn(N,1);

% channel: perturb T randomly chosen entries
q = randperm(M);
e = zeros(M,1);
e(q(1:T)) = randsrc(T,1);

% Small noise
x0 = randn(N,1);
Ax0 = A*x0;
sigma = median(abs(Ax0))/20; % control the power in small noise
q_y = randn(M,1)*sigma;

% Received codeword
y = A*x+e+q_y;
figure(1); clf; plot([y e q_y]); %hold on; plot(e,'g-.'); plot(q_y,'r--');
title('Codeword, sparse errors and small noise');
legend('Codeword', 'Sparse', 'Noise'); shg

% Regularization parameter
tau = 0.01*max(abs(Q*y)); % l1_ls
% if sigma>0
%     tau = sigma * sqrt(log(N)*2); % BPDN
%     % tau = max(abs(Q'*q_y)); % ideal ???
% end

% Data recovery
t1_s = cputime;
Qy = Q*y;
[ep, gamma_e, ep_iter] =  BPDN_homotopy_function(Q, Qy, tau, 4*M);
t1_e = cputime;
t1 = t1_e-t1_s;
% xp = inv(A'*A)*A'*(y-ep);
xp = AiAtA'*(y-ep);

figure(2); clf; plot(ep); hold on; plot(e,'ro');
title(['Estimated and original sparse errors with ', num2str(M),' measurements']);
legend('Estimated', 'Original'); shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup for adding m_u new observations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No. of new observations
m_u =10;

b = randn(m_u,N)/sqrt(M); % new rows in coding matrix
T_new = randint(1,1,m_u/5); % new gross/sparse errors
d = zeros(m_u,1);
q_new = randperm(m_u);
d(q_new(1:T_new)) = randsrc(T_new,1);
q_w = randn(m_u,1)*sigma;

w = b*x+d+q_w;
F = [A; b];
s = [y;w];
c = [e; d];
q_yw = [q_y; q_w];

iAtAbt = iAtA*b';
biAtAbt = b*iAtAbt;
S_biAtAbt = inv(eye(m_u)+biAtAbt);
AiAtAbt = AiAtA*b';
iAtAbt_S = (iAtAbt*S_biAtAbt);
AiAtAbt_S = A*iAtAbt_S;

FtF = AtA+b'*b; %F'*F;
iFtF = iAtA - iAtAbt_S*iAtAbt'; % inv(FtF);
FiFtF = [AiAtA - AiAtAbt_S*iAtAbt'; iAtAbt'-biAtAbt*S_biAtAbt*iAtAbt']; % F*iFtF;
FiFtFFt = [AiAtAAt - AiAtAbt_S*AiAtAbt' AiAtAbt-AiAtAbt_S*biAtAbt; AiAtAbt'-biAtAbt*AiAtAbt_S' biAtAbt-biAtAbt*S_biAtAbt*biAtAbt]; %FiFtF*F';
P = eye(M+m_u)-FiFtFFt;

dp = w-b*xp;
z_d = sign(dp);
cp_h = [ep; dp];
gamma_n = M+find(abs(dp)>2*eps);
gamma_n_old = gamma_n;
gamma_h = [gamma_e; gamma_n];
epsilon = 0;
e0 = 0;

QgQ = Q(gamma_e,gamma_e);
PgP = P(gamma_h,gamma_h); 
uQ = AiAtAbt_S(gamma_e,:);
vQ = AiAtAbt(gamma_e,:);
QgQ_update = QgQ + uQ*vQ';
PgP_update = [[QgQ_update; P(gamma_n,gamma_e)] P(gamma_h,gamma_n)];
PgP = PgP_update;

iQgQ = inv(QgQ);
iQgQ_update = iQgQ - (iQgQ*uQ)*(inv(eye(m_u)+vQ'*iQgQ*uQ))*(vQ'*iQgQ);

P11 = QgQ_update;
P12 = P(gamma_e,gamma_n);
P21 = P12';
P22 = P(gamma_n,gamma_n);
S_P = inv(P22-P21*iQgQ_update*P12);
iP11_P12 = iQgQ_update*P12;
iPgP_update = [iQgQ_update+(iP11_P12*S_P)*iP11_P12' -iP11_P12*S_P; -S_P*iP11_P12' S_P];
iPgP = iPgP_update;

pk_old = P*(cp_h-s); % last indices will be zero.

% REC homotopy
[cp_h gamma_h cp_h_iter th] = DynamicSeq_REC_function(P, s, iPgP, cp_h, gamma_h, gamma_n, pk_old, tau, M, m_u, 4*M);
figure(3); clf; plot(cp_h); hold on; plot(c,'ro');
title(['Estimated and original sparse errors after new ', num2str(m_u),' measurements']);
legend('Estimated', 'Original'); shg

% % Check the solution using cvx
% cvx_begin
%     cvx_precision high
%     variables cp_cvx(M+m_u) zp_cvx(M+m_u) ;
%     minimize( tau*(norm(cp_cvx,1))+ .5*(sum(square(zp_cvx))))
%     subject to 
%     zp_cvx == P*(cp_cvx-s)
% cvx_end
% 
% % figure(103); clf; plot(cp_h); hold on; plot(cp2,'kx'); plot([cp_cvx; dp_cvx],'r.')
% figure(10); hold on; plot(cp_cvx-cp_h); shg

% Check the solution using homotopy from scratch.
[cp2, gamma_c2, cp2_iter, t2] = BPDN_homotopy_function(P, P*s, tau, 4*M);


pk2 = P*(cp2-s);
pk_h = P*(cp_h-s);

% [old_homotopy new_homotopy update_homotopy]
disp(' ');
disp('old_homotopy , new_homotopy , update_homotopy')
time_table = [t1 t2 th]
iter_table = [ep_iter cp2_iter cp_h_iter]

figure(4); subplot(2,1,1); hold on; plot(cp_h-cp2);
title('Difference between the solution of BPDN homotopy and REC update')
subplot(2,1,2); hold on; plot(pk2-pk_h); drawnow
title('Difference between the primal constraints of BPDN homotopy and REC update')
