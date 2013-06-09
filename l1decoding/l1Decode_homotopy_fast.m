% l1Decode_homotopy_fast.m
%
% Test out l1 decoding with homotopy update scheme
% Keeps adding new measurements until the original signal is recovered
% exactly.
%
% In this script, we start the decoding with $m=n$ measurements and
% sequentially add m_u new measurements and solve using homotopy.
%
% Created by Salman Asif @ Georgia Tech
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+


% Modification history :
% 09/29/08
% 10/02/08 : multiple measurements
% 02/15/09 : fast update of inverse matrix

close all
clear; clc

% % load fixed random states
% load RandomStates
% rand('state', s_rand);
% randn('state', s_randn);

% source length
N = 256;

% codeword length
M = 4*N;

% number of sparse errors
T = round(.2*N);

% Coding matrix
Orth_mat = randn(M,M);
Orth_mat = randsrc(M,M);
% Orth_mat = orth(Orth_mat);
G_mat = Orth_mat(:,1:N);

% Homotopy setup
m_st = N;   % number of measurements to start with
A = G_mat(1:m_st,:);
B_stack = G_mat(m_st+1:end,:);

% Source word
x = randn(N,1);

% Channel: perturb T randomly chosen entries
q = randperm(m_st-1);
c = zeros(m_st,1);
c(q(1:T)) = randn(T,1);
y = A*x-c;
err_loc = (rand(M-N,1)>=.9);
err_new = err_loc.*randn(M-N,1);
w_vec = B_stack*x-err_new;


% initial estimate of x
% l1_min or direct solution (LS) depending on n_st
iA = inv(A);
x0 = iA*y;
e0 = A*x0-y;
lambda_0 = zeros(m_st,1);

A_old = A;
y_old = y;
x_old = x0;
c_old = e0;
lambda_old = lambda_0;
gamma_old = []; %find(A*x0-y) or find(abs(lambda_0)==1);
m_old = m_st;
all_done = 0;
cond_list = [];
cond_iter = 1;
gamma_kc_new = [1:N]';
iGg_new = iA';

% number of new measurements added at a time.
m_u = 8; 

total_iter = 0;

for iter_mat = 1:m_u:M-m_st
    m_new = m_old + m_u;
    x_k = x_old;
    lambda_k = lambda_old;
    gamma_k = gamma_old;
    %     gamma_kc = gamma_kc_old;

    B = B_stack(iter_mat:iter_mat+m_u-1,:);
    G = [A_old' B'];
    w = w_vec(iter_mat:iter_mat+m_u-1);
    q = [y_old; w];
    
    c_original = [c; err_new(1:iter_mat+m_u-1)];
    c_old = [c_old; B*x_k-w];

    if abs(abs(B*x_k-w))>=10*eps
        nu = sign(B*x_k-w);
        i_nu = [m_new-m_u+1:m_new]';
        gammak_new = [gamma_k; i_nu];
    else
        stp = 1; % Something is wrong! (degenracy)
        break;
    end

    lambda_old = [lambda_k; nu]; % xi_old

    done = 0;
    epsilon_old = 0;
    iter_in = 0;
    while ~done
        iter_in = iter_in+1;
        %[iter_mat iter_in]
    
        gamma_k = gammak_new;
        x_k = x_old;
        c_k = c_old;
        lambda_k = lambda_old; % xi_k

        gamma_kc = gamma_kc_new;

        % INVERSE UPDATE
        %iGg = inv(G(:,gamma_kc));
        iGg = iGg_new;
        del_lambda = -iGg*(G(:,i_nu)*sign(lambda_k(i_nu))); % del_xi

        del_lambda_vec = zeros(m_new,1); % del_xi_vec
        del_lambda_vec(gamma_kc) = del_lambda;

        % find incoming elements in support of error vector
        constr1 = (1-lambda_k(gamma_kc))./del_lambda_vec(gamma_kc);
        constr1_pos = constr1(constr1>2*eps);
        delta1 = min(constr1_pos);
        idelta1 = gamma_kc(find(constr1==delta1));

        constr2 = -(1+lambda_k(gamma_kc))./del_lambda_vec(gamma_kc);
        constr2_pos = constr2(constr2>2*eps);
        delta2 = min(constr2_pos);
        idelta2 = gamma_kc(find(constr2==delta2));

        if delta1>delta2
            delta = delta2;
            idelta = idelta2;
        else
            delta = delta1;
            idelta = idelta1;
        end

        if epsilon_old+delta <=1
            lambdak_1 = lambda_k + delta*del_lambda_vec;
            lambdak_1(idelta) = sign(lambdak_1(idelta));
            gammak_temp = [gamma_k; idelta];

            % outgoing element from the support of
            uz_vec = zeros(N,1);
            idelta_gamma_kc = find(idelta==gamma_kc);
            uz_vec(idelta_gamma_kc) = sign(lambdak_1(idelta));
            delx = iGg'*uz_vec;

            delc= G'*delx;
            constr3 = -(c_k(gamma_k)./delc(gamma_k));
            constr3_pos = constr3(constr3>2*eps);
            delta3 = min(constr3_pos);
            idelta3 = gamma_k(find(constr3==delta3));
            xk_1 = x_k + delta3*delx;
            ck_1 = c_k + delta3*delc;
            ck_1(idelta3) = 0;

            out_lambda = idelta3;

            % See if signal is decoded exactly...
            if [m_new-N~=length(find(abs(ck_1)>=1e-10))]
                gamma_final = setdiff(gammak_temp,find(abs(G'*xk_1-q)<=100*eps));
                x_old = xk_1;
                c_old = ck_1;
                lambda_old = lambdak_1;
                stp = 1; %% WHAT TO DO NOW ???
                all_done = 1;
                break; break; break;
            end
            gammak_new = setdiff(gammak_temp,idelta3);
            gamma_kc_new = gamma_kc;
            gamma_kc_new(idelta_gamma_kc) = idelta3;

            % UPDATE INVERSE OF G(gamma_kc)
            % Gg = Gg + (g_new - g_old)*1_idelta_gamma_kc';
            diff_G_in_out = G(:,idelta3)-G(:,idelta);
            iGg_new = iGg - (iGg*diff_G_in_out/(1+iGg(idelta_gamma_kc,:)*diff_G_in_out))*(iGg(idelta_gamma_kc,:));

            x_old = xk_1;
            c_old = ck_1;
            lambda_old = lambdak_1;
            gamma_old = gammak_new;

            if ~isempty(find(idelta3 == i_nu, 1))
                lambda_old(idelta3) = (epsilon_old+delta)*lambdak_1(idelta3);
                i_nu = setdiff(i_nu,idelta3);
                if isempty(i_nu)
                    epsilon = 1;
                else
                    epsilon = delta + epsilon_old;
                end
            else
                epsilon = delta + epsilon_old;
            end

            if epsilon >=1
                done = 1;

                A_old = G';
                y_old = q;
                m_old = m_new;
            end
            epsilon_old = epsilon;

        else
            new_delta = 1-epsilon_old;
            lambdak_1 = lambda_k + new_delta*del_lambda_vec;
            lambda_old = [lambdak_1];
            gamma_old = gammak_new;

            done = 1;

            A_old = G';
            y_old = q;
            m_old = m_new;
        end
    end
    figure(2);  clf;
    subplot(211);plot(x_old); hold on; plot(x,'r.');
    title(['Estimated signal with ', num2str(iter_mat+m_u-1+m_st), ' measurements in the presence of ',num2str(T),  '  errors']);
    legend('Estimation', 'Original');
    subplot(212);plot(c_old); hold on; plot(c_original,'r.');
    title(['Estimated and original sparse error vectors']);
    legend('Estimation', 'Original'); drawnow

    total_iter = total_iter+iter_in;
    disp(['No. of errors = ',num2str(nnz(c_original)), ', Current no. of measurements = ',num2str(iter_mat+m_u-1+m_st), ', Homotopy iterations = ', num2str(iter_in), ',  Running sum of iterations = ', num2str(total_iter)]);
            
    if all_done == 1
        % Need L1-magic package for this verification.
        % x02 = G'\q; %inv(G*G')*G*q;
        % [xp2 lamf1 lamf2] = l1decode_pd2(x02, G', [], q, 1e-12, 30);
        break;
    end
end
x_out = x_old;
figure(2);  clf;
subplot(211);plot(x_old); hold on; plot(x,'r.');
title(['Estimated signal with ', num2str(iter_mat+m_u-1+m_st), ' measurements in the presence of ',num2str(T),  '  errors']);
legend('Estimation', 'Original');
subplot(212);plot(c_old); hold on; plot(c_original,'r.');
title(['Estimated and original sparse error vectors']);
legend('Estimation', 'Original'); drawnow
