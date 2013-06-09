% Test out l1 decoding with homotopy update scheme
%
% Created by Salman Asif @ Georgia Tech

close all
clear

randn('state',0);
rand('state',0);

% source length
N = 128;

% Number of errors introduced at random
T_list = [10 20 30 40 50];
P_list = [1 2 4 8];

Stack_meas = [];
Stack_iter = [];

no_exp = 20; % number of simulations

for P = P_list
    % Number of new elements added at a time
    No_new_element = P;
    
    table_iter = zeros(no_exp, length(T_list));
    table_error = zeros(no_exp, length(T_list));
    
    fprintf('N = %d, P = %d, sim_count = %d...\n',N,No_new_element,no_exp);
    
    for iT = 1:length(T_list)
        for inn_iter = 1:no_exp
            
            % Maximum allowed length of message
            M = 20*N;
            
            % number of perturbations
            T = T_list(iT); %round(.3*N);
            
            % coding matrix
            Orth_mat = randn(M,M);
            % Orth_mat = randsrc(M,M);
            % [Orth_mat R_mat] = qr(Orth_mat,0);
            G_mat = Orth_mat(:,1:N);
            
            % Homotopy setup
            m_st = N;   % number of measurements to start with
            A = G_mat(1:m_st,:);
            B_stack = G_mat(m_st+1:end,:);
            
            % source word
            x = randn(N,1);
            
            % channel: perturb T randomly chosen entries
            q = randperm(m_st-1);
            c = zeros(m_st,1);
            c(q(1:T)) = randn(T,1);
            
            y = A*x-c;
            
            % Introduce sparse errors
            err_loc = (rand(M-N,1)>1);
            err_new = err_loc.*randn(M-N,1);
            w_vec = B_stack*x+err_new;
            
            % initial estimate of x
            % l1_min or direct inversion depending on n_st
            iA = inv(A);
            x0 = iA*y;
            e0 = A*x0-y;
            lambda_0 = zeros(m_st,1);
            
            A_old = A;
            y_old = y;
            x_old = x0;
            c_old = e0;
            lambda_old = lambda_0;
            
            m_u = No_new_element;
            gamma_old = []; %find(A*x0-y) or find(abs(lambda_0)==1);
            m_old = m_st;
            all_done = 0;
            cond_list = [];
            cond_iter = 1;
            gamma_kc_new = [1:N]';
            iGg_new = iA';
            
            iter_total = 0;
            iter_in = 0;
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
                
                c_old = [c_old; B*x_k-w];
                
                if abs(abs(B*x_k-w))>=10*eps
                    nu = sign(B*x_k-w);
                    i_nu = [m_new-m_u+1:m_new]';
                    gammak_new = [gamma_k; i_nu];
                else
                    stp = 1;
                end
                
                lambda_old = [lambda_k; nu]; % xi_old
                
                done = 0;
                epsilon_old = 0;
                iter_in = 0;
                while ~done
                    iter_in = iter_in+1;
                    % [iter_mat iter_in]
                    
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
                        
                        if [m_new-N~=length(find(abs(G'*xk_1-q)>=1e5*eps))]
                            gamma_final = setdiff(gammak_temp,find(abs(G'*xk_1-q)<=100*eps));
                            x_old = xk_1;
                            c_old = ck_1;
                            lambda_old = lambdak_1;
                            stp = 1; %% WHAT TO DO NOW ???
                            all_done = 1;
                            iter_total = iter_total+iter_in; 
                            break; break; break;
                        end 
                        gammak_new = setdiff(gammak_temp,idelta3);
                        gamma_kc_new = gamma_kc;
                        gamma_kc_new(idelta_gamma_kc) = idelta3;
                        
                        % UPDATE INVERSE OF G(gamma_kc)
                        % Gg = Gg + (g_new - g_old)'*1_idelta_gamma_kc;
                        diff_G_in_out = G(:,idelta3)-G(:,idelta);
                        iGg_new = iGg - (iGg*diff_G_in_out/(1+iGg(idelta_gamma_kc,:)*diff_G_in_out))*(iGg(idelta_gamma_kc,:));
                        
                        x_old = xk_1;
                        c_old = ck_1;
                        lambda_old = lambdak_1;
                        gamma_old = gammak_new;
                        
                        if length(find(idelta3 == i_nu))
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
                            iter_total = iter_total+iter_in;
                        end
                        epsilon_old = epsilon;
                        
                    else
                        new_delta = 1-epsilon_old;
                        lambdak_1 = lambda_k + new_delta*del_lambda_vec;
                        lambda_old = [lambdak_1];
                        gamma_old = gammak_new;
                        
                        iter_total = iter_total+iter_in;
                        done = 1;
                        
                        A_old = G';
                        y_old = q;
                        m_old = m_new;
                    end
                end
                % Current
                % disp(['No. of errors = ',num2str(T_list(iT)), ', Simulation run = ', num2str(inn_iter),', Homotopy iterations = ', num2str(iter_total), ', No. of measurements = ',num2str(iter_mat+m_u-1+m_st),', Error in recovery = ',num2str(norm(x-xk_1))]);
                
                if all_done == 1
                    break;
                end
            end
            % Total number of homotopy iterations.
            table_iterT(inn_iter,iT) = iter_total;
            % Total number of measurements required for perfect recovery!
            table_meas(inn_iter,iT) = iter_mat+m_u-1+m_st;
        end
    end
    % With few measurements smaller number of steps taken to update per
    % measuremnt and as we approach final solution number of steps per
    % measurement increase slightly.
    disp('Final results');
    Error_count = [T_list] % Total number of sparse errors
    Total_measurements = mean(table_meas,1) % Total number of measurements required for perfect recovery!
    Total_iterations = mean(table_iterT,1)  % Total number of homotopy iterations.
    Average_iterations = Total_iterations./(Total_measurements-N) % Avg homotopy iterations per new measurement
    
    Stack_meas = [Stack_meas; Total_measurements];
    Stack_iter = [Stack_iter; Total_iterations];    
end

%%
plot_results = 1;
if plot_results
exp_name = 'l1decoding_N64_Exp20';
exp_name = 'l1decoding_N128_Exp20';

load(exp_name);

addpath ../
addpath ../utils/utils_fig 

% figure setup
marker = {'bx','r+','k*','gd'};
linewidth = 1.5;


axis_prop = {};
axis_prop{1,4} = 'FontSize'; axis_prop{2,4} = 14;
axis_prop{1,5} = 'FontWeight'; axis_prop{2,5} = 'normal';
axis_prop{1,8} = 'XGrid'; axis_prop{2,8} = 'on';
axis_prop{1,9} = 'GridLineStyle'; axis_prop{2,9} = '--';

ydata = Stack_iter./(Stack_meas-N); 
figure(1)
clf; hold on;
set(gca,'FontSize',14,'FontWeight','normal');
plot(T_list, ydata(1,:),'-bx', 'LineWidth',1.5,'MarkerSize',10);
plot(T_list, ydata(2,:),'-.r+', 'LineWidth',1.5,'MarkerSize',10);
plot(T_list, ydata(3,:),':k*', 'LineWidth',1.5,'MarkerSize',10);
plot(T_list, ydata(4,:),'--md', 'LineWidth',1.5,'MarkerSize',10);
xlabel('Number of errors');
ylabel('Number of homotopy iterations');
% title('Average number of homotopy per measurement');
legend('P=1','P=2','P=4','P=8','Location','NorthWest');
axis tight;

img_name = sprintf('%s-AvgIter',exp_name);
set(gcf, 'Color', 'w');
set(gcf,'Position',[300 0 500 600]);

eval(sprintf('export_fig %s.png',img_name));
eval(sprintf('export_fig %s.pdf',img_name));

ydata = Stack_iter; 
figure(2)
clf; hold on;
set(gca,'FontSize',14,'FontWeight','normal');
plot(T_list, ydata(1,:),'-bx', 'LineWidth',1.5,'MarkerSize',10);
plot(T_list, ydata(2,:),'-.r+', 'LineWidth',1.5,'MarkerSize',10);
plot(T_list, ydata(3,:),':k*', 'LineWidth',1.5,'MarkerSize',10);
plot(T_list, ydata(4,:),'--md', 'LineWidth',1.5,'MarkerSize',10);
xlabel('Number of errors');
ylabel('Number of homotopy iterations');
% title('Total number of homotopy iteration');
legend('P=1','P=2','P=4','P=8','Location','NorthWest');
axis tight;

img_name = sprintf('%s-TotalIter',exp_name);
set(gcf, 'Color', 'w');
set(gcf,'Position',[300 0 500 600]);

eval(sprintf('export_fig %s.png',img_name));
eval(sprintf('export_fig %s.pdf',img_name));
end