% BPDN_homotopy_Weighted
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \e_i |x_i| + 1/2*||y-Ax||_2^2
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
%
%-------------------------------------------+
% Copyright (c) 2010.  Muhammad Salman Asif 
%-------------------------------------------+

clear
clc

% load RandomStates
rseed = sum(100*clock);
rseed = 0;
% rand('state',rseed);
% randn('state',rseed);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));

N = 128;   % signal length
T = 20;    % sparsity level
M = 100;    % no. of measurements

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
sigma = 0.01;
e = randn(M,1)*sigma;
y = A*x+e;


% Regularization parameters 
epsilon_vec = [ones(N/2,1); ones(N/2,1)*1e4]*1e-2;

tau = 1e-2;
W_new = ones(N,1);
W_new(abs(x)>0) = 1e-5./abs(x(abs(x)>0));

epsilon_vec = tau.*W_new;
unique_eps = sort(unique(epsilon_vec),'descend');
thresh = min(epsilon_vec);

maxiter = 4*N;

in = [];
in.x_orig = x;
in.epsilon_vec = epsilon_vec;
in.maxiter = maxiter;
in.delx_mode = 'mil';
out = wtBPDN_function(A,y,in);
xh_wt = out.x_out;
gamma_wt = out.gamma;
iter_wt = out.iter;

%%

% Initialization of primal sign and support
z_x = zeros(N,1);
gamma_x = [];       % Primal support

% Initial step
Primal_constrk = -A'*y; 
constr_mask = abs(Primal_constrk)>epsilon_vec;
[c i] = max(abs(Primal_constrk.*constr_mask));
eps_iter = sum(unique_eps>c)+1;

gamma_xk = i;

z_x(gamma_xk) = -sign(Primal_constrk(gamma_xk));
epsilon = c;
Primal_constrk(gamma_xk) = sign(Primal_constrk(gamma_xk))*epsilon;
xk_1 = zeros(N,1);

z_xk = z_x;

% loop parameters
done = 0;
iter = 0;
data_precision = eps;   % floating point precision

old_delta = 0;
out_x = [];
count_delta_stop = 0;

% Turn the plots on and off
constraint_plots = 1; % 1 - draw the plots after every homotopy iteration.

AtAgx = A(:,gamma_xk)'*A(:,gamma_xk);
iAtAgx = inv(AtAgx);

while ~done
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')

    gamma_x = gamma_xk;
    z_x = z_xk;
    x_k = xk_1;

    %%%%%%%%%%%%%%%%%%%%%
    %%%% update on x %%%%
    %%%%%%%%%%%%%%%%%%%%%

    % Update direction
%     del_x = -inv(A(:,gamma_lambda)'*A(:,gamma_x))*z_lambda(gamma_lambda);
%     diff_inv = [ max(max(abs(inv(A(:,gamma_lambda)'*A(:,gamma_x))-iAtglAgx)))]
%     if diff_inv > 1e-8
%         stp =1;
%     end
    indicator_temp = epsilon>epsilon_vec(gamma_x);
    del_x = iAtAgx*(indicator_temp.*z_x(gamma_x));
    del_x_vec = zeros(N,1);
    del_x_vec(gamma_x) = del_x;

    pk = Primal_constrk;
    % dk = AtA*del_x_vec;
    dk = A'*(A(:,gamma_x)*del_x);
    
    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW. 
    pk_temp = Primal_constrk;
    gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,1e-12));
%     pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
    
    xk_temp = x_k;
    gammaX_temp = find(abs(x_k)<1*eps);
    %%%---
    
    temp_gamma = zeros(N,1);
    temp_gamma(gamma_x) = gamma_x;
    gamma_xc = find([1:N]' ~= temp_gamma);
    
    % For different values of regularization parameters
    eps_temp = epsilon.*(epsilon>epsilon_vec(gamma_xc))+epsilon_vec(gamma_xc).*(epsilon<=epsilon_vec(gamma_xc));
    one_temp = epsilon>epsilon_vec(gamma_xc);
    if find(one_temp-1)
        stp = 1;
    end
    if find(eps_temp-epsilon)
        stp = 1;
    end
    delta1_constr = (eps_temp-pk(gamma_xc))./(one_temp+dk(gamma_xc));
    delta1_pos_ind = find(delta1_constr>0);
    delta1_pos = delta1_constr(delta1_pos_ind);
    [delta1 i_delta1] = min(delta1_pos);
    if isempty(delta1)
        delta1 = inf;
    end
    delta2_constr = (eps_temp+pk(gamma_xc))./(one_temp-dk(gamma_xc));
    delta2_pos_ind = find(delta2_constr>0);
    delta2_pos = delta2_constr(delta2_pos_ind);
    [delta2 i_delta2] = min(delta2_pos);
    if isempty(delta2)
        delta2 = inf;
    end
    
    if delta1>delta2
        delta = delta2;
        i_delta = gamma_xc(delta2_pos_ind(i_delta2));
    else
        delta = delta1;
        i_delta = gamma_xc(delta1_pos_ind(i_delta1));
    end
    
    delta3_constr = (-x_k(gamma_x)./del_x_vec(gamma_x));
    delta3_pos_index = find(delta3_constr>0);
    [delta3 i_delta3] = min(delta3_constr(delta3_pos_index));
    out_x_index = gamma_x(delta3_pos_index(i_delta3));
    
    chk_x = 0;
    out_x = [];
    if delta3 > 0 & delta3 <= delta
        chk_x = 1;
        delta = delta3;
        out_x = out_x_index;
    end

    [i_delta out_x delta -chk_x]
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
    
    if epsilon <= unique_eps(eps_iter)
        epsilon = unique_eps(eps_iter);
        delta_end = epsilon_old-epsilon;
        Primal_constrk = pk+delta_end*dk;
        epsilon_temp = epsilon.*(epsilon>epsilon_vec(gamma_xk))+epsilon_vec(gamma_xk).*(epsilon<=epsilon_vec(gamma_xk));
        Primal_constrk(gamma_x) = sign(Primal_constrk(gamma_x)).*epsilon_temp;
 
        xk_1 = x_k + delta_end*del_x_vec;
        eps_iter = eps_iter+1;
        if eps_iter > length(unique_eps)
            disp('done!');
            break;
        else
            disp('switch epsilon!');
            continue; 
        end
    end
    if epsilon <= min(epsilon_vec); %sqrt(2*log(N))*sigma; %1e-7 %|| iter > 5*T || (length(gamma_lambda) == K)
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
        gx_old = gamma_x;

        outx_index = find(gamma_x==out_x);
        gamma_x(outx_index) = gamma_x(end);
        gamma_x(end) = out_x;
        gamma_xk = gamma_x(1:end-1);

        rowi = outx_index; % ith row of A is swapped with last row (out_x)
        colj = outx_index; % jth column of A is swapped with last column (out_lambda)
        AtAgx_ij = AtAgx;
        temp_row = AtAgx_ij(rowi,:);
        AtAgx_ij(rowi,:) = AtAgx_ij(end,:);
        AtAgx_ij(end,:) = temp_row;
        temp_col = AtAgx_ij(:,colj);
        AtAgx_ij(:,colj) = AtAgx_ij(:,end);
        AtAgx_ij(:,end) = temp_col;
        
        iAtAgx_ij = iAtAgx;
        temp_row = iAtAgx_ij(colj,:);
        iAtAgx_ij(colj,:) = iAtAgx_ij(end,:);
        iAtAgx_ij(end,:) = temp_row;
        temp_col = iAtAgx_ij(:,rowi);
        iAtAgx_ij(:,rowi) = iAtAgx_ij(:,end);
        iAtAgx_ij(:,end) = temp_col;

        AtAgx = AtAgx_ij(1:end-1,1:end-1);
        iAtAgx = update_inverse(AtAgx_ij, iAtAgx_ij,2);
        xk_1(out_x) = 0;
    else
        % If an element is added to gamma_x
        gamma_xk = [gamma_x; i_delta];
        i_theta = i_delta;
        new_x = i_delta;
        AtAnl = A(:,gamma_x)'*A(:,new_x);
        AtAgx_mod = [AtAgx A(:,gamma_x)'*A(:,i_theta); AtAnl' A(:,new_x)'*A(:,i_theta)];

        AtAgx = AtAgx_mod;
        iAtAgx = update_inverse(AtAgx, iAtAgx,1);
        xk_1(i_theta) = 0;
        
        gamma_x  = gamma_xk;
    end

    z_xk = zeros(N,1);
    z_xk(gamma_xk) = -sign(Primal_constrk(gamma_xk));
    epsilon_temp = epsilon.*(epsilon>epsilon_vec(gamma_x))+epsilon_vec(gamma_x).*(epsilon<=epsilon_vec(gamma_x));
    Primal_constrk(gamma_x) = sign(Primal_constrk(gamma_x)).*epsilon_temp;
    
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
            plot(new_x, Primal_constrk(new_x),'or','MarkerSize',18,'LineWidth',2);
            text(new_x, Primal_constrk(new_x)*1.1, ['Incoming \gamma = ',num2str(new_x)],'FontSize',14);
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

        figure(1);
        subplot(2,1,2)
        hold off
        plot(x_k,'.r','MarkerSize',14); hold on;
        plot(xk_1,'LineWidth',1);
        if chk_x == 1
            plot(out_x, 0,'ok', 'MarkerSize',18,'LineWidth',2);
        else
            plot(new_x, 0,'or', 'MarkerSize',18,'LineWidth',2);
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
            %pause
        end

        %drawnow;
        %print(gcf,'-dbmp','-r200',['BPDNPath_',num2str(iter)])
    end
%     pause
%     [iter epsilon delta]
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
plot(x_k,'.r','MarkerSize',14); hold on;
plot(xk_1,'LineWidth',1);
set(gca,'FontSize',16,'XLim',[1 N]);
title(['BPDN estimate at \tau = ',num2str(thresh), ', iter = ', num2str(iter)]);

%%
digit = 6; if sigma > 0; digit = 7; end
opts = [];
opts.tol = 10^(-digit);
W_new = ones(N,1);
opts.weights = W_new;
opts.print = 0;
opts.maxit = maxiter;

opts.nu = 0; opts.rho = tau;
tic;
[x_yall1,Out] = yall1(A,y,opts);

%% CVX check
cvx_begin
    cvx_precision best
    variable xc(N);
    minimize( sum(epsilon_vec.*abs(xc))+ 1/2*sum(square(A*xc-y)))
cvx_end

%%
figure(101); plot([xc xk_1])
figure(102); plot([xc-xk_1])
