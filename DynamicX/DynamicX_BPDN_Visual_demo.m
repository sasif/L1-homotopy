% Test file for LASSO (BPDN) update with dynamic change in x

% Author: Muhammad Salman Asif @ Georgia Tech
% E-mail: sasif@ece.gatech.edu
% Created: November 2008
% Modified: February 2009
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+


clear

% % load fixed random states
% load RandomStates
% rand('state', rand_state);
% randn('state', randn_state);

% signal length
N = 128;
% number of spikes to put down
T = 20;
% number of observations to make
M = 64;

% random signal supported on q(1:T)
x = zeros(N,1);
q = randperm(N);

%N_S = 1: Gaussian signal, N_S = 0: Random signed sequence
N_S = 0;
if N_S == 1
    x(q(1:T)) = randn(T,1);
else
    x(q(1:T)) = sign(randn(T,1));
end

% measurement matrix
%disp('Creating measurment matrix...');
G = randn(M,N)/sqrt(M);
A = (orth(G'))';
% A = G;

% observations
sigma = .01;
e = randn(M,1)*sigma;
y = A*x + e;

% regularization parameter
tau = 1e-2;

% Solve for measurements y
t0 = clock;
[xp, gamma_x, xp_iter, xp_time] = BPDN_homotopy_function(A, y, tau,4*N); %BPDN
t1 = clock;
t0_1 = etime(t1,t0);

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

epsilon = 0;
e0 = 0;
xp2_h = xp;
gamma_xh = gamma_x;
itheta = [];

done = 0;
iter = 0;
pk_old = A'*(A*xp-y);
AtAgx = A(:,gamma_xh)'*A(:,gamma_xh);
iAtAgx = inv(AtAgx);

tu0 = clock;
while ~done
    iter = iter +1 ;

    %%% NEED UPDATE
    % iAtAgx = inv(A(:,gamma_xh)'*A(:,gamma_xh));

    d = A(:,gamma_xh)'*(-y+yt);
    delx = iAtAgx*d;
    pk = pk_old;
    dk = A'*(A(:,gamma_xh)*delx+y-yt);
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;

    temp_gamma = zeros(N,1);
    temp_gamma(gamma_xh) = gamma_xh;
    gamma_xc = find([1:N]' ~= temp_gamma);
    % gamma_xc = setdiff([1:N],[gamma_xh]);

    b_constr1 = (tau-pk(gamma_xc))./dk(gamma_xc);
    b_constr2 = (tau+pk(gamma_xc))./-dk(gamma_xc);
    b_constr3 = (-xp2_h(gamma_xh)./delx_vec(gamma_xh));
    itheta_1 = find(b_constr1>2*eps);
    itheta_2 = find(b_constr2>2*eps);
    itheta_3 = find(b_constr3>2*eps);
    theta1 = min(b_constr1(itheta_1));
    theta2 = min(b_constr2(itheta_2));
    theta3 = min(b_constr3(itheta_3));
    if theta1>theta2
        theta = theta2;
        itheta = gamma_xc(find(b_constr2==theta2));
        flag = 1;
    else
        theta = theta1;
        itheta = gamma_xc(find(b_constr1==theta1));
        flag = 1;
    end
    if theta3 < theta
        theta = theta3;
        itheta = gamma_xh(find(b_constr3==theta3));
        flag = 0;

    end
    %     if flag == 0
    %         [epsilon+theta theta itheta flag xp2_h(itheta) delx_vec(itheta)]
    %     else
    %         [epsilon+theta theta itheta flag pk(itheta) dk(itheta)]
    %     end
    
    epsilon = (theta+e0);
    if epsilon < 0
        epsilon= inf;
    end
    gamma_xh_old = gamma_xh;
    xp2_old = xp2_h;
    xp2_h = xp2_h+theta*delx_vec;
    if flag == 1
        AtgxAnx = A(:,gamma_xh)'*A(:,itheta);
        AtAgx_mod = [AtAgx AtgxAnx; AtgxAnx' A(:,itheta)'*A(:,itheta)];

        iAtAgx = update_inverse(AtAgx_mod, iAtAgx,1);
        AtAgx = AtAgx_mod;

        gamma_xh = [gamma_xh; itheta];
        xp2_h(itheta) = 0;
    else
        outx_index = find(gamma_xh==itheta);
        gamma_xh(outx_index) = gamma_xh(end);
        gamma_xh(end) = itheta;
        gamma_xh = gamma_xh(1:end-1);

        rowi = outx_index; % ith row of A is swapped with last row (out_x)
        colj = outx_index; % jth column of A is swapped with last column (out_x)
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

        xp2_h(itheta) = 0;
    end

    %[epsilon theta itheta flag]
    %     itheta_stack = [itheta_stack; itheta];
    %     flag_stack = [flag_stack; flag];

    if epsilon >=1
        theta_end = (1-e0);
        xp2_h = xp2_old + theta_end*delx_vec;
        gamma_xh = gamma_xh_old;
        epsilon = 1;
        disp('done!');
        break;
    end

    e0 = (theta+e0);

    pk_old = pk+theta*dk;
    pk_old([gamma_xh; itheta]) = sign(pk_old([gamma_xh; itheta]))*tau;

    figure(2);
    subplot(2,1,1);
    hold off;
    plot(pk,'.r','MarkerSize',14);
    hold on;
    plot(pk_old, 'LineWidth',1);

    if flag == 1
        plot(itheta, pk_old(itheta),'or','MarkerSize',18, 'LineWidth',2);
        text(itheta,pk_old(itheta)*1.2, ['Incoming \gamma = ',num2str(itheta)],'FontSize',14);
    else
        plot(itheta, pk_old(itheta),'*k','MarkerSize',18, 'LineWidth',2);
        text(itheta,pk_old(itheta)*1.2, ['Outgoing \gamma = ',num2str(itheta)],'FontSize',14);
    end
    set(gca,'FontSize',16, 'XLim',[1 N] );
    title('BPDN shrinkage constraints');
    %legend('Old constraint', 'New constraint', 'Element change','Location','NorthEast');
    title(['BPDN constraints, n = ',num2str(N), ', m = ', num2str(M), ', T = ', num2str(T), ', Tin = ', num2str(T_in)]);
    plot(1:N, tau*ones(1,N),'--k');
    plot(1:N, -tau*ones(1,N), '--k');
    axis([1 N -tau*1.5 tau*1.5])

    figure(2);
    subplot(2,1,2)
    hold off
    plot(xp2_old,'.r','MarkerSize',14); hold on;
    plot(xp2_h,'LineWidth',1);
    set(gca,'FontSize',16,'XLim',[1 N]);
    title({['BPDN estimate at \epsilon = ',num2str(epsilon),', Iteration count = ',num2str(iter)]});
    if flag == 0
        plot(itheta, 0,'ok', 'MarkerSize',18,'LineWidth',2);
    else
        plot(itheta, 0,'or', 'MarkerSize',18,'LineWidth',2);
    end
    drawnow;
    % pause;

    %print(gcf,'-dbmp','-r200',['BPDNPath_',num2str(iter)])
end
figure(2);
subplot(2,1,1);
hold off;
plot(pk,'.r','MarkerSize',14);
hold on;
plot(pk_old, 'LineWidth',1);


set(gca,'FontSize',16, 'XLim',[1 N] );
title('BPDN shrinkage constraints');
%legend('Old constraint', 'New constraint', 'Element change','Location','NorthEast');
title(['BPDN constraints, n = ',num2str(N), ', m = ', num2str(M), ', T = ', num2str(T), ', Tin = ', num2str(T_in)]);
plot(1:N, tau*ones(1,N),'--k');
plot(1:N, -tau*ones(1,N), '--k');
axis([1 N -tau*1.5 tau*1.5])

figure(2);
subplot(2,1,2)
hold off
plot(xp2_old,'.r','MarkerSize',14); hold on;
plot(xp2_h,'LineWidth',1);
set(gca,'FontSize',16,'XLim',[1 N]);
title({['BPDN estimate at \epsilon = ',num2str(epsilon),', Iteration count = ',num2str(iter)]});
drawnow;
