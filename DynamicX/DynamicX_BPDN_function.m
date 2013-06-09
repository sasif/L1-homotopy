% DynamicX_BPDN_function.m
% 
% Solves  
% min_x  \tau ||x||_1 + 1/2*||Ax-(1-\epsilon)y-(\epsilon)yt||_2^2
%
% using homotopy update method while changing epsilon from 0 to 1.
%
% [xp_h, gamma_xh, iter, cpu] = DynamicX_BPDN_function(A, AtAgx, iAtAgx, y, yt, xp, gamma_x, pk, tau, maxiter)
%
% A: mxn matrix
% AtAgx: AtA submatrix at row and column indices at gamma_x.
% iAtAgx: inv(AtAgx)
% y: old set of measurements
% yt: new measurements
% xp: old solution
% gamma_x: old support
% pk: old constraints
% tau: threshold parameter. 
% cpu: CPU time

% LASSO update with dynamic change in x
% Author: Muhammad Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
% Created: July 2008
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

function [xp_h, gamma_xh, iter,cpu] = DynamicX_BPDN_function(A, AtAgx, iAtAgx, y, yt, xp, gamma_x, pk, tau, maxiter)
t0 = cputime;
N = length(xp);
xp_h = xp;
gamma_xh = gamma_x;
pk_old = pk;

iter = 0;
done = 0;

% if ~isa(A, 'function_handle')
%   Atf = @(x) A'*x;
%   Af = @(x) A*x;
%   AtAf = @(x)Atf(Af(x));
% end
warning('off','MATLAB:divideByZero');

e0 = 0;
itheta = [];
d = A'*(y-yt);

while iter < maxiter
    iter = iter +1 ;

    delx = -iAtAgx*d(gamma_xh);
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;

    pk = pk_old;
    %dk = A'*(A*delx_vec+y-yt); 
    % Agx = A(:,gamma_xh);
    Ag_delx = A(:,gamma_xh)*delx;
    AtAg_delx = A'*Ag_delx;
    dk = AtAg_delx+d;
    %dk = AtAf(delx_vec)+d;
    
    temp_gamma = zeros(N,1);
    temp_gamma(gamma_xh) = gamma_xh;
    gamma_xc = find([1:N]' ~= temp_gamma);
    % gamma_xc = setdiff([1:N],[gamma_xh]);

    b_constr1 = (tau-pk(gamma_xc))./dk(gamma_xc);
    b_constr2 = (tau+pk(gamma_xc))./-dk(gamma_xc);
    b_constr3 = (-xp_h(gamma_xh)./delx_vec(gamma_xh));
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

    epsilon = (theta+e0);
    if epsilon < 0
        epsilon= inf;
    end
    
    if epsilon >=1
        theta_end = (1-e0);
        xp_h = xp_h + theta_end*delx_vec;
        cpu = cputime-t0;
        break;
    end
    
    gamma_xh_old = gamma_xh;
    xp_old = xp_h;
    xp_h = xp_h+theta*delx_vec;
    
    if flag == 1
        AtgxAnx = A(:,gamma_xh)'*A(:,itheta);
        AtAgx_mod = [AtAgx AtgxAnx; AtgxAnx' A(:,itheta)'*A(:,itheta)];

        iAtAgx = update_inverse(AtAgx_mod, iAtAgx,1);
        AtAgx = AtAgx_mod;
        
        AtAgx_mod = [];
        gamma_xh = [gamma_xh; itheta];
        xp_h(itheta) = 0;
    else
        g_old = gamma_xh;
        outx_index = find(gamma_xh==itheta);
        len_gamma = length(gamma_xh);
        gamma_xh(outx_index) = gamma_xh(len_gamma);
        gamma_xh = gamma_xh(1:len_gamma-1);
        
        rowi = outx_index; % ith row of A is swapped with last row (out_x)
        colj = outx_index; % jth column of A is swapped with last column (out_lambda)
        
        AtAgx_ij = AtAgx;
        temp_row = AtAgx_ij(rowi,:);
        AtAgx_ij(rowi,:) = AtAgx_ij(len_gamma,:);
        AtAgx_ij(len_gamma,:) = temp_row;
        temp_col = AtAgx_ij(:,colj);
        AtAgx_ij(:,colj) = AtAgx_ij(:,len_gamma);
        AtAgx_ij(:,len_gamma) = temp_col;
        
        iAtAgx_ij = iAtAgx;
        temp_row = iAtAgx_ij(colj,:);
        iAtAgx_ij(colj,:) = iAtAgx_ij(len_gamma,:);
        iAtAgx_ij(len_gamma,:) = temp_row;
        temp_col = iAtAgx_ij(:,rowi);
        iAtAgx_ij(:,rowi) = iAtAgx_ij(:,len_gamma);
        iAtAgx_ij(:,len_gamma) = temp_col;
        
        AtAgx = AtAgx_ij(1:len_gamma-1,1:len_gamma-1);
        iAtAgx = update_inverse(AtAgx_ij, iAtAgx_ij,2);
        
        xp_h(itheta) = 0;
    end

    e0 = (theta+e0);

    pk_old = pk+theta*dk;
    pk_old([gamma_xh; itheta]) = sign(pk_old([gamma_xh; itheta]))*tau;
    %pk_old(gamma_xh) = sign(pk_old(gamma_xh))*tau;
end
cpu = cputime-t0;

% 
% function iAtB = update_inverse(AtB, iAtB_old,flag);
% 
% n = size(AtB,1);
% 
% %A11 = AtB(1:n-1,1:n-1);
% A12 = AtB(1:n-1,n);
% A21 = AtB(n,1:n-1);
% A22 = AtB(n,n);
% 
% % add columns
% if flag == 1
%     iA11 = iAtB_old;
%     iA11A12 = iA11*A12;
%     A21iA11 = A21*iA11;
%     S = A22-A21*iA11A12;
%     Q11_right = iA11A12*(A21iA11/S); 
% %     Q11 = iA11+ Q11_right;
% %     Q12 = -iA11A12/S;
% %     Q21 = -A21iA11/S;
% %     Q22 = 1/S;
% 
%     iAtB = zeros(n);
%     %iAtB = [Q11 Q12; Q21 Q22]; 
%     iAtB(1:n-1,1:n-1) = iA11+ Q11_right;
%     iAtB(1:n-1,n) = -iA11A12/S; 
%     iAtB(n,1:n-1) =  -A21iA11/S;
%     iAtB(n,n) = 1/S;
% %delete columns
% else if flag == 2
%         Q11 = iAtB_old(1:n-1,1:n-1);
%         Q12 = iAtB_old(1:n-1,n);
%         Q21 = iAtB_old(n,1:n-1);
%         Q22 = iAtB_old(n,n);
%         Q12Q21_Q22 = Q12*(Q21/Q22);
%         iAtB = Q11 - Q12Q21_Q22;
%         %iAtB = iA11;
%     end
% end