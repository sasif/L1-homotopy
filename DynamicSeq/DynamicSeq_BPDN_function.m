% DynamicLasso_update_function.m
% 
% Solves  
% min_x  \epsilon ||x||_1 + 1/2*||Ax-y||_2^2 + 1/2* ||bx-w||_2^2
%
% using homotopy update method
%
% [xp_h, gamma_xh, iter] = DynamicSeq_BPDN_function(A, b, AtAgx, iAtAgx, y, w, xp, gamma_x, pk, tau, chk_e0, maxiter)
%
% A: mxn matrix
% b: new rows in measurement matrix
% AtAgx: AtA submatrix at row and column indices at gamma_x.
% iAtAgx: inv(AtAgx)
% y: old set of measurements
% w: new measurements
% xp: old solution
% gamma_x: old support
% pk: old constraints
% tau: threshold parameter. 
% chk_e0: 0 or 1: This selects if we want to take U = (A'A+b'b) (1) or (A'A+e0*b'b) (0)...!
% In simulation U = (A'A+b'b) works better

% BPDN solution update with dynamic measurements 
% Author: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
% Created: February 2008
%
%-------------------------------------------+
% Copyright (c) 2008.  Muhammad Salman Asif 
%-------------------------------------------+

function [xp_h, gamma_xh, iter, th] = DynamicSeq_BPDN_function(A, b, AtAgx, iAtAgx, y, w, xp, gamma_x, pk, tau, chk_e0, maxiter)

t0 = cputime;
N = length(xp);
xp_h = xp;
gamma_xh = gamma_x;
itheta = [];
e0 = 0;

% Gram matrix update
% AtAgx = A(:,gamma_xh)'*A(:,gamma_xh)+chk_e0*b(:,gamma_xh)'*b(:,gamma_xh);
pk_old = pk;

notdone = 1;
iter = 0;

while iter < maxiter
    iter = iter +1 ;

    e0h = e0*(1-chk_e0)+1*chk_e0; 
    % this way we can set e0h to 1 or e0 which changes the definition of theta and U accordingly

    %%% UPDATE
    %     U = A(:,gamma_xh)'*A(:,gamma_xh)+e0h*b(:,gamma_xh)'*b(:,gamma_xh);
    %     iU = inv(U);
    U = AtAgx;
    iU = iAtAgx;

    d = b(:,gamma_xh)'*(b(:,gamma_xh)*xp_h(gamma_xh)-w);
    u = (b(:,gamma_xh)*iU*b(:,gamma_xh)');
    dx = -iU*d;
    pk = pk_old; %A'*(A*xp_h-y)+e0*b'*(b*xp_h-w);
    %dk = b'*(b*xp_h-w)+(Atf(A(:,gamma_xh)*dx)+b'*(e0h*(b(:,gamma_xh)*dx)));
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = dx;
    Adelx = A(:,gamma_xh)*dx;
    AtAdelx = (A')*Adelx;
    dk = b'*(b*xp_h-w)+AtAdelx+b'*(e0h*(b*delx_vec));

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
    if isempty(theta1)
        theta1 = inf;
    end
    if isempty(theta2)
        theta2 = inf;
    end
    if isempty(theta3)
        theta3 = inf;
    end

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
    epsilon = e0+(theta)/(1-theta*u);
    epsilon = (theta+e0-e0h*u*theta)/(1-theta*u);
    
    if epsilon < 0
        epsilon= inf;
    end

    if epsilon >=1
        theta_end = (1-e0)/(1+u-e0h*u);
        xp_h = xp_h + theta_end*delx_vec;
        th = cputime-t0;
        break;
    end

    xp_h_old = xp_h;
    gamma_xh_old = gamma_xh;
    xp_h = xp_h+theta*delx_vec;
    
    if flag == 1
        if e0h == e0
            AtgxAnx = A(:,gamma_xh)'*A(:,itheta)+b(:,gamma_xh)'*(epsilon*b(:,itheta));
            % btbgx = b(:,gamma_xh)'*b(:,gamma_xh);
            AtAgx_new = AtAgx + ((epsilon-e0)*b(:,gamma_xh)')*b(:,gamma_xh);
            iAtAgx_bt = iAtAgx*b(:,gamma_xh)';
            % iAtAgx_new = iAtAgx - (epsilon-e0)*(iAtAgx_bt)*(inv(1+(epsilon-e0)*b(:,gamma_xh)*iAtAgx_bt)*iAtAgx_bt');
            iAtAgx_new = iAtAgx - (((epsilon-e0)/(1+(epsilon-e0)*(b(:,gamma_xh)*iAtAgx_bt)))*iAtAgx_bt)*iAtAgx_bt';
            AtAgx_mod = [AtAgx_new AtgxAnx; AtgxAnx' A(:,itheta)'*A(:,itheta)+epsilon*b(:,itheta)'*b(:,itheta)];

            iAtAgx = update_inverse(AtAgx_mod, iAtAgx_new,1);
            AtAgx = AtAgx_mod;
            
            AtAgx_mod = [];
            AtAgx_new = [];
            iAtAgx_new = [];
            
            gamma_xh = [gamma_xh; itheta];
            xp_h(itheta) = 0;
        else
            AtgxAnx = A(:,gamma_xh)'*A(:,itheta)+b(:,gamma_xh)'*b(:,itheta);
            % btbgx = b(:,gamma_xh)'*b(:,gamma_xh);
            iAtAgx_bt = iAtAgx*b(:,gamma_xh)';
            AtAgx_mod = [AtAgx AtgxAnx; AtgxAnx' A(:,itheta)'*A(:,itheta)+b(:,itheta)'*b(:,itheta)];

            iAtAgx = update_inverse(AtAgx_mod, iAtAgx,1);
            AtAgx = AtAgx_mod;
            
            AtAgx_mod = [];
            gamma_xh = [gamma_xh; itheta];
            xp_h(itheta) = 0;
        end
    else
        outx_index = find(gamma_xh==itheta);
        gamma_xh(outx_index) = gamma_xh(end);
        gamma_xh(end) = itheta;
        gamma_xh = gamma_xh(1:end-1);

        if e0h == e0
            % btbgx = b(:,gamma_xh_old)'*b(:,gamma_xh_old);
            AtAgx_new = AtAgx + ((epsilon-e0)*b(:,gamma_xh_old)')*b(:,gamma_xh_old);
            iAtAgx_bt = iAtAgx*b(:,gamma_xh_old)';
            % iAtAgx_new = iAtAgx - (epsilon-e0)*(iAtAgx_bt)*(inv(1+(epsilon-e0)*b(:,gamma_xh_old)*iAtAgx_bt)*iAtAgx_bt');
            iAtAgx_new = iAtAgx - (((epsilon-e0)/(1+(epsilon-e0)*(b(:,gamma_xh_old)*iAtAgx_bt)))*iAtAgx_bt)*iAtAgx_bt';

            rowi = outx_index; % ith row of A is swapped with last row (out_x)
            colj = outx_index; % jth column of A is swapped with last column (out_x)
            AtAgx_ij = AtAgx_new;
            temp_row = AtAgx_ij(rowi,:);
            AtAgx_ij(rowi,:) = AtAgx_ij(end,:);
            AtAgx_ij(end,:) = temp_row;
            temp_col = AtAgx_ij(:,colj);
            AtAgx_ij(:,colj) = AtAgx_ij(:,end);
            AtAgx_ij(:,end) = temp_col;
            iAtAgx_ij = iAtAgx_new;
            temp_row = iAtAgx_ij(colj,:);
            iAtAgx_ij(colj,:) = iAtAgx_ij(end,:);
            iAtAgx_ij(end,:) = temp_row;
            temp_col = iAtAgx_ij(:,rowi);
            iAtAgx_ij(:,rowi) = iAtAgx_ij(:,end);
            iAtAgx_ij(:,end) = temp_col;

            AtAgx = AtAgx_ij(1:end-1,1:end-1);
            iAtAgx = update_inverse(AtAgx_ij, iAtAgx_ij,2);
            
            AtAgx_mod = [];
            AtAgx_new = [];
            iAtAgx_new = [];
            AtAgx_ij = [];
            iAtAgx_ij = [];
        else

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
            
            AtAgx_ij = [];
            iAtAgx_ij = [];
        end

        xp_h(itheta) = 0;
    end

    e0 = epsilon;
    pk_old = pk+theta*dk;
    pk_old([gamma_xh; itheta]) =  sign(pk_old([gamma_xh; itheta]))*tau;
end


