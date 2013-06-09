% DynamicX_update_function.m
%
% Solves
% min_x  ||x||_1 s.t. \|A'*(Ax-(1-\epsilon)y-(\epsilon)yt)\|_\infty <= \tau
%
% using homotopy update method
%
% [xp_h, lambda_h, gamma_xh, gamma_lh, iter, th] = DynamicX_DS_function(A, Q_glgx, R_glgx, y, yt, xp, lame, gamma_x, gamma_lambda, pk, ak, tau, maxiter)
%
% A: mxn matrix
% [Q_glgx, R_glgx] = qr factors for A_glgx
% y: old set of measurements
% yt: new measurements
% xp: old primal solution
% lame: old dual solution
% gamma_x: old primal support
% gamma_lambda: old dual support
% pk: old primal constraints
% ak: old dual constraints
% tau: threshold parameter. 

% DS update with dynamic change in x
% Author: Muhammad Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
% Created: February 2009
% Modified; June 2009
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

function [xp_h, lambda_h, gamma_xh, gamma_lh, iter, th] = DynamicX_DS_function(A, Q_glgx, R_glgx, y, yt, xp, lame, gamma_x, gamma_lambda, pk, ak, tau, maxiter)

t0 = cputime;
N = length(xp);
e0 = 0;
gamma_xh = gamma_x;
gamma_lh = gamma_lambda;
xp_h = xp;
lambda_h = lame;
idelta = [];
itheta = [];
epsilon = 0;
pk_old = pk;
pk_old(gamma_lh) = sign(pk_old(gamma_lh))*tau;
ak_old = ak;
ak_old(gamma_xh) = sign(ak_old(gamma_xh));
out_lambda = [];
out_x = [];
in_x = [];
in_lambda = [];

iter = 0;

d = A'*(y-yt);

% figure(10); clf; plot(A'*(A*xp-y)); hold on
% figure(11); clf; plot(A'*A*lame); hold on

while iter < maxiter
    iter = iter +1 ;
    in_lambda = [];
    out_lambda = [];
    in_x = [];
    out_x = [];
    
    %%%%%%%%%%%%%%%%%%%%%
    %%% Primal update %%%
    %%%%%%%%%%%%%%%%%%%%%

    %delx = -inv(A(:,gamma_lh)'*A(:,gamma_xh))*d(gamma_lh);
    delx = -R_glgx\(Q_glgx'*d(gamma_lh));
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;

    pk = pk_old;
    %dk = A'*(A*delx_vec+y-yt);
    dk = A'*(A(:,gamma_xh)*delx)+d;
    %dk = AtAf(delx_vec)+d;

    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW. 
    pk_temp = pk;
    gammaL_temp = find(abs(abs(pk)-tau)<min(epsilon,1e-12));
    pk_temp(gammaL_temp) = sign(pk(gammaL_temp))*tau;
    
    temp_gamma = zeros(N,1);
    temp_gamma(gamma_lh) = gamma_lh;
    gamma_lc = find([1:N]' ~= temp_gamma);
    % gamma_lc = setdiff([1:N],[gamma_lh]);

    b_constr1 = (tau-pk_temp(gamma_lc))./dk(gamma_lc);
    b_constr2 = (tau+pk_temp(gamma_lc))./-dk(gamma_lc);
    b_constr3 = (-xp_h(gamma_xh)./delx_vec(gamma_xh));
    itheta_1 = find(b_constr1>2*eps);
    itheta_2 = find(b_constr2>2*eps);
    itheta_3 = find(b_constr3>2*eps);
    theta1 = min(b_constr1(itheta_1));
    theta2 = min(b_constr2(itheta_2));
    theta3 = min(b_constr3(itheta_3));
    if theta1>theta2
        theta = theta2;
        itheta = gamma_lc(find(b_constr2==theta2));
        flag_theta = 1;
    else
        theta = theta1;
        itheta = gamma_lc(find(b_constr1==theta1));
        flag_theta = 1;
    end
    if theta3 < theta
        theta = theta3;
        itheta = gamma_xh(find(b_constr3==theta3));
        flag_theta = 0;
    end
    % [theta itheta flag_theta]
    epsilon = (theta+e0);
    if epsilon < 0
        epsilon= inf;
    end
    gamma_lh_old = gamma_lh;
    gamma_xh_old = gamma_xh;
    xp_old = xp_h;
    xp_h = xp_h+theta*delx_vec;
    if flag_theta == 1
        % a new element is added to the dual support
        new_lambda = itheta;
        gamma_lh_app = [gamma_lh; itheta];
        in_lambda = itheta;
    else
        % an element is removed from the primal support
        outx_index = find(gamma_xh==itheta);
        gamma_xh = gamma_xh([1:outx_index-1 outx_index+1:length(gamma_xh)]);

        outl_index = length(gamma_lh);
        new_lambda = gamma_lh(outl_index);
        gamma_lh = gamma_lh([1:outl_index-1 outl_index+1:length(gamma_lh)]);
        gamma_lh_app = [gamma_lh; new_lambda];

        out_x = itheta;
        xp_h(itheta) = 0;

        %%% UPDATE THE MATRIX AND ITS INVERSE HERE
        % out_x and new_lambda removed
        Q0 = Q_glgx;
        R0 = R_glgx;
        loc_vector = zeros(length(gamma_lh_old),1);
        loc_vector(outl_index) = 1;
        rep_row = (A(:,new_lambda)'*A(:,gamma_xh_old))';
        [Q1 R1] = qrupdate(Q0, R0, loc_vector,-rep_row);
        loc_vector = zeros(length(gamma_xh_old),1);
        loc_vector(outx_index) = 1;
        rep_vec = A(:,gamma_lh_old)'*A(:,out_x);
        rep_vec(outl_index) = 1;
        [Q2 R2] = qrupdate(Q1, R1, -rep_vec, loc_vector);
        Q_glgx = Q2([1:outl_index-1 outl_index+1:length(gamma_lh_old)],[1:outx_index-1 outx_index+1:length(gamma_xh_old)]);
        R_glgx = R2([1:outx_index-1 outx_index+1:length(gamma_xh_old)],[1:outx_index-1 outx_index+1:length(gamma_xh_old)]);
    end

    if epsilon >=1
        theta_end = (1-e0);
        xp_h = xp_old + theta_end*delx_vec;
        gamma_xh = gamma_xh_old;
        th = cputime-t0;
        break;
    end

    %%%%%%%%%%%%%%%%%%%
    %%% Dual update %%%
    %%%%%%%%%%%%%%%%%%%
    z_l = sign(pk(new_lambda)+theta*dk(new_lambda));
    dl = A(:,gamma_xh)'*A(:,new_lambda)*z_l;

    %del_lambda = -inv(A(:,gamma_xh)'*A(:,gamma_lh))*dl;
    del_lambda =  -Q_glgx*(R_glgx'\dl);
    del_lam_vec = zeros(N,1);
    del_lam_vec(gamma_lh) = del_lambda;
    del_lam_vec(new_lambda) = z_l;

    ak = ak_old;
    bk = A'*(A*del_lam_vec);
    if flag_theta == 0
        if sign(ak(out_x))==sign(bk(out_x))
            z_l = -z_l;
            del_lam_vec = -del_lam_vec;
            bk = -bk;
        end
    end
    ak(out_x) = sign(ak(out_x));
    
    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW.
    ak_temp = ak;
    gammaX_temp = find(abs(abs(ak)-1)<1e-12);
    ak_temp(gammaX_temp) = sign(ak(gammaX_temp));
    
    % gamma_xc = setdiff([1:N],[gamma_xh]);
    temp_gamma = zeros(N,1);
    temp_gamma(gamma_xh) = gamma_xh;
    gamma_xc = find([1:N]' ~= temp_gamma);

    b_constr1 = (1-ak_temp(gamma_xc))./bk(gamma_xc);
    b_constr2 = (1+ak_temp(gamma_xc))./-bk(gamma_xc);
    b_constr3 = (-lambda_h(gamma_lh_app)./del_lam_vec(gamma_lh_app));
    idelta_1 = find(b_constr1>2*eps);
    idelta_2 = find(b_constr2>2*eps);
    idelta_3 = find(b_constr3>2*eps);
	delta1 = min(b_constr1(idelta_1));
    delta2 = min(b_constr2(idelta_2));
    delta3 = min(b_constr3(idelta_3));

    if delta1>delta2
        delta = delta2;
        idelta = gamma_xc(find(b_constr2==delta2));
        flag_delta = 1;
    else
        delta = delta1;
        idelta = gamma_xc(find(b_constr1==delta1));
        flag_delta = 1;
    end
    if delta3 < delta
        delta = delta3;
        idelta = gamma_lh_app(find(b_constr3==delta3));
        flag_delta = 0;
    end

    lambda_old = lambda_h;
    lambda_h = lambda_h + delta*del_lam_vec;
    if flag_delta == 1
        gamma_lh_old = gamma_lh;
        gamma_xh_old = gamma_xh;
        gamma_xh = [gamma_xh; idelta];
        gamma_lh = gamma_lh_app;
        in_x = idelta;

        %%% UPDATE THE MATRIX AND ITS INVERSE HERE
        % new_lambda and in_x added
        Q0 = Q_glgx;
        R0 = R_glgx;
        Q0t = [Q0 zeros(length(gamma_xh_old),1); zeros(1,length(gamma_lh_old)) 1];
        R0t = [R0 zeros(length(gamma_xh_old),1); zeros(1,length(gamma_lh_old)) 1];
        loc_vector = zeros(length(gamma_lh),1);
        loc_vector(end) = 1;
        rep_row = (A(:,new_lambda)'*A(:,gamma_xh))';
        [Q2t R2t] = qrupdate(Q0t, R0t, loc_vector, rep_row);
        rep_vec = A(:,gamma_lh)'*A(:,in_x);
        rep_vec(end) = -1;
        [Q_glgx R_glgx] = qrupdate(Q2t, R2t, rep_vec, loc_vector);
    else
        % gamma_lh = setdiff(gamma_lh_app,idelta);
        outl_index = find(gamma_lh==idelta);
        out_lambda = idelta;

        out_lambda = idelta;
        lambda_h(idelta) = 0;
        if ~isempty(outl_index) % Because otherwise new_lambda is removed
            gamma_lh_old = gamma_lh;
            gamma_xh_old = gamma_xh;
            gamma_lh(outl_index) = new_lambda;

            %%% UPDATE THE MATRIX AND ITS INVERSE HERE
            % new_lambda and out_lambda swapped
            Q0 = Q_glgx;
            R0 = R_glgx;
            loc_vector = zeros(length(gamma_xh),1);
            loc_vector(outl_index) = 1;
            row_replaced =  A(:,new_lambda)'*A(:,gamma_xh)-A(:,out_lambda)'*A(:,gamma_xh);
            [Q_glgx, R_glgx] = qrupdate(Q0,R0, loc_vector,row_replaced');
        end

    end

    e0 = epsilon;
    pk_old = pk+theta*dk;
    pk_old([gamma_lh; out_lambda]) = sign(pk_old([gamma_lh; out_lambda]))*tau;
    ak_old = ak+delta*bk;
    ak_old(gamma_xh) = sign(ak_old(gamma_xh));
    
    [itheta theta flag_theta idelta delta flag_delta iter epsilon]

%     figure(10);
%     plot(pk_old,'r.'); shg
%     figure(11);
%     plot(ak_old,'r.'); shg; drawnow;

end

th = cputime-t0;
