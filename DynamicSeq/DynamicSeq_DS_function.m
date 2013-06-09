% Dantzig Selector with sequential measurements using homotopy.
% 
% Solves  
% min_x  ||x||_1 s.t. ||A'(Ax-y)+b'(bx-w)||_\infty <= tau
%
% using homotopy update method
% [xp_h, lambda_h, gamma_xh, gamma_lh, iter, th] =
% DynamicSeq_DS_function(A, b, Q_glgx, R_glgx, y, w, xp, lame, gamma_x, gamma_lambda, pk, ak, tau, maxiter)
%
% A: mxn matrix
% b: new rows in measurement matrix
% [Q_glgx, R_glgx] = qr factors for A_glgx
% y: old set of measurements
% w: new measurements
% xp: old primal solution
% lame: old dual solution
% gamma_x: old primal support
% gamma_lambda: old dual support
% pk: old primal constraints
% ak: old dual constraints
% tau: threshold parameter. 

% Dantzig selector solution dynamica update with sequential measurements
% Rank one update is based on qr factorization.

% Author: Salman Asif, Georgia Tech.
% Email: sasif@ece.gatech.edu
% Created: June 2009
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+

function [xp_h, lambda_h, gamma_xh, gamma_lh, iter, th] = DynamicSeq_DS_function(A, b, Q_glgx, R_glgx, y, w, xp, lame, gamma_x, gamma_lambda, pk, ak, tau, maxiter)

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

while iter < maxiter
    iter = iter+1;
   
    e0h = e0; % making e0h==1 workes with BPDN but NOT with DS. 
    % The reason for this is because the value of step size is not
    % guaranteed to stay positive. 
    
    dx = b(:,gamma_lh)'*(b(:,gamma_xh)*xp_h(gamma_xh)-w);
    
    %%% NEED FAST UPDATE
    % iAglgx = inv(A(:,gamma_lh)'*A(:,gamma_xh)+e0h*b(:,gamma_lh)'*b(:,gamma_xh));
    % ux = (b(:,gamma_xh)*iAglgx*b(:,gamma_lh)');
    % delx = -iAglgx*dx;
    
    iAb_vec = R_glgx\(Q_glgx'*b(:,gamma_lh)');
    ux = (b(:,gamma_xh)*iAb_vec);
    delx = -R_glgx\(Q_glgx'*dx);

    pk = pk_old; % A'*(A*xp_h-y)+e0*b'*(b*xp_h-w);
    dk = b'*(b*xp_h-w)+(A'*(A(:,gamma_xh)*delx)+e0h*b'*(b(:,gamma_xh)*delx));
    delx_vec = zeros(N,1);
    delx_vec(gamma_xh) = delx;

    dl = b(:,gamma_xh)'*(b(:,gamma_lh)*lambda_h(gamma_lh));
    ul = ux;
   
    %%% NEED FAST UPDATE
    % iAgxgl = iAglgx';
    % del_lam = -iAgxgl*dl;
     
    del_lam = -Q_glgx*(R_glgx'\dl);
        
    ak = ak_old; %A'*A*lambda_h+e0*b'*b*lambda_h;
    bk = b'*(b*lambda_h)+(A'*(A(:,gamma_lh)*del_lam)+e0h*b'*(b(:,gamma_lh)*del_lam));
    del_lam_vec = zeros(N,1);
    del_lam_vec(gamma_lh) = del_lam;

    % gamma_lc = setdiff([1:N],[gamma_lh]);
    temp_gamma = zeros(N,1);
    temp_gamma(gamma_lh) = gamma_lh;
    gamma_lc = find([1:N]' ~= temp_gamma);

    d_constr1 = (tau-pk(gamma_lc))./dk(gamma_lc);
    d_constr2 = (tau+pk(gamma_lc))./-dk(gamma_lc);
    d_constr3 = (-xp_h(gamma_xh)./delx_vec(gamma_xh));
    idelta_1 = find(d_constr1>2*eps);
    idelta_2 = find(d_constr2>2*eps);
    idelta_3 = find(d_constr3>2*eps);
    delta1 = min(d_constr1(idelta_1));
    delta2 = min(d_constr2(idelta_2));
    delta3 = min(d_constr3(idelta_3));
    if isempty(delta1)
        delta1 = inf;
    end
    if isempty(delta2)
        delta2 = inf;
    end
    if isempty(delta3)
        delta3 = inf;
    end

    if delta1>delta2
        delta = delta2;
        idelta = gamma_lc(find(d_constr2==delta2));
        flag_delta = 1;
    else
        delta = delta1;
        idelta = gamma_lc(find(d_constr1==delta1));
        flag_delta = 1;
    end
    if delta3 < delta
        delta = delta3;
        idelta = gamma_xh(find(d_constr3==delta3));
        flag_delta = 0;
    end
    
    % gamma_xc = setdiff([1:N],[gamma_xh]);
    temp_gamma = zeros(N,1);
    temp_gamma(gamma_xh) = gamma_xh;
    gamma_xc = find([1:N]' ~= temp_gamma);

    b_constr1 = (1-ak(gamma_xc))./bk(gamma_xc);
    b_constr2 = (1+ak(gamma_xc))./-bk(gamma_xc);
    b_constr3 = (-lambda_h(gamma_lh)./del_lam_vec(gamma_lh));
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
        flag_theta = 1;
    else
        theta = theta1;
        itheta = gamma_xc(find(b_constr1==theta1));
        flag_theta = 1;
    end
    if theta3 < theta
        theta = theta3;
        itheta = gamma_lh(find(b_constr3==theta3));
        flag_theta = 0;
    end

    if 1-theta*ul <=0 & ul >0
        stp = 1;
    end
    if 1-delta*ux <=0 & ux >0
        stp = 1;
    end

    epsilon_x = (delta+e0-e0h*delta*ux)/(1-delta*ux);
    epsilon_l = (theta+e0-e0h*theta*ul)/(1-theta*ul);

    dteta = min(delta, theta);

    if epsilon_l <0
        epsilon_l = inf;
    end
    if epsilon_x < 0
        epsilon_x = inf;
    end

    epsilon = min(epsilon_x,epsilon_l);

    if epsilon >=1
        dteta_end = ((1-e0)/(1+ux*(1-e0h)));
        xp_h = xp_h + dteta_end*delx_vec;
        lambda_h = lambda_h + dteta_end*del_lam_vec;
        break;
    end

    gamma_xh_old = gamma_xh;
    gamma_lh_old = gamma_lh;
    xp_h_old = xp_h;
    xp_h = xp_h+dteta*delx_vec;
    lambda_h_old = lambda_h;
    lambda_h = lambda_h+dteta*del_lam_vec;
    
    pk_t = pk+dteta*dk;
    ak_t = ak+dteta*bk;

    in_lambda = [];
    out_lambda = [];
    in_x = [];
    out_x = [];
    if epsilon_x > epsilon_l
        %%%%%%%%%%%%%%%%%%%%%
        %%% x post-update %%%
        %%%%%%%%%%%%%%%%%%%%%
        if flag_theta == 1
            % when an element is added to Gamma_x by selected
            % epsilon_dual
            new_x = itheta;
            gamma_xh_app = [gamma_xh; new_x];
            in_x = itheta;
            
            %%% UPDATE THE MATRIX AND ITS INVERSE HERE
            % just update with new epsilon
            Q0 = Q_glgx;
            R0 = R_glgx;
            [Q_glgx R_glgx] = qrupdate(Q0, R0, (epsilon-e0)*b(:,gamma_lh)',b(:,gamma_xh)');
        else
            % when an element is removed from Gamma_lambda by selected
            % epsilon_dual
            outl_index = find(gamma_lh == itheta);
            % gamma_lh = setdiff(gamma_lh,itheta);
            gamma_lh = gamma_lh([1:outl_index-1 outl_index+1:length(gamma_lh)]);
            
            % [max_val outx_index] = max(abs(iAglgx(:,outl_index)));
            outx_index = length(gamma_xh);
            new_x = gamma_xh(outx_index);
            % gamma_xh = setdiff(gamma_xh,new_x);
            gamma_xh = gamma_xh([1:outx_index-1 outx_index+1:length(gamma_xh)]);
            gamma_xh_app = [gamma_xh; new_x];
            
            lambda_h(itheta) = 0;
            out_lambda = itheta;
            
            %%% UPDATE THE MATRIX AND ITS INVERSE HERE
            % out_lambda and new_x removed, and new epsilon
            Q0 = Q_glgx;
            R0 = R_glgx;
            [Q1 R1] = qrupdate(Q0, R0, (epsilon-e0)*b(:,gamma_lh_old)',b(:,gamma_xh_old)');
            loc_vector = zeros(length(gamma_lh_old),1);
            loc_vector(outl_index) = 1;
            rep_row = (A(:,out_lambda)'*A(:,gamma_xh_old)+epsilon*b(:,out_lambda)'*b(:,gamma_xh_old))';
            [Q2t R2t] = qrupdate(Q1, R1, loc_vector,-rep_row);
            loc_vector = zeros(length(gamma_xh_old),1);
            loc_vector(outx_index) = 1;
            rep_vec = (A(:,gamma_lh_old)'*A(:,new_x)+epsilon*b(:,gamma_lh_old)'*b(:,new_x));
            rep_vec(outl_index) = 1;
            [Q3t R3t] = qrupdate(Q2t, R2t, -rep_vec, loc_vector);
            Q_glgx = Q3t([1:outl_index-1 outl_index+1:length(gamma_lh_old)],[1:outx_index-1 outx_index+1:length(gamma_xh_old)]);
            R_glgx = R3t([1:outx_index-1 outx_index+1:length(gamma_xh_old)],[1:outx_index-1 outx_index+1:length(gamma_xh_old)]);
        
            % [Q2 R2] = qrdelete(Q1, R1, outl_index,'row');
            % [Q_glgx R_glgx] = qrdelete(Q2, R2, outx_index, 'col');
        end
        z_x = -sign(ak_t(new_x));
    
        dx = ((A(:,gamma_lh)'*A(:,new_x)+epsilon*b(:,gamma_lh)'*b(:,new_x))*z_x);
        
        %%% NEED FAST UPDATE
        % iAglgx_btb = inv(A(:,gamma_lh)'*A(:,gamma_xh)+epsilon*b(:,gamma_lh)'*b(:,gamma_xh));
        % del_x_hat = -iAglgx_btb*dx;

        del_x_hat = -R_glgx\(Q_glgx'*dx);

        del_xhat_vec = zeros(N,1);
        del_xhat_vec(gamma_xh) = del_x_hat;
        del_xhat_vec(new_x) = z_x;

        dk_t = A'*(A*del_xhat_vec)+epsilon*b'*(b*del_xhat_vec);
        if flag_theta == 0
            if sign(pk_t(itheta))== sign(dk_t(itheta))
                z_x = -z_x;
                del_xhat_vec = -del_xhat_vec;
                dk_t = -dk_t;
            end
        end
        pk_t(out_lambda) = sign(pk_t(out_lambda))*tau;
        % gamma_lc = setdiff([1:N],[gamma_lh]);
        temp_gamma = zeros(N,1);
        temp_gamma(gamma_lh) = gamma_lh;
        gamma_lc = find([1:N]' ~= temp_gamma);

        d_constr1 = (tau-pk_t(gamma_lc))./dk_t(gamma_lc);
        d_constr2 = (tau+pk_t(gamma_lc))./-dk_t(gamma_lc);
        d_constr3 = (-xp_h(gamma_xh_app)./del_xhat_vec(gamma_xh_app));
        idelta_1 = find(d_constr1>2*eps);
        idelta_2 = find(d_constr2>2*eps);
        idelta_3 = find(d_constr3>2*eps);
        delta1 = min(d_constr1(idelta_1));
        delta2 = min(d_constr2(idelta_2));
        delta3 = min(d_constr3(idelta_3));

        if delta1>delta2
            delta = delta2;
            idelta = gamma_lc(find(d_constr2==delta2));
            flag_delta = 1;
        else
            delta = delta1;
            idelta = gamma_lc(find(d_constr1==delta1));
            flag_delta = 1;
        end
        if delta3 < delta
            delta = delta3;
            idelta = gamma_xh_app(find(d_constr3==delta3));
            flag_delta = 0;
        end
        
        xp_h_old = xp_h;
        xp_h = xp_h + delta*del_xhat_vec;

        if flag_delta == 1
            gamma_xh_old = gamma_xh;
            gamma_lh_old = gamma_lh;
            gamma_lh = [gamma_lh; idelta];
            gamma_xh = gamma_xh_app;
            in_lambda = idelta;
            
            %%% UPDATE THE MATRIX AND ITS INVERSE HERE FOR
            % new_x and new_lambda added;
            Q0 = Q_glgx;
            R0 = R_glgx;
            Q0t = [Q0 zeros(length(gamma_xh_old),1); zeros(1,length(gamma_lh_old)) 1];
            R0t = [R0 zeros(length(gamma_xh_old),1); zeros(1,length(gamma_lh_old)) 1];
            loc_vector = zeros(length(gamma_lh),1);
            loc_vector(end) = 1;
            rep_row = (A(:,in_lambda)'*A(:,gamma_xh)+epsilon*b(:,in_lambda)'*b(:,gamma_xh))';
            [Q2t R2t] = qrupdate(Q0t, R0t, loc_vector, rep_row);
            rep_vec = A(:,gamma_lh)'*A(:,new_x)+epsilon*b(:,gamma_lh)'*b(:,new_x);
            rep_vec(end) = -1;
            [Q_glgx R_glgx] = qrupdate(Q2t, R2t, rep_vec, loc_vector);

            % [Q1, R1] = qrinsert(Q0, R0, length(gamma_lh), A(:,gamma_lh_old)'*A(:,new_x)+epsilon*b(:,gamma_lh_old)'*b(:,new_x),'col');
            % [Q_glgx, R_glgx] = qrinsert(Q1,R1, length(gamma_lh), A(:,in_lambda)'*A(:,gamma_xh)+epsilon*b(:,in_lambda)'*b(:,gamma_xh),'row');
        else
            gamma_xh_old = gamma_xh;
            gamma_lh_old = gamma_lh;
            
            % gamma_xh = setdiff(gamma_xh_app,idelta);
            outx_index = find(gamma_xh==idelta);
            out_x = idelta;
            xp_h(out_x) = 0;
            
            if ~isempty(outx_index) % Because otherwise new_x is removed
                gamma_xh(outx_index) = new_x;
    
                %%% UPDATE THE MATRIX AND ITS INVERSE HERE
                % new_x and out_x swapped
                Q0 = Q_glgx;
                R0 = R_glgx;
                loc_vector = zeros(length(gamma_xh),1);
                loc_vector(outx_index) = 1;
                [Q_glgx, R_glgx] = qrupdate(Q0,R0, A(:,gamma_lh)'*A(:,new_x)+epsilon*b(:,gamma_lh)'*b(:,new_x)-(A(:,gamma_lh)'*A(:,out_x)+epsilon*b(:,gamma_lh)'*b(:,out_x)),loc_vector);
            end
        end

        pk_old = pk_t+delta*dk_t;
        pk_old(gamma_lh) = sign(pk_old(gamma_lh))*tau;
        ak_old = ak_t;
        ak_old(gamma_xh_app) = sign(ak_old(gamma_xh_app));
    else 
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% lambda post-update %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if flag_delta == 1
            % when an element is added to Gamma_lambda by selected
            % epsilon_primal
            new_lambda = idelta;
            gamma_lh_app = [gamma_lh; new_lambda];
            in_lambda = idelta;
            %%% UPDATE THE MATRIX AND ITS INVERSE HERE
            % just update the epsilon
            Q0 = Q_glgx;
            R0 = R_glgx;
            [Q_glgx R_glgx] = qrupdate(Q0, R0, (epsilon-e0)*b(:,gamma_lh)',b(:,gamma_xh)');
        else
            gamma_lh_old = gamma_lh;
            gamma_xh_old = gamma_xh;
            
            % when an element is removed from Gamma_x by selected
            % epsilon_primal
            outx_index = find(gamma_xh==idelta);
            % gamma_xh = setdiff(gamma_xh,idelta);
            gamma_xh = gamma_xh([1:outx_index-1 outx_index+1:length(gamma_xh)]);
            
            % [max_val outl_index] = max(abs(iAgxgl(:,outx_index)));
            outl_index = length(gamma_lh);
            new_lambda = gamma_lh(outl_index);
            % gamma_lh = setdiff(gamma_lh,new_lambda);
            gamma_lh = gamma_lh([1:outl_index-1 outl_index+1:length(gamma_lh)]);
            gamma_lh_app = [gamma_lh; new_lambda];
            xp_h(idelta) = 0;
            out_x = idelta;
            
            %%% UPDATE THE MATRIX AND ITS INVERSE HERE
            % out_x and new_lambda removed and new epsilon
            Q0 = Q_glgx;
            R0 = R_glgx;
            [Q1 R1] = qrupdate(Q0, R0, (epsilon-e0)*b(:,gamma_lh_old)',b(:,gamma_xh_old)');
            loc_vector = zeros(length(gamma_lh_old),1);
            loc_vector(outl_index) = 1;
            rep_row = (A(:,new_lambda)'*A(:,gamma_xh_old)+epsilon*b(:,new_lambda)'*b(:,gamma_xh_old))';
            [Q2t R2t] = qrupdate(Q1, R1, loc_vector,-rep_row);
            loc_vector = zeros(length(gamma_xh_old),1);
            loc_vector(outx_index) = 1;
            rep_vec = (A(:,gamma_lh_old)'*A(:,out_x)+epsilon*b(:,gamma_lh_old)'*b(:,out_x));
            rep_vec(outl_index) = 1; 
            [Q3t R3t] = qrupdate(Q2t, R2t, -rep_vec, loc_vector);
            Q_glgx = Q3t([1:outl_index-1 outl_index+1:length(gamma_lh_old)],[1:outx_index-1 outx_index+1:length(gamma_xh_old)]);
            R_glgx = R3t([1:outx_index-1 outx_index+1:length(gamma_xh_old)],[1:outx_index-1 outx_index+1:length(gamma_xh_old)]);
            
            % [Q2 R2] = qrdelete(Q1, R1, outl_index,'row');
            % [Q_glgx R_glgx] = qrdelete(Q2, R2, outx_index, 'col'); 
        end
        z_l = sign(pk_t(new_lambda));
        
        dl = ((A(:,gamma_xh)'*A(:,new_lambda)+epsilon*b(:,gamma_xh)'*b(:,new_lambda))*z_l);
        
        %%% NEED FAST UPDATE
        % iAgxgl_btb = inv(A(:,gamma_xh)'*A(:,gamma_lh)+epsilon*b(:,gamma_xh)'*b(:,gamma_lh));
        % del_lam_hat = -iAgxgl_btb*dl;
         
        del_lam_hat = -Q_glgx*(R_glgx'\dl);
        
        del_lhat_vec = zeros(N,1);
        del_lhat_vec(gamma_lh) = del_lam_hat;
        del_lhat_vec(new_lambda) = z_l;

        bk_t = A'*(A*del_lhat_vec)+epsilon*b'*(b*del_lhat_vec);
        if flag_delta == 0
            if sign(ak_t(idelta))==sign(bk_t(idelta))
                z_l = -z_l;
                del_lhat_vec = -del_lhat_vec;
                bk_t = -bk_t;
            end
        end
        ak_t(out_x) = sign(ak_t(out_x));
        % gamma_xc = setdiff([1:N],[gamma_xh]);
        temp_gamma = zeros(N,1);
        temp_gamma(gamma_xh) = gamma_xh;
        gamma_xc = find([1:N]' ~= temp_gamma);

        b_constr1 = (1-ak_t(gamma_xc))./bk_t(gamma_xc);
        b_constr2 = (1+ak_t(gamma_xc))./-bk_t(gamma_xc);
        b_constr3 = (-lambda_h(gamma_lh_app)./del_lhat_vec(gamma_lh_app));
        itheta_1 = find(b_constr1>2*eps);
        itheta_2 = find(b_constr2>2*eps);
        itheta_3 = find(b_constr3>2*eps);
        theta1 = min(b_constr1(itheta_1));
        theta2 = min(b_constr2(itheta_2));
        theta3 = min(b_constr3(itheta_3));

        if theta1>theta2
            theta = theta2;
            itheta = gamma_xc(find(b_constr2==theta2));
            flag_theta = 1;
        else
            theta = theta1;
            itheta = gamma_xc(find(b_constr1==theta1));
            flag_theta = 1;
        end
        if theta3 < theta
            theta = theta3;
            itheta = gamma_lh_app(find(b_constr3==theta3));
            flag_theta = 0;
        end

        lambda_h_old = lambda_h;
        lambda_h = lambda_h + theta*del_lhat_vec;
        if flag_theta == 1
            gamma_lh_old = gamma_lh;
            gamma_xh_old = gamma_xh;
            gamma_xh = [gamma_xh; itheta];
            gamma_lh = gamma_lh_app;
            in_x = itheta;
            
            %%% UPDATE THE MATRIX AND ITS INVERSE HERE
            % new_lambda and in_x added
            Q0 = Q_glgx;
            R0 = R_glgx;
            Q0t = [Q0 zeros(length(gamma_xh_old),1); zeros(1,length(gamma_lh_old)) 1];
            R0t = [R0 zeros(length(gamma_xh_old),1); zeros(1,length(gamma_lh_old)) 1];
            loc_vector = zeros(length(gamma_lh),1);
            loc_vector(end) = 1;
            rep_row = (A(:,new_lambda)'*A(:,gamma_xh)+epsilon*b(:,new_lambda)'*b(:,gamma_xh))';
            [Q2t R2t] = qrupdate(Q0t, R0t, loc_vector, rep_row);
            rep_vec = A(:,gamma_lh)'*A(:,in_x)+epsilon*b(:,gamma_lh)'*b(:,in_x);
            rep_vec(end) = -1;
            [Q_glgx R_glgx] = qrupdate(Q2t, R2t, rep_vec, loc_vector);

            % [Q1, R1] = qrinsert(Q0, R0, length(gamma_xh), A(:,gamma_lh_old)'*A(:,in_x)+epsilon*b(:,gamma_lh_old)'*b(:,in_x),'col');
            % [Q_glgx, R_glgx] = qrinsert(Q1,R1, length(gamma_xh), A(:,new_lambda)'*A(:,gamma_xh)+epsilon*b(:,new_lambda)'*b(:,gamma_xh),'row');
        else
            % gamma_lh = setdiff(gamma_lh_app,itheta);
            outl_index = find(gamma_lh==itheta);
            out_lambda = itheta;
            
            out_lambda = itheta;
            lambda_h(itheta) = 0;
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
                row_replaced =  (A(:,new_lambda)'*A(:,gamma_xh)+epsilon*b(:,new_lambda)'*b(:,gamma_xh)-(A(:,out_lambda)'*A(:,gamma_xh)+epsilon*b(:,out_lambda)'*b(:,gamma_xh)));
                [Q_glgx, R_glgx] = qrupdate(Q0,R0, loc_vector,row_replaced');
            end

        end
        ak_old = ak_t+theta*bk_t;
        ak_old(gamma_xh) = sign(ak_old(gamma_xh));
        pk_old = pk_t;
        pk_old(gamma_lh_app) = sign(pk_old(gamma_lh_app))*tau;
    end
    e0 = epsilon;
end
th = cputime-t0;