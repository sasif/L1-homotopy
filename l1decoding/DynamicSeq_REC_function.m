% DynamicSeq_REC_function.m
% Robust ell_1 decoding with homotopy update scheme for new measurements.
%
% [cp_h gamma_h iter th] = DynamicSeq_REC_function(P, s, iPgP, cp_old, gamma_h, gamma_n, pk, tau, M, m_u, maxiter);
%
% Outputs:
% cp_h: sparse error estimate
% gamma_h: support of sparse errors
% iter: number of homotopy steps taken
%
% Inputs
% P: I-F*inv(F'F)*F'
% s: [y; w]
% iPgP: inv(PgP) where PgP = P(gamma_xh,gamma_xh)
% cp_old: Old estimate of error (includes the estimate for new measurements)
% gamma_h: support of cp_old
% gamma_n: support of estimated errors in new measurements
% pk: P*(cp_old-s)
% tau: regularization parameter
% M: Number of old measurements
% m_u: Number of new measurements
% maxiter: Maximum allowed iterations
%
% Created by Salman Asif @ Georgia Tech
% February 2009
%
%-------------------------------------------+
% Copyright (c) 2009.  Muhammad Salman Asif 
%-------------------------------------------+


function [cp_h gamma_h iter th] = REC_homotopy_function(P, s, iPgP, cp_old, gamma_h, gamma_n, pk, tau, M, m_u, maxiter);

t0 = cputime;
cp_h = cp_old;
epsilon = 0;
e0 = 0;
pk_old = pk; % last indices will be zero at the start.

pk_old([gamma_h]) = sign(pk_old([gamma_h]))*tau;
pk_old(gamma_n) = pk_old(gamma_n)*e0*tau;

done = 0;
iter = 0;
tu0 = cputime;
outc = [];

while iter < maxiter
    iter = iter +1 ;

    PgP = P(gamma_h,gamma_h);

    c_k = cp_h;
    len_e = length(gamma_h)-length(gamma_n);
    z_d = sign(c_k(gamma_n));
    igamma_n_gamma_h = zeros(length(gamma_n),1);
    for ii = 1:length(gamma_n)
        igamma_n_gamma_h(ii) = find(gamma_n(ii) == gamma_h);
    end
    temp_vec = zeros(length(gamma_h),1);
    temp_vec(igamma_n_gamma_h) = z_d;
    delc = -iPgP*temp_vec;
    % figure(102); imagesc(iPgP-inv(P(gamma_h,gamma_h))); colorbar; shg
    % delc = -inv(P(gamma_h,gamma_h))*temp_vec;
    
    delc_vec = zeros(M+m_u,1);
    delc_vec(gamma_h) = delc;
    pk = pk_old;
    dk = P(:,gamma_h)*delc;
    outc = [];
    temp_gamma = zeros(M+m_u,1);
    temp_gamma(gamma_h) = gamma_h;
    gamma_hc = find([1:M+m_u]' ~= temp_gamma);

    b_constr1 = (tau-pk(gamma_hc))./dk(gamma_hc);
    b_constr2 = (tau+pk(gamma_hc))./-dk(gamma_hc);
    b_constr3 = (-c_k(gamma_h)./delc_vec(gamma_h));
    itheta_1 = find(b_constr1>2*eps);
    itheta_2 = find(b_constr2>2*eps);
    itheta_3 = find(b_constr3>2*eps);
    theta1 = min(b_constr1(itheta_1));
    theta2 = min(b_constr2(itheta_2));
    theta3 = min(b_constr3(itheta_3));
    if theta1>theta2
        theta = theta2;
        itheta = gamma_hc(find(b_constr2==theta2));
        flag = 1;
    else
        theta = theta1;
        itheta = gamma_hc(find(b_constr1==theta1));
        flag = 1;
    end

    if theta3 < theta
        theta = theta3;
        itheta = gamma_h(find(b_constr3==theta3));
        outc = itheta;
        flag = 0;
    end
    
    if length(gamma_h)>=M+m_u-128 & flag ==1
        stp = 1;
    end

    epsilon = theta/tau + e0;

    if epsilon < 0
        epsilon= inf;
    end
    gamma_h_old = gamma_h;
    gamma_n_old = gamma_n;
    cp_old = c_k;
    cp_h = c_k+theta*delc_vec;

    if epsilon >=1
        theta_end = tau*(1-e0);
        cp_h = cp_old + theta_end*delc_vec;
        gamma_h = gamma_h_old;
        pk_old = pk+theta_end*dk;
        break;
    end

    if flag == 1
        gamma_h = [gamma_h; itheta];
        cp_h(itheta) = 0;
        PgP_mod = P(gamma_h,gamma_h);

        iPgP = update_inverse(PgP_mod, iPgP,1);
        PgP = PgP_mod;

    else
        outc_index = find(gamma_h==itheta);
        gamma_h(outc_index) = gamma_h(end);
        gamma_h(end) = itheta;
        gamma_h = gamma_h(1:end-1);

        rowi = outc_index; % ith row of P is swapped with last row (out_c)
        colj = outc_index; % jth column of P is swapped with last column (out_c)
        PgP_ij = PgP;
        temp_row = PgP_ij(rowi,:);
        PgP_ij(rowi,:) = PgP_ij(end,:);
        PgP_ij(end,:) = temp_row;
        temp_col = PgP_ij(:,colj);
        PgP_ij(:,colj) = PgP_ij(:,end);
        PgP_ij(:,end) = temp_col;
        iPgP_ij = iPgP;
        temp_row = iPgP_ij(colj,:);
        iPgP_ij(colj,:) = iPgP_ij(end,:);
        iPgP_ij(end,:) = temp_row;
        temp_col = iPgP_ij(:,rowi);
        iPgP_ij(:,rowi) = iPgP_ij(:,end);
        iPgP_ij(:,end) = temp_col;

        PgP = PgP_ij(1:end-1,1:end-1);
        iPgP = update_inverse(PgP_ij, iPgP_ij,2);

        cp_h(itheta) = 0;
        if length(find(itheta == gamma_n))
            outn_index = find(gamma_n==itheta);
            gamma_n(outn_index) = gamma_n(end);
            gamma_n(end) = itheta;
            gamma_n = gamma_n(1:end-1);
            % gamma_n = setdiff(gamma_n,itheta);
            if isempty(gamma_n)
                e0_old = epsilon;
                epsilon = 1;
                pk_old = pk+theta*dk;
                break;
            end
        end
    end

    e0 = (theta/tau+e0);

    pk_old = pk+theta*dk;
    pk_old([gamma_h; itheta]) = sign(pk_old([gamma_h; itheta]))*tau;
    pk_old(gamma_n_old) = pk_old(gamma_n_old)*e0;

end
th= cputime-t0;