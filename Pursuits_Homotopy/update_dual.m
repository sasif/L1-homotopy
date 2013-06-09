% update_dual.m
%
% This function computes the minimum step size in the dual update direction and
% finds change in the primal or dual support with that step.
% 
% Inputs: 
% gamma_x - current support of x
% gamma_lambda - current support of lambda
% z_x - sign sequence of x
% z_lambda - sign sequence of lambda
% del_lambda_p - dual update direction
% ak 
% bk
% new_lambda - element entered in the support of lambda during primal update
% out_x - element of x shrunk to zero during primal update phase.
%
% Outputs: 
% i_theta - index corresponding to newly active dual constraint (new_x)
% out_lambda - element in lambda shrunk to zero
% theta - dual step size 
% chk_lambda - 1  an element is removed from support of lambda
%              0  a new element enters the support of x
% 
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu

function [i_theta, out_lambda, theta, chk_lambda] = update_dual(gamma_x, gamma_lambda, z_lambda, lambda_k, del_lambda_p, ak, bk, new_lambda, out_x);

N = length(lambda_k);

% gamma_xc = setdiff([1:N]', [gamma_x; out_x]); WRONG
% check out_x as well, that is if outgoing x switches sign in just one step
temp_gamma = zeros(N,1);
temp_gamma(gamma_x) = gamma_x;
gamma_xc = find([1:N]' ~= temp_gamma);

theta1_constr = (1-ak(gamma_xc))./bk(gamma_xc);
theta1_pos_index = find(theta1_constr>0);
theta1_pos = theta1_constr(theta1_pos_index);
[theta1 i_theta1] = min(theta1_pos);
if isempty(theta1)
    theta1 = inf;
end
theta2_constr = -(1+ak(gamma_xc))./bk(gamma_xc);
theta2_pos_index = find(theta2_constr>0);
theta2_pos = theta2_constr(theta2_pos_index);
[theta2 i_theta2] = min(theta2_pos);
if isempty(theta2)
    theta2 = inf;
end

if theta1 > theta2
    theta = theta2;
    i_theta = gamma_xc(theta2_pos_index(i_theta2));
else
    theta = theta1;
    i_theta = gamma_xc(theta1_pos_index(i_theta1));
end

gamma_lambda_app = [gamma_lambda; new_lambda];
theta3_constr = (-lambda_k(gamma_lambda_app)./del_lambda_p(gamma_lambda_app));
theta3_pos_index = find(theta3_constr>0);
[theta3 i_theta3] = min(theta3_constr(theta3_pos_index));
out_lambda_index = gamma_lambda_app(theta3_pos_index(i_theta3));

chk_lambda = 0;
out_lambda = [];
if theta3 > 0 & theta3<theta
    chk_lambda = 1;
    theta = theta3;
    out_lambda = out_lambda_index;
end


%%% THESE ARE PROBABLY UNNECESSARY 
%%% NEED TO REMOVE THEM. 

% This one is ONLY for those indices which are zero. And we don't know where
% will its dlambda point in next steps, so after we calculate dlambda and its in opposite
% direction to z_lambda, we will have to remove that index from the support.

lambdak_1 = lambda_k+theta*del_lambda_p;
lambdak_1(out_lambda) = 0;
wrong_sign = find(sign(lambdak_1(gamma_lambda)).*z_lambda(gamma_lambda)==-1);
if ~isempty(gamma_lambda(wrong_sign))
    chk_lambda = 1;
    theta = 0;
    out_lambda = gamma_lambda(wrong_sign(1));
end

% The following checks are just to deal with degenerate cases when more
% than one elements want to enter or leave the support at any step 
% (e.g., Bernoulli matrix with small number of measurements)

% This happens if more than one dual constraints become active in one step
% so some of the new elements in support of x got missed, here we check if
% they are still active.
i_theta_temp = gamma_xc(find(abs(ak(gamma_xc)+theta*bk(gamma_xc))-1 >= 10*eps));
if ~isempty(out_x)
    i_theta_more = i_theta_temp(find(i_theta_temp~=out_x));
else
    i_theta_more = i_theta_temp;
end
if ~isempty(i_theta_more)
    theta = 0;
    i_theta = i_theta_more(1);
    out_lambda=[];
    chk_lambda=0;
end
