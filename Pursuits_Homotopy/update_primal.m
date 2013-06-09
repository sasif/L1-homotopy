% update_primal.m
%
% This function computes the minimum step size in the primal update direction and
% finds change in the primal or dual support with that step.
%
% Inputs:
% gamma_x - current support of x
% gamma_lambda - current support of lambda
% z_x - sign sequence of x
% z_lambda - sign sequence of lambda
% del_x_vec - primal update direction
% pk
% dk
% epsilon - current value of epsilon
% out_lambda - element removed from support of lambda in previous step (if any)
%
% Outputs:
% i_delta - index corresponding to newly active primal constraint (new_lambda)
% out_x - element in x shrunk to zero
% delta - primal step size
% chk_x - 1  an element is removed from support of x
%         0  a new element enters the support of lambda
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu

function [i_delta, out_x, delta, chk_x] = update_primal(gamma_x, gamma_lambda, z_x, x_k, del_x_vec, pk, dk, epsilon, out_lambda);

N = length(x_k);

% gamma_lc = setdiff([1:N]', [gamma_lambda; out_lambda]); WRONG
% check out_lambda as well, that is if outgoing lambda switches sign in just one step
temp_gamma = zeros(N,1);
temp_gamma(gamma_lambda) = gamma_lambda;
gamma_lc = find([1:N]' ~= temp_gamma);

delta1_constr = (epsilon-pk(gamma_lc))./(1+dk(gamma_lc));
delta1_pos_ind = find(delta1_constr>0);
delta1_pos = delta1_constr(delta1_pos_ind);
[delta1 i_delta1] = min(delta1_pos);
if isempty(delta1)
    delta1 = inf;
end
delta2_constr = (epsilon+pk(gamma_lc))./(1-dk(gamma_lc));
delta2_pos_ind = find(delta2_constr>0);
delta2_pos = delta2_constr(delta2_pos_ind);
[delta2 i_delta2] = min(delta2_pos);
if isempty(delta2)
    delta2 = inf;
end

if delta1>delta2
    delta = delta2;
    i_delta = gamma_lc(delta2_pos_ind(i_delta2));
else
    delta = delta1;
    i_delta = gamma_lc(delta1_pos_ind(i_delta1));
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


% If program goes in any of the following loops, there is some problem!!!
%
%%% THESE ARE PROBABLY UNNECESSARY 
%%% NEED TO REMOVE THEM. 
%
% The following checks are just to deal with degenerate cases when more
% than one elements want to enter or leave the support at any step
% (e.g., Bernoulli matrix with small number of measurements)

% This one is ONLY for those indices which are zero. And we don't know where
% will its dx point in next steps, so after we calculate dx and its in opposite
% direction to z_x, we will have to remove that index from the support.
xk_1 = x_k+delta*del_x_vec;
xk_1(out_x) = 0;
wrong_sign = find(sign(xk_1(gamma_x)).*z_x(gamma_x)==-1);
if ~isempty(gamma_x(wrong_sign))
    disp('Sign mismatch!'); 
    % keyboard;
    
    chk_x = 1;
    delta = 0;
    % can also choose specific element which became non-zero first but all
    % that matters here is AtA(gx,gl) doesn't become singular.
    % [val_wrong_x ind_wrong_x] =  sort(abs(del_x_vec(gamma_x(wrong_sign))),'descend');
    out_x = gamma_x(wrong_sign(1));
end

% If more than one primal constraints became active in previous iteration i.e.,
% more than one elements wanted to enter the support and we added only one.
% So here we need to check if those remaining elements are still active.
i_delta_temp = gamma_lc(abs(pk(gamma_lc)+delta*dk(gamma_lc))-(epsilon-delta) >= 10*eps);
if ~isempty(i_delta_temp)
    if ~isempty(out_lambda)
        i_delta_more = i_delta_temp;%(find(i_delta_temp~=out_lambda));
    else
        i_delta_more = i_delta_temp;
    end
    if length(i_delta_more)>=1 & ~sum((i_delta_temp==i_delta))
        
        disp('Degenerate indices in previous iteration!'); 
        % keyboard;
        
        % ideal way would be to check that incoming element doesn't make AtA
        % singular!
        [v_temp i_temp] = max(-pk(i_delta_more)./dk(i_delta_more));
        i_delta = i_delta_more(i_temp);
        delta = 0;
        chk_x = 0;
        out_x = [];
    end
end