function out = compute_delta(in)
%
% computes smallest positive delta such that
% delta1 = min_{Gamma^c} ( {(ak-pk)/(dk-bk)}_+, {(-ak-pk)/(dk+bk)}_+)
% delta2 = min_{Gamma} ( -delx/xk )
% delta = min(delta1,delta2)
% 
% input: 
%   required: pk, dk, ak, bk, x, delx_vec
%   optional: 
%       shrinkage_flag (select stepsize by checking 
%           0 -- only if an element in the active set shrinks to zero (outgoing)
%           1 -- only if an inactive constraint becomes active (incoming)
%           2 -- both shrinking elements in the active set
%               and constraint violations of the inactive set
%       nonneg (apply positivity constraint)
%           0 -- no sign constraint on the solution
%           1 -- consider only those constraints that become negative
%           (i.e., causes signal amplitude to become positive)

% output: 
%   delta   -- stepsize 
%   idelta  -- index of the element that enters or leave the support 
%   flag    -- 0    idelta leaves the support
%              1    idelta enters the support
%

gamma_x = in.gamma; gamma_xc = in.gamma_c;
pk = in.pk(gamma_xc); dk = in.dk(gamma_xc);
if length(in.ak) == 1
    ak = in.ak*ones(length(gamma_xc),1);
else
    ak = in.ak(gamma_xc);
end
if length(in.bk) == 1
    bk = in.bk*ones(length(gamma_xc),1);
else
    bk = in.bk(gamma_xc);
end
delx_vec = in.delx_vec; x_k = in.x;

% select delta flag... determine the type of shrinkage
if isfield(in,'shrinkage_flag'); shrinkage_flag = in.shrinkage_flag; else shrinkage_flag = 2; end
% non-negativity constraint? 
if isfield(in,'nonneg'); nonneg = in.nonneg; else nonneg = 0; end

out = [];
       
% Constraint violation (new element to add)
if shrinkage_flag ~= 0 
    % For different values of regularization parameters
    
    % violating constraints are positive
    delta1_constr = (ak-pk)./(dk-bk);
    % delta1_constr = delta1_constr(gamma_xc);
    delta1_pos_ind = find(delta1_constr>0);
    delta1_pos = delta1_constr(delta1_pos_ind);
    [delta1 i_delta1] = min(delta1_pos);
    if isempty(delta1)
        delta1 = inf;
    end
    
    % violating constraints are negative
    delta2_constr = (-ak-pk)./(dk+bk);
    % delta2_constr = delta2_constr(gamma_xc);
    delta2_pos_ind = find(delta2_constr>0);
    delta2_pos = delta2_constr(delta2_pos_ind);
    [delta2 i_delta2] = min(delta2_pos);
    if isempty(delta2)
        delta2 = inf;
    end
    
    if delta1>delta2 || nonneg 
        % infact if we want non-negative solution, then there is no need to
        % compute delta1 at all.
        delta = delta2;
        idelta = gamma_xc(delta2_pos_ind(i_delta2));
    else
        delta = delta1;
        idelta = gamma_xc(delta1_pos_ind(i_delta1));
    end
    
    out.delta_in = delta; out.idelta_in = idelta;
    if shrinkage_flag == 1
        out.delta = delta; out.idelta = idelta;
        out.flag = 1;
        return;
    end    
end

if shrinkage_flag ~= 1
    delta3_constr = (-x_k(gamma_x)./delx_vec(gamma_x));
    delta3_pos_index = find(delta3_constr>0);
    [delta3 i_delta3] = min(delta3_constr(delta3_pos_index));
    out_x = gamma_x(delta3_pos_index(i_delta3));
    if isempty(delta3)
        delta3 = inf;
    end
    
    out.delta_out = delta3; out.idelta_out = out_x;
    
    if shrinkage_flag == 0
        out.delta = delta3; out.idelta = out_x;
        out.flag = 0;
        return;
    end
end

flag = 1;
if delta3 > 0 && delta3 <= delta
    flag = 0;
    delta = delta3;
    idelta = out_x;
end

out.delta = delta; out.idelta = idelta;
out.flag = flag;

