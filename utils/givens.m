function [c s] = givens(a,b)
% Givens plane rotation [c s;-s c]. Entries c and s
% are computed using numbers x1 and x2.
if b==0
    c = 1; s=0;
else
    if abs(b)> abs(a)
        tau = -a/b; s=1/sqrt(1+tau^2); c=s*tau;
    else
        tau = -b/a; c=1/sqrt(1+tau^2); s = c*tau;
    end
end
