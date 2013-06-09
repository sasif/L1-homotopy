% betaedge.m
%
% Transition function for a lapped orthogonal transform.
%
% beta(t) = 0   t < 0
%         = sqrt(35t^4-84t^5+70t^6-20t^7)   0<t<1
%         = 1   t > 1
%
% (if -1<t<1, then beta(t) needs to be shifted by 1 and scaled by 1/2, i.e, t --> (t+1)/2)

function b = betaedge(t)

b = zeros(length(t),1);
ti = find((t>=0)&(t<=1));
b(ti) = sqrt(35*t(ti).^4 - 84*t(ti).^5 + 70*t(ti).^6 - 20*t(ti).^7);
b(t>1) = 1;

