% lotwin.m
%
% Window function for the lapped orthgonal transform
%

function g0 = lotwin(t, eta, varargin)

if nargin > 1
    eta = varargin{1};
else
    eta = 0.25;
end
eta1 = eta; 
if nargin > 2
    eta1 = varargin{1};
end

% length
%L = 1;

% We should eventually generalize to this parameterization
% edgewidth
%eta = 1/4;
% g0 = zeros(size(t));
% ti = find((t>=-1/4)&(t<=1/4));
% g0(ti) = betaedge(2*t(ti) + 1/2);
% ti = find((t>1/4)&(t<=3/4));
% g0(ti) = 1;
% ti = find((t>3/4)&(t<=5/4));
% g0(ti) = betaedge(5/2-2*t(ti));


% General value of eta
g0 = zeros(size(t));
tl = find((t>=-eta)&(t<=eta));
g0(tl) = betaedge((t(tl)+eta)/(2*eta));
tc = find((t>=eta)&(t<=1-eta1));
g0(tc) = 1;
tr = find((t>=1-eta1)&(t<=1+eta1));
g0(tr) = betaedge((-t(tr)+1+eta1)/(2*eta1));

% verify
% figure; plot([g0(tl) g0(tr) g0(tl).^2+g0(tr).^2])