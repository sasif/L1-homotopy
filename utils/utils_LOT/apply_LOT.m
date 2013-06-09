function alpha = apply_LOT(x,varargin)

% Input
%   x -- input signal 
%   L -- length of the window
%   eta -- length of overlapping window

N = length(x);

% length of the window
if nargin > 1; L = varargin{1}; else L = 128; end
% rollover/overlap 
if nargin > 2; eta = varargin{2}; else eta = 1/2*L; end

if N < L+2*eta
    error('Length of the signal is smaller than the LOT window.');
end

in_Psi = []; in_Psi.L = L; in_Psi.eta = eta; 
Psi = create_LOT(in_Psi);

alpha = [];
tl = 0;
for p = 1:floor((N-2*eta)/L)
    sigt = x(tl+1:tl+L+2*eta);
    alpha = [alpha; Psi'*sigt];
    tl = tl+L;
    %     fig(10); plot(alpha); drawnow;
end
