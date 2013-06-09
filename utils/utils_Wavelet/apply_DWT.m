function alpha = apply_DWT(x,varargin)

% Apply DWT on a streaming signal
% Input
%   x -- input signal 
%   N -- length of the window
%   wType -- type of wavelet functions
%   J   -- depth of wavelet scales (finest scale)

Nx = length(x);

% length of the window
if nargin > 1; N = varargin{1}; else N = 128; end
% wavelet function type
if nargin > 2; wType = varargin{2}; else wType = 'daub4'; end
% depth of wavelets
if nargin > 3; J = varargin{3}; else J = 3; end
% type of signal extension
if nargin > 4; sym = varargin{4}; else sym = 0; end

in_Psi = []; in_Psi.N = N; 
in_Psi.wType = wType; in_Psi.J = J; in_Psi.sym = sym;

Psi = create_DWT(in_Psi);
L = size(Psi,1);

alpha = [];
tl = 0;
for p = 1:floor((Nx-2*(L-N))/N)
    sigt = x(tl+1:tl+L);
    if sym == 1 || sym == 2
        alpha = [alpha; inv(Psi'*Psi)*(Psi'*sigt)];
    else
        alpha = [alpha;  (Psi'*sigt)];
    end
    tl = tl+N;
    %     fig(10); plot(alpha); drawnow;
end
