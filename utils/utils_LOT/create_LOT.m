% Creates the LOT representation matrix

function Psi = create_LOT(in)

% lapped window modulated by DCT-IV 
%
% inputs: 
%   L -- lenght of LOT window
%   eta -- transition width on the left edge (a_p,eta_p)
%
% optional input parameters:
%   eta1 -- transition width on the right edge (a_{p+1}, eta_{p+1}) 
%       default (eta)
%   extend -- use overcomplete samples of DCT bases

L = in.L;
eta = in.eta;

% optional parameters
if isfield(in,'eta1'); eta1 = in.eta1; else eta1 = eta; end
if isfield(in,'extend'); extend = in.extend; else extend = 1; end

D = zeros(eta+L+eta1,L);

% LOT window.*cosine
%t = t(:);
% pt = sqrt(2/L)*lotwin(t/L).*cos(pi*(w+1/2)*t/L);
t = [-eta:L+eta1-1]'+0.5;
gn = lotwin(t/L,eta/L,eta1/L);

if eta1 == 0
    % DCT-1 modulation for the LOT window on the right border
    % modulated by DCT-1 (\S 8.4 Mallat) to avoid a sharp discontinuity
    % DCT-I provides even-symmetric extension on the right border
    for nn = 0:L-1
        % D(:,nn+1) = cos(pi*(nn+1/2)*t/L); % DCT-IV
        D(:,nn+1) = cos(pi*nn*t/L); % DCT-I
    end
    D(:,1) = D(:,1)/sqrt(2);
else
    for nn = 0:L-1
        D(:,nn+1) = cos(pi*(nn+1/2)*t/L); % DCT-IV
        % D(:,nn+1) = cos(pi*nn*t/L); % DCT-I
    end
end
Psi = sqrt(2/L).*diag(gn)*D;

% debug
% figure; imagesc(Psi'*Psi-eye(L))
% figure; imagesc([Psi(1:2*eta*L,:)'*Psi(end-2*eta*L+1:end,:)])