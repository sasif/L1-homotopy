function PSI = create_PSI_DWT(in)
% Arrange DWT syntheis matrices over P overlapping windows (edges ignored)

% Input
%   P -- desired number of interval
%   Psi -- a fixed matrix or a function handle to create Psi over one
%   interval
%   Np -- length of intervals 
%           Example: Np = [N, N, N, ... ]
%   wType -- type of wavelets
%   Jp   -- depth of wavelet scales at pth interval...
%   sym -- type of signal extension at the border

P = in.P;

if ~isfield(in,'Np');
    Psi = in.Psi; 
    [L N] = size(Psi); 
    
    PSI = [];
    for p = 1:P
        PSI((p-1)*N+1:(p-1)*N+L,(p-1)*N+1:p*N) = Psi;
    end
else
    Np = in.Np; 
    Jp = in.Jp; 
    wType = in.wType;
    sym = in.sym;
    if length(Np) ~= length(Jp)
        error('number of wavelet parameters is not correct');
    end
    PSI = [];
    pr = 1; pc = 1;
    for p = 1:length(Np)
        N = Np(p);
        J = Jp(p); 
        in_Psi = []; in_Psi.N = N; in_Psi.wType = wType; in_Psi.J = J; in_Psi.sym = sym;
        Psi = create_DWT(in_Psi); % DWT synthesis matrix over a window
        L = size(Psi,1);
        PSI(pr:pr+L-1,pc:pc+N-1) = Psi;
        pr = pr+N;
        pc = pc+N;
    end
end 
