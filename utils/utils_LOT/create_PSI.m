function PSI = create_PSI(in)
% Arrange LOT syntheis matrices over P overlapping windows (edges ignored)

% Input
%   P -- desired number of LOT window
%   Psi -- a fixed matrix or a function handle to create Psi over one LOT window
%   Lp -- length of pth interval 
%           Example: Lp = [N, N, N, ... ]
%   ETA -- list of eta values on the Lp

P = in.P;

if ~isfield(in,'Lp');
    Psi = in.Psi; 
    L = size(Psi,2); 
    eta = (size(Psi,1)-L)/(2);
    
    % time at the edge of the matrices...
    %   tle -- time at the left edge of the first window
    %   tre -- time at the right edge of the last window
    % if isfield(in,'tle'); tle = in.tle; else tle = -eta; end
    % if isfield(in,'tre'); tre = in.tre; else tre = L+eta; end
    % 
    % Psi_1 = Psi; Psi_e = Psi;
    % if tle > -eta; Psi_1 = Psi(tle+eta+1:end,:); end
    % if tre < L+eta; Psi_e = Psi(1:tre,:); end
    % Pr = P*L+2*eta-(tle+eta)-(L+eta-tre);
    
    Pr = P*L+2*eta;
    PSI = zeros(Pr,P*L);
    for p = 1:P
        PSI((p-1)*L+1:p*L+2*eta,(p-1)*L+1:p*L) = Psi;
        % if p == 1
        %     PSI(1:size(Psi_1,1),(p-1)*L+1:p*L) = Psi_1;
        % elseif p==P
        %     PSI(end-size(Psi_e,1)+1:end,(p-1)*L+1:p*L) = Psi_e;
    end
else
    Lp = in.Lp; 
    ETA = in.ETA; 
    if length(Lp)+1 ~= length(ETA)
        error('number of LOT parameters is not correct');
    end
    Pr = sum(Lp)+ETA(1)+ETA(end);
    Pc = sum(Lp);
    PSI = zeros(Pr,Pc);
    pr = 1; pc = 1;
    for p = 1:length(Lp)
        L = Lp(p);
        if L < ETA(p)+ETA(p+1)
            error('LOT interval partition parameters are not workable');
        end
        in_Psi = []; in_Psi.L = L; in_Psi.eta = ETA(p); in_Psi.eta1 = ETA(p+1);
        Psi = create_LOT(in_Psi); % LOT synthesis matrix over a window
        PSI(pr:pr+L+ETA(p)+ETA(p+1)-1,pc:pc+L-1) = Psi;
        pc = pc+L;
        pr = pr+L+ETA(p)-ETA(p+1);
    end
end 
