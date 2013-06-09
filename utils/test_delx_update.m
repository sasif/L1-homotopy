clear;

N = 1000;
S = 500;

rseed = sum(100*clock);
% rseed = 0;
% rand('state',rseed);
% randn('state',rseed);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));

A = (randn(N,S));
AtA = A'*A;
z = sign(randn(S,1));

% R-CGS initialize
% Ri = zeros(S,S);
ata = A(:,1)'*A(:,1);
Ri(1,1) = 1/sqrt(ata);
delx_cgs = z(1,1)/ata;
t_cgs = 0;

% MIL initialize
iAtA = zeros(1,1);
%iAtA = sparse(iAtA);
iAtA(1,1) = 1/ata;
delx_mil = z(1,1)/ata;
t_mil = 0;

% MGS qr update initialize
Q = zeros(N,S); R = zeros(S,S);
Q(:,1) = A(:,1)/norm(A(:,1)); R(1,1) = norm(A(:,1));
delx_qr = z(1,1)/ata;
t_qr = 0;

% qrupdate matlab initialize
[Q0 R0] = qr(A(:,1)); % Q0 = 1; R0 = ata;
delx_qrM = z(1,1)/ata;
t_qrM = 0;

% Run the active ones
qr_cgs = 0;
mil = 1;
qr_gs = 1; 
qrM = 0;
error = 0;

for ii = 2:S
    %% True direction
    if error
        AtA = (A(:,1:ii)'*A(:,1:ii));
        delx = AtA\z(1:ii);
    end
    
    %% R factor CGS
    if qr_cgs
        tic; t = cputime;
        ai = A(:,ii);
        Ata = A(:,1:ii-1)'*ai;
        % Pai = Ri(1:ii-1,1:ii-1)*(Ri(1:ii-1,1:ii-1)'*Ata);
        Pai = Ri*(Ri'*Ata);
        Si = ai'*ai-Ata'*Pai;
        at = [-Pai; 1]/sqrt(Si);
        % Ri(1:ii,ii) = at;
        Ri = [[Ri; zeros(1,ii-1)] at];
        alpha = at'*z(1:ii);
        delx_cgs = [delx_cgs;0]+ alpha*at;
        t1 = cputime-t;
        t1 = toc;
        t_cgs = t_cgs+t1;
    end
    %     RiRit = Ri(1:ii,1:ii)*Ri(1:ii,1:ii)';
    %     fprintf('iter=%d, ||delx-delx-cg||_2=%3.4g, ||iAtA-RiRiT||_2=%3.4g.\n',ii, norm(delx-delx_cgs), norm(iAtAgx-RiRit));
    
    %% MIL update
    if mil
        tic; t = cputime;
        ai = A(:,ii);
        Ata = A(:,1:ii-1)'*ai; 
        Pai = iAtA*Ata;
        Si = ai'*ai-Ata'*Pai;
        at = [-Pai; 1]/sqrt(Si);
        alpha = at'*z(1:ii);
        delx_mil = [delx_mil;0]+ alpha*at;
        iAtA = [iAtA zeros(ii-1,1); zeros(1,ii)] + at*at';
        
        % Ata = zeros(S,1);
        % Gamma = 1:ii-1;
        % Ata(Gamma) = A(:,1:ii-1)'*ai;
        % at = -Pai/sqrt(Si);
        % at(ii) = 1/sqrt(Si);
        % at = sparse(at);
        % alpha = at'*z;
        % delx_mil = delx_mil+ alpha*at;
        % iAtA  = iAtA+at*at';
        
        t1 = cputime-t;
        t1 = toc;
        t_mil = t_mil+t1;
    end
    
    %% qr update with MGS for appending a column. 
    if qr_gs
        % compute QR using Gram-Schmidt
        tic; t = cputime;
        ai = A(:,ii);
        v = ai;
        for i=1:ii-1
            R(i,ii) = Q(:,i)'*v; % change ai to v for MGS...
            v = v - R(i,ii)*Q(:,i); % subtract the projection (q_j'a_j)q_j = (q_j'v)q_j!
        end
        % pj = (Q(:,1:ii-1)'*ai);
        % R(1:ii-1,ii) = pj;
        % w = ai-Q(:,1:ii-1)*pj;
        R(ii,ii) = norm(v);
        Q(:,ii) = v/R(ii,ii);
         % Rc = R(1:ii,1:ii);
        % Rct = Rc';
        % delx_qr = Rc\(Rct\z(1:ii));
        ri = R(1:ii,1:ii)\[zeros(ii-1,1); 1];
        delx_qr = [delx_qr; 0]+(ri'*z(1:ii))*ri;
        
        % Can we store inverse of Rc too? If yes, then we can update the
        % direction instead of solving the whole problem. 
        t1 = cputime-t;
        t1 = toc;
        t_qr = t_qr+t1;
    end
    
    %% matlab qr update
    if qrM
        tic; t = cputime;
%         Q0t = [Q0 zeros(ii-1,1); zeros(1,ii-1) 1];
%         R0t = [R0 zeros(ii-1,1); zeros(1,ii-1) 1];
%         loc_vector = zeros(ii,1);
%         loc_vector(end) = 1;
%         ai = A(:,ii);
%         Ata = A(:,1:ii-1)'*ai;
%         rep_row = [Ata; ai'*ai]; % adds the new row to A'A
%         [Q2t R2t] = qrupdate(Q0t, R0t, loc_vector, rep_row);
%         rep_vec = [Ata; -1]; % adds remaining part of the column in A'A and removes the 1 introduced in the beginning.
%         [Q0 R0] = qrupdate(Q2t, R2t, rep_vec, loc_vector);
%         delx_qrM = R0\(Q0'*z(1:ii));
        
        % [Qt, Rt] = qrinsert(Q0, R0, ii, Ata,'col');
        [Q0, R0] = qrinsert(Q0,R0, ii, A(:,ii),'col');
        delx_qrM = (R0\(R0'\z(1:ii)));
        t1 = cputime-t;
        t1 = toc;
        t_qrM = t_qrM+t1;
    end
    
end
if error
    fprintf('iter=%d, t_cgs=%3.4g, t_mil=%3.4g, t_qr=%3.4g, t_qrM=%3.4g, errors: cgs=%3.4g, mil=%3.4g, qr=%3.4g, qrM=%3.4g.\n',ii, t_cgs, t_mil, t_qr, t_qrM, norm(delx-delx_cgs), norm(delx-delx_mil),norm(delx-delx_qr), norm(delx-delx_qrM));
else
    fprintf('iter=%d, t_cgs=%3.4g, t_mil=%3.4g, t_qr=%3.4g, t_qrM=%3.4g.\n',ii, t_cgs, t_mil, t_qr, t_qrM);
end

