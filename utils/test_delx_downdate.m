clear;

N = 1000;
S = 300;

rseed = sum(100*clock);
% rseed = 0;
% rand('state',rseed);
% randn('state',rseed);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));

A = (randn(N,S));
AtA = A'*A;
z = sign(randn(S,1));


% MIL initialize
iAtA = inv(A'*A);
delx_mil = iAtA*z;
t_mil = 0;

% Givens-rotation qr update initialize
[Q R] = qr(A,0);
% [Q1 R1] = Givred(A);
% [Q2 R2] = qr_givens(A);
delx_qr = R\(R'\z);
t_qr = 0;

% qrupdate matlab infitialize
[Q0 R0] = qr(A); % Q0 = 1; R0 = ata;
delx_qrM = 0;
t_qrM = 0;

% Run the active ones
mil = 1;
qr_giv = 1;
qrM = 0;
error = 1;

Gamma = [1:S]';
Gamma_new = Gamma;

for ii = S-1:-1:2
    
    Gamma = Gamma_new;
    out_i = floor(rand(1,1)*(ii))+1;
    out_x = Gamma(out_i);
    
    % How to remove the column
    swap_shift = 0;
    if swap_shift
        % Swap last column with outgoing column
        Gamma(out_i) = Gamma(end);
        % Gamma = Gamma(1:end-1);
        new_order = [1:out_i-1, ii+1, out_i+1:ii, out_i]';
    else
        % Shift columns to the left (preferred)
        out_loc = Gamma==out_x;
        Gamma = Gamma(~out_loc);
        Gamma = [Gamma; out_x];
        new_order = [1:out_i-1, out_i+1:ii+1, out_i]';
    end
    Gamma_new = Gamma(1:end-1);
    %% True direction
    if error
        AtA = (A(:,Gamma(1:end-1))'*A(:,Gamma(1:end-1)));
        delx = AtA\z(Gamma(1:end-1));
    end
    
    %% MIL update
    if mil
        tic; t = cputime;
        iAtA = iAtA(new_order,new_order);
        Si = 1/iAtA(ii+1,ii+1);
        at = [iAtA(1:ii,ii+1)*sqrt(Si); 1/sqrt(Si)];
        ab = at(1:ii);
        alpha = at'*z(Gamma);
        delx_mil = delx_mil(new_order);
        delx_mil = delx_mil(1:ii)-alpha*ab;
        iAtA = iAtA(1:ii,1:ii)-ab*ab';
        % delx_mil = iAtA*z(Gamma);
        t1 = cputime-t;
        t1 = toc;
        t_mil = t_mil+t1;
    end
    
    %% qr update with Givens rotation for deleting a column
    if qr_giv
        % Update QR using Givens roation
        R_old = R; Q_old = Q;
        Q = Q';
        tic; t = cputime;
        R = [R(:,new_order)];
        
        % if swap_shift
        %     % Make into Hessenberg and then perform another updating
        %     % (Do not swap the columns, work directly with the Hessenberg form below)
        %     for i = ii+1:-1:out_i
        %         j = i-1;
        %         [c, s] = givens(R(j,out_i), R(i,out_i));
        %         R([j i],out_i:end) = [c s; -s c]'*R([j i],out_i:end);
        %         %Q([j i],:) = [c s; -s c]'*Q([j i],:);
        %     end
        % end
        for j = out_i:ii
            i = j+1;
            [c, s] = givens(R(j,j), R(i,j));
            R([j i],j:end) = [c -s; s c]*R([j i],j:end);
            Q([j i],:) = [c -s; s c]*Q([j i],:);
        end
        ri = R\[zeros(ii,1); 1];
        delx_qr = delx_qr(new_order);
        delx_qr_simple = delx_qr-(ri'*z(Gamma))*ri;
        delx_qr = delx_qr_simple(1:end-1);
        
        R = R(1:ii,1:ii);
        % Rt = R';
        % delx_qr = R\(Rt\z(Gamma));
        
        t1 = cputime-t;
        t1 = toc;
        t_qr = t_qr+t1;
        R = triu(R); Q = Q(1:ii,:); Q = Q';
        
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
        %         delx_qrM = R0\(Q0'*z(1:ii));
        
        %         [Q0 R0] = qrupdate(Q2t, R2t, rep_vec, loc_vector);
        [Q0, R0] = qrdelete(Q0,R0, out_i, 'col');
        R0t = R0';
        delx_qrM = (R0\(R0t\z(Gamma(1:end-1))));
        t1 = cputime-t;
        t1 = toc;
        t_qrM = t_qrM+t1;
    end
    
end

if error
    fprintf('iter=%d, t_mil=%3.4g, t_qr=%3.4g, t_qrM=%3.4g, errors: mil=%3.4g, qr=%3.4g, qrM=%3.4g.\n',ii, t_mil, t_qr, t_qrM, norm(delx-delx_mil), norm(delx-delx_qr), norm(delx-delx_qrM));
else
    fprintf('iter=%d, t_mil=%3.4g, t_qr=%3.4g, t_qrM=%3.4g.\n',ii, t_mil, t_qr, t_qrM);
end