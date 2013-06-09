function out = compute_delx(in,A)
% This function solves a system of the form A^TA delx = z.

persistent count_factor count_rec

% Recompute delx directly after max_rec iterations (otherwise update delx
% recursively)
if ~isfield(in, 'max_rec'), max_rec = 100; else max_rec = in.max_rec; end

% Recompute factorization/inverse directly after max_factor iterations
if ~isfield(in,'max_factor'), max_factor = 500; else max_factor = in.max_factor; end

if isempty(count_factor), count_factor = 0; end
if isempty(count_rec), count_rec = 0; end
count_rec = count_rec+1;
count_factor = count_factor+1;

Gamma = in.Gamma;   % Old signal support
% add: no change
% remove: outgoing index moved to the end
rhs = in.rhs; % Full dimensional vector
delx = in.delx;
T = length(Gamma);
M = in.M;
itr = in.itr;

switch in.type
    case 'mil'
        %% Matrix inversion lemma update
        iAtA = in.iAtA;
        
        switch in.add_rem
            case 1
                % add a column
                new_x = in.new_x;
                
                ai = A(:,new_x);
                % Ata = A(:,Gamma)'*ai;
                % Ata = mvprod(A,ai,Gamma,1);
                if issparse(A)
                    Ata = A(:,Gamma)'*ai;
                else
                    Ata = mvprod(A,ai,Gamma,1); % elements of column k in A_Gamma'A_Gamma matrix
                end                
                % Pai = iAtA(1:T,1:T)*Ata;
                Pai = iAtA*Ata;
                
                Si = ai'*ai-Ata'*Pai;
                
                at = [-Pai; 1]/sqrt(Si);
                alpha = at'*rhs([Gamma; new_x]);
                if mod(itr, max_rec)
                    delx = [delx;0]+ alpha*at;
                end
                
                % iAtA(1:ii,1:ii) = iAtA(1:ii,1:ii)+at*at';
                % at = sparse([at; zeros(M-T-1,1)]);
                % iAtA = iAtA + at*at';
                
                iAtA(T+1,T+1) = 0;
                iAtA = iAtA+at*at';
                
                % iAtA = [iAtA zeros(T,1); zeros(1,T+1)] + at*at';
                
                Gamma = [Gamma; new_x]; T = T+1;                                
                
            case 0
                % Remove a column
                new_order = in.new_order; % New order for indices after putting outgoing index at the end
                out_i = in.out_i;
                
                % With full matrix (any order)
                iAtA(1:T,1:T) = iAtA(new_order,new_order);
                
                %         % With sparse matrix (two rows/columns swapped)
                %         tmp = iAtA(:,out_i);
                %         iAtA(:,out_i) = iAtA(:,T);
                %         iAtA(:,T) = tmp;
                %         tmp = iAtA(out_i,:);
                %         iAtA(out_i,:) = iAtA(T,:);
                %         iAtA(T,:) = tmp;
                
                
                Si = 1/iAtA(T,T);
                at = [iAtA(1:T-1,T)*sqrt(Si); 1/sqrt(Si)];
                ab = at(1:T-1);
                alpha = at'*rhs(Gamma);
                if mod(itr, max_rec)
                    delx = delx(new_order);
                    delx = delx(1:T-1)-alpha*ab;
                end
                
                %         iAtA(1:T-1,1:T-1) = iAtA(1:T-1,1:T-1)-ab*ab';
                %         iAtA(T:end,:) = 0;
                %         iAtA(:,T:end) = 0;
                % at = sparse([at; zeros(M-T,1)]);
                iAtA = iAtA - at*at';
                
                iAtA(:,T) = []; iAtA(T,:) = [];
                
                Gamma = Gamma(1:end-1); T = T-1;                
                
            case 2
                % swap a column (standard rank-1 update)
                
                % First remove outgoing column then add incoming column
                new_order = in.new_order; % New order for indices after putting outgoing index at the end
                out_i = in.out_i;
                new_x = in.new_x;
                Gamma_new = Gamma;
                Gamma_new(end) = new_x;
                
                % With full matrix (any order)
                iAtA(1:T,1:T) = iAtA(new_order,new_order);
                
                % remove a column
                Si = 1/iAtA(T,T);
                at = [iAtA(1:T-1,T)*sqrt(Si); 1/sqrt(Si)];
                alpha = at'*rhs(Gamma);
                % delx = delx(new_order);
                % delx  = delx-alpha*at;
                % delx(end) = 0;
                
                % at = sparse([at; zeros(M-T,1)]);
                iAtA = iAtA - at*at';
                iAtA(T,:) = []; iAtA(:,T) = [];
                
                Gamma = Gamma(1:end-1);
                T = T-1;
                
                % add column
                ai = A(:,new_x);
                % Ata = A(:,Gamma)'*ai;
                % Ata = mvprod(A,ai,Gamma,1);
                if issparse(A)
                    Ata = A(:,Gamma)'*ai;
                else
                    Ata = mvprod(A,ai,Gamma,1); % elements of column k in A_Gamma'A_Gamma matrix
                end                
                
                Pai = iAtA*Ata;
                Si = ai'*ai-Ata'*Pai;
                at = [-Pai; 1]/sqrt(Si);
                alpha = at'*rhs([Gamma; new_x]);
                % delx = delx+ alpha*at;
                
                % at = sparse([at; zeros(M-T-1,1)]);
                % iAtA = iAtA + at*at';
                iAtA = [iAtA zeros(T,1); zeros(1,T+1)] + at*at';
                
                Gamma = [Gamma; new_x]; T = T+1;
                delx = iAtA*rhs(Gamma);
        end
        if mod(itr, max_factor) == 0
            iAtA = pinv(A(:,Gamma)'*A(:,Gamma));
        end
        if mod(itr, max_rec) == 0;
            delx = iAtA(1:T,1:T)*rhs(Gamma);
        end
        if abs(Si) < 1e-12
            disp('Gram matrix severely ill-conditioned. Computing pseudo inverse!');
            iAtA = pinv(A(:,Gamma)'*A(:,Gamma));
            delx = iAtA(1:T,1:T)*rhs(Gamma);
        end
        out.iAtA = iAtA;
        out.delx = delx;
        
    case 'qr'
        %% qr update (my code)
        % for when we are sensitive to the numerical sensitivity (I mean stability :p)
        % max_rec controls tradeoff between speed and accuracy
        % Set max_rec = 1, for the best accuracy
        % Set max_rec = inf, for the best speed
        
        Q = in.Q;
        R = in.R;
        M = size(Q,1);
        
        switch in.add_rem
            case 1
                % add a column using modified Gram-Schmidt process
                new_x = in.new_x;
                if isfield(in,'nre')
                    % Number of reorthogonalization steps
                    NRE = in.nre;
                else
                    NRE = 1;
                end
                
                % Q = [Q zeros(size(Q,1),1)];
                R = [R zeros(T,1); zeros(1,T) 0];
                ai = A(:,new_x);
                v = ai;
                
                % See reorth comments in Stephen Becker's SVT code for
                % reorthogonalization of Gram-Schmidt.
                alpha = 0.5; normr = norm(v); normr_old = 0; nre = 0;
                while normr < alpha*normr_old || nre == 0
                    if nre == 0
                        for i=1:T
                            R(i,T+1) = Q(:,i)'*v; % change ai to v for MGS...
                            v = v - R(i,T+1)*Q(:,i); % subtract the projection (q_j'a_j)q_j = (q_j'v)q_j!
                        end
                        R(T+1,T+1) = norm(v);
                    else
                        for i=1:T
                            t = Q(:,i)'*v; % change ai to v for MGS...
                            v = v - t*Q(:,i); % subtract the projection (q_j'a_j)q_j = (q_j'v)q_j!
                        end
                    end
                    normr_old = normr;
                    normr = norm(v);
                    nre = nre + 1;
                    if nre >= NRE
                        break;
                    end
                end
                
                % pj = (Q(:,1:ii-1)'*ai);
                % R(1:ii-1,ii) = pj;
                % w = ai-Q(:,1:ii-1)*pj;
                
                v = v/norm(v);
                Q(:,T+1) = v;
                % Rc = R(1:ii,1:ii);
                % Rct = Rc';
                % delx_qr = Rc\(Rct\z(1:ii));
                if mod(itr, max_rec)
                    ri = R\[zeros(T,1); 1];
                    delx = [delx; 0]+(ri'*rhs([Gamma; new_x]))*ri;
                end
                Gamma = [Gamma; new_x];
            case 0
                % Remove a column using Givens rotations (slower than MGS, but
                % there is no other way... or is it? please tell me if you know one)
                
                new_order = in.new_order; % New order for indices after putting outgoing index at the end
                out_i = in.out_i;
                
                R = R(:,new_order);
                
                for j = out_i:T-1
                    i = j+1;
                    [c, s] = givens(R(j,j), R(i,j));
                    % R([j i],j:end) = [c s; -s c]'*R([j i],j:end);
                    % Qt([j i],:) = [c s; -s c]'*Qt([j i],:);
                    R([j i],j:end) = [c -s; s c]*R([j i],j:end);
                    % Qt([j i],:) = [c -s; s c]*Qt([j i],:);
                    Q(:,[j i]) = Q(:,[j i])*[c s; -s c];
                end
                R = triu(R);
                
                if mod(itr, max_rec)
                    ri = R\[zeros(T-1,1); 1];
                    delx = delx(new_order);
                    delx = delx-(ri'*rhs(Gamma))*ri;
                    delx = delx(1:T-1);
                end
                
                R(:,T) = [];
                R(T,:) = [];
                % R = R(1:T-1,1:T-1);
                % Q = Q(:,1:T-1);
                Q(:,T) = [];
                
                Gamma = Gamma(1:end-1);
                
            case 2
                % swap two columns (standard rank-1 update)
                new_order = in.new_order; % New order for indices after putting outgoing index at the end
                new_x = in.new_x;
                out_i = in.out_i;
                
                R = R(:,new_order);
                for j = out_i:T-1
                    i = j+1;
                    [c, s] = givens(R(j,j), R(i,j));
                    % R([j i],j:end) = [c s; -s c]'*R([j i],j:end);
                    % Qt([j i],:) = [c s; -s c]'*Qt([j i],:);
                    R([j i],j:end) = [c -s; s c]*R([j i],j:end);
                    % Qt([j i],:) = [c -s; s c]*Qt([j i],:);
                    Q(:,[j i]) = Q(:,[j i])*[c s; -s c];
                end
                
                R = triu(R);
                
                %             if mod(itr, max_rec)
                %                 ri = R\[zeros(T-1,1); 1];
                %                 delx = delx(new_order);
                %                 delx_rec = delx-(ri'*rhs(Gamma))*ri;
                %             end
                
                % ONLY APPLICABLE IF A_GAMMA IS A SQUARE MATRIX
                % ri_old = ri; R_old = R;
                R(:,end) = Q'*A(:,new_x);
                Gamma(end) = new_x;
                %             if mod(itr, max_rec)
                %                 ri = R\[zeros(T-1,1); 1];
                %                 delx = delx_rec+(ri'*rhs(Gamma))*ri;
                %             end
                
                % figure(1); plot(delx_orig-delx); shg;
                
                Rt = R';
                delx = R\(Rt\rhs(Gamma));
        end
        if mod(itr, max_factor) == 0
            [Q R] = qr(A(:,Gamma),0);
        end
        if mod(itr, max_rec) == 0
            Rt = R';
            delx = R\(Rt\rhs(Gamma));
        end
        out.Q = Q;
        out.R = R;
        out.delx = delx;
        
    case 'qrM'
        %% qr update (MATLAB)
        Q0 = in.Q0;
        R0 = in.R0;
        
        switch in.add_rem
            case 1
                % add a column
                new_x = in.new_x;
                ai = A(:,new_x);
                [Q0, R0] = qrinsert(Q0,R0, T+1, ai,'col');
                R0t = R0';
                Gamma = [Gamma; new_x];
            case 0
                % Remove a column
                out_i = in.out_i;
                
                [Q0, R0] = qrdelete(Q0,R0, out_i, 'col');
                R0t = R0';
                Gamma = Gamma(1:end-1);
            case 2
                % swap two columns (rank-1 update)
                out_i = in.out_i;
                [Q0, R0] = qrdelete(Q0,R0, out_i, 'col');
                
                new_x = in.new_x;
                ai = A(:,new_x);
                [Q0, R0] = qrinsert(Q0,R0, T, ai,'col');
                R0t = R0';
                
                Gamma(end) = new_x;
        end
        delx = R0\(R0t\rhs(Gamma));
        out.Q0 = Q0;
        out.R0 = R0;
        out.delx = delx;
        
    case 'chol'
        %% Cholesky update
       
        % Fast Cholesky insert and remove functions
        % Adopted from the PMTK toolbox
        % https://code.google.com/p/pmtksupport/source/browse/trunk/lars/lars.m
        %
        % PMTKauthor Karl Skoglund
        % PMTKurl http://www2.imm.dtu.dk/pubdb/views/publication_details.php?id=3897
        % IMM, DTU, kas@imm.dtu.dk

        R = in.R;
        % Updates R in a Cholesky factorization R'R = X'X of a data matrix X. R is
        % the current R matrix to be updated. x is a column vector representing the
        % variable to be added and X is the data matrix containing the currently
        % active variables (not including x).
        %
        % Suppose we have Cholesky factorization for A'A = R'R;  
        % Now we want to update R factors for [A a]'[A a] = [A'A A'a; a'A a'a]
        % To do this write [A a]'[A a] = [R r]'[R r]
        % Now observe that ([R r]')^{-1}[A'a; a'a] = r;
        % Since [R r]' and R' are lower-triangular matrices, their inverses
        % are also lower-triangular. Moreover, the top part of r can be 
        % computed by using inverse of R' alone, as (R')^{-1}A'a. Last 
        % element in r can be computed by looking at norm(r) because 
        % r'r = a'a.
        
        switch  in.add_rem
            case 1
                new_x = in.new_x;
                ai = A(:,new_x);
                diag_k = ai'*ai; % diagonal element k in A'A matrix
                if isempty(R)
                    R = sqrt(diag_k);
                else
                    if issparse(A)
                        Ata = A(:,Gamma)'*ai;
                    else
                        Ata = mvprod(A,ai,Gamma,1); % elements of column k in A_Gamma'A_Gamma matrix
                    end
                    
                    R_k = R'\Ata; % R'R_k = (A'A)_k, solve for R_k
                    R_kk = sqrt(diag_k - R_k'*R_k); % norm(x'x) = norm(R'*R), find last element by exclusion
                    R = [R R_k; [zeros(1,T) R_kk]]; % update R
                end
                
                if mod(itr, max_rec)
                    ri = R\[zeros(T,1); 1];
                    delx = [delx; 0]+(ri'*rhs([Gamma; new_x]))*ri;
                end
                Gamma = [Gamma; new_x];
                
            case 0
                % Deletes a variable from the X'X matrix in a Cholesky factorisation R'R =
                % X'X. Returns the downdated R. This function is just a stripped version of
                % Matlab's qrdelete.
                % function R = choldelete(R,j)
                
                j = in.out_i;                
                R(:,j) = []; % remove column j
                n = size(R,2);
                for k = j:n
                    p = k:k+1;
                    [G,R(p,k)] = planerot(R(p,k)); % remove extra element in column
                    if k < n
                        R(p,k+1:n) = G*R(p,k+1:n); % adjust rest of row
                    end
                end
                R(end,:) = []; % remove zero'ed out row
                Gamma = Gamma(1:end-1);
                
                Rt = R';
                delx = R\(Rt\rhs(Gamma));
                
            case 2
                disp('Not implemented yet');
        end
        if mod(itr, max_factor) == 0
            [Q R] = qr(A(:,Gamma),0);
        end
        if mod(itr, max_rec) == 0
            Rt = R';
            delx = R\(Rt\rhs(Gamma));
        end
        out.R = R;
        out.delx = delx;   
        
    case 'cg'
        %% conjugate gradient iterative update
        W = diag(in.W);
        A_f = @(z) A(:,Gamma)*z;
        At_f = @(z)  A(:,Gamma)'*z;
        AtA_f = @(z) At_f(A_f(z));
        
        b = rhs(Gamma);%-AtA_f([delx;0]); b(1:end-1) = 0;
        x0 = in.x0(Gamma);
        cg_tol = 1e-12; cg_maxit = 500;
        [y,flag,relres,iter,resvec] = pcg(AtA_f,b,cg_tol,cg_maxit,W,[],x0);
        
        out.delx = y;
    otherwise
        disp('NOA');
end
