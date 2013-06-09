clear;
addpath qr/

N = 512;
M = 256;
maxiter = 500;

rseed = sum(100*clock);
rseed = 0;
% rand('state',rseed);
% randn('state',rseed);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));

O = randn(M,N); 
W = diag(rand(N,1)); W = diag(ones(N,1));
A = O*W;
AtA = A'*A;
z = (randn(N,1));

T_init = 250;
q = randperm(N)';
Gamma = q(1:T_init); 
delx = (A(:,Gamma)'*A(:,Gamma))\z(Gamma);

% MIL initialize
iAtA = zeros(T_init,T_init);
iAtA(1:T_init,1:T_init) = inv(A(:,Gamma)'*A(:,Gamma));
delx_mil = delx;
t_mil = 0;

% MGS qr update initialize
% Q = zeros(M,M); R = zeros(M,M);
Q = zeros(M,M); 
[Q R] = qr(A(:,Gamma),0);
delx_qr = delx;
t_qr = 0;

% qrupdate matlab initialize
[Q0 R0] = qr(A(:,Gamma)); % Q0 = 1; R0 = ata;
delx_qrM = delx;
t_qrM = 0;

% CG initialize
delx_cg = delx;
t_cg = 0;

% Run with following options
mil = 1;
qr_g = 1;
qrM = 1;
cg = 1;
error = 1;
in = [];
Gamma_new = Gamma;

for itr = 1:maxiter;
    
    Gamma = Gamma_new;
    GammaC = setdiff(1:N,Gamma);
    T = length(Gamma);
    if T == 1
        add_rem = 1;
    elseif T == M
        toss = rand;
        add_rem = (toss>0.5)*2;
    else
        add_rem = rand>0;
    end
    
    switch add_rem
        case 1
            ind_i = floor(rand(1,1)*(N-T))+1;
            new_x = GammaC(ind_i);
            Gamma_new = [Gamma; new_x];
            
            in.new_x = new_x;
            in.add_rem = 1;
        case 0
            out_i = floor(rand(1,1)*(T))+1;
            out_x = Gamma(out_i);
            
            % Shift columns to the left (preferred)
            out_loc = Gamma==out_x;
            Gamma = Gamma(~out_loc);
            Gamma = [Gamma; out_x];
            new_order = [1:out_i-1, out_i+1:T, out_i]';
            
            
            %         % Swap the outgoing column with the last one
            %         Gamma(out_loc) = Gamma(end);
            %         Gamma(end) = out_x;
            %         if out_i == T
            %             new_order = [1:out_i-1, T]';
            %         else
            %             new_order = [1:out_i-1, T, out_i+1:T-1, out_i]';
            %         end
            
            Gamma_new = Gamma(1:end-1);
            
            in.new_order = new_order;
            in.out_i = out_i;
            in.add_rem = 0;
        case 2
            out_i = floor(rand(1,1)*(T))+1;
            out_x = Gamma(out_i);
            
            ind_i = floor(rand(1,1)*(N-T))+1;
            new_x = GammaC(ind_i);
            
            % Shift columns to the left (preferred)
            out_loc = Gamma==out_x;
            Gamma = Gamma(~out_loc);
            Gamma = [Gamma; out_x];
            new_order = [1:out_i-1, out_i+1:T, out_i]';
            
            Gamma_new = [Gamma(1:end-1); new_x];

            in.new_order = new_order;
            in.out_i = out_i;
            in.new_x = new_x;
            in.add_rem = 2;
    end
    
    
    %% True direction
    if error
        AtA = (A(:,Gamma_new)'*A(:,Gamma_new));
        delx = AtA\z(Gamma_new);
    end
    
    in.Gamma = Gamma;
    in.rhs = z;
    in.itr = itr;
    
    %% MIL update
    if mil
        t = tic; % t = cputime;
        
        if in.add_rem == 0
            
        end
        
        in.iAtA = iAtA;
        in.delx = delx_mil;
        in.type = 'mil';
        in.M = M;
        out = compute_delx(in,A);
        delx_mil = out.delx;
        iAtA = out.iAtA;
        
        t1 = toc(t); % t1 = cputime-t;
        t_mil = t_mil+t1;
    end
    
    %% qr update with Givens rotation for deleting a column
    if qr_g
        % Update QR using Givens roation
        t=tic; % t = cputime;
        
        in.R = R;
        in.Q = Q;
        in.nre = 2;
        in.delx = delx_qr;
        in.type = 'qr';
        out = compute_delx(in,A);
        R = out.R;
        Q = out.Q;
        delx_qr = out.delx;
        
        t1 = toc(t); % t1 = cputime-t;
        t_qr = t_qr+t1;
        
        orth_err = Q'*Q-eye(size(Q,2));
        % figure(1); plot(vec(orth_err));  shg
        if norm(orth_err(:)) > 1e-3
            error('reorthogonalize...');
        end
    end
    
    %% matlab qr update
    if qrM
        t=tic; % t = cputime;
        
        in.R0 = R0;
        in.Q0 = Q0;
        in.delx = delx_qrM;
        in.type = 'qrM';
        out = compute_delx(in,A);
        R0 = out.R0;
        Q0 = out.Q0;
        delx_qrM = out.delx;
        
        t1 = toc(t); % t1 = cputime-t;
        t_qrM = t_qrM+t1;
    end
    
    %% CG update
    if cg
        t=tic; % t = cputime;
        
        in.delx = delx_cg;
        in.type = 'cg';
        out = compute_delx(in,A);
        delx_cg = out.delx;
        
        t1 = toc(t); % t1 = cputime-t;
        t_cg= t_cg+t1;
    end
    
    
    if error
        fprintf('iter=%d, T=%d, t_mil=%3.4g, t_qr=%3.4g, t_qrM=%3.4g, t_cg=%3.4g, errors: mil=%3.4g, qr=%3.4g, qrM=%3.4g, cg=%3.4g.\n',itr, T,t_mil, t_qr, t_qrM, t_cg, norm(delx-delx_mil), norm(delx-delx_qr), norm(delx-delx_qrM),norm(delx-delx_cg));
    else
        fprintf('iter=%d, T=%d, t_mil=%3.4g, t_qr=%3.4g, t_qrM=%3.4g, t_cg=%3.4g. \n',itr, T, t_mil, t_qr, t_qrM, t_cg);
    end
    pause(1/60)
end
