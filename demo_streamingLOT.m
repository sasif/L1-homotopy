% demo_streamingLOT
%
% Solves the following streaming problem
% min_x  \|x\|_1 + 1/2*||A  x - y||_2^2
%
% by adding and removing one set of measurements at every time instant
%
% Applications:
%           streaming signal recovery using lapped orthogonal transform (LOT)
%
% We can also add any other regularization operator in the reconstruction
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: August 2012

clear
close all force

% Limit the number of computational threads (for profiling)
% maxNumCompThreads(1);

%% Setup path
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath utils/
addpath utils/utils_Wavelet
addpath utils/utils_LOT
addpath solvers/

fprintf(['----------',datestr(now),'-------%s------------------\n'],mname)

% load RandomStates
%
rseed = 2013;
rand('state',rseed);
randn('state',rseed);

% simulation parameters
mType = 'sign'; % {'randn','orth','rdct','streaming','subsample'};
sType = 'LinChirp'; % % menu('Choose function:', 'LinChirp', 'TwoChirp', 'QuadChirp', 'MishMash','HypChirps','LinChirps', 'Chirps');
SNR = 35;       % additive Gaussian noise

N = 256;   % signal length
R = 8;      % compression rate
M = round(N/R);    % no. of measurements

% signal length
sig_length = 2^15; % 128*128;

% streaming window
P = 5; % size of the working window is P*N

% signal extension for prediction of new coefficients.
eType = 'sym';

% add snapshots of the signal in streaming window and average them before comitting to the output.
avg_output = 0;

% rank-1 update mode
delx_mode = 'mil'; % mil or qr

verbose = 0;

%% LOT setup
% Length of each window is L = N. (extend to adaptive/dyadic windows later?)
eta = 1/2*N; % overlapping factor on the left edge (2*eta samples are shared by adjacent windows)
eta1 = eta; % overlapping factor on the right edge
lp = N; % length of the first interval

in_Psi = []; in_Psi.L = lp; in_Psi.eta = eta; in_Psi.eta1 = eta1;
Psi = create_LOT(in_Psi); % LOT synthesis matrix over a window
% Psi_right = create_LOT_right(N,eta);
in_Psi.eta1 = 0;
Psi_right = create_LOT(in_Psi);

ETA = [eta; eta1];
Lp = [lp];
% Need to control them in a better way
% Right now I divide signal into equal intervals and assign them same
% transition parameters.

%% SM: Sampling modes
% 'universal' sampling scheme
% (align before the overlapping regions of ?2 LOT windows that are measured)
%
SM = 3; % Other sampling modes removed from this release...
SM3_overlap = 0; % overlapping or block diagonal system for measurements!
if SM3_overlap == 0
    % disjoint measurement windows
    LM = N;
    LeftEdge_trunc = 0;
else
    % LM: length of measurement window (e.g., LM = 2*N-2*eta)
    % overlapping measurement windows
    LM = 2*N;
    LeftEdge_trunc = 0;
end
LeftProj_cancel = 1;

fprintf('Streaming signal measureed and reconstructed with overlapping windows using LOT ... \n');
% str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,R) = %d, %d, %d, SM = %d, P = %d, eta_frac = %4.3g, LM = %d, signal length = %d.', mType, sType, SNR, N, M, R, SM, P,eta/N, LM, sig_length);
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,R) = %d, %d, %d, P = %d, SM = %d, LM = %d, \n eta_frac = %4.3g, eType = %s; specified signal-length = %d.', mType, sType, SNR, N, round(N/R), R, P, SM, LM, eta/N, eType, sig_length);
disp(str0);

%% Read signal
in = []; in.type = sType; in.take_fwt = 0; in.Np = 4*N;
% [sigW sig wave_struct] = genSignal(sig_length,in);
sig = genSignal(sig_length,in);
sig = [zeros(2*eta+1,1); sig];
sig_length = length(sig);

% view LOT coefficients...
alpha_vec = apply_LOT(sig,N,eta);
figure(123); imagesc(reshape(alpha_vec,N,length(alpha_vec)/N));
axis xy;
% view spectrogram...
% N = 256; figure(N); spectrogram(sig,N*1,0,N,length(sig),'yaxis');

%% Save results..
streaming_iter = ceil(length(sig)/N);
SIM_stack = cell(streaming_iter,1);
SIM_memory = cell(streaming_iter,1);

x_vec = zeros(N*streaming_iter,1);
xh_vec = zeros(N*streaming_iter,3);
sig_vec = zeros(length(sig),1);
sigh_vec = zeros(length(sig),3);

%% Create analysis/synthesis matrix explicitly
in = [];
in.P = P; in.Psi = Psi;
ETA = eta*ones(P+1,1);
Lp = N*ones(P,1);
eta = ETA(1);
eta1 = ETA(end);
if length(ETA)~=P+1
    error('number of LOT parameters not correct');
end
in.ETA = ETA; in.Lp = Lp;
PSI = create_PSI(in);

%% Setup sensing matrices
in = []; in.type = mType; in.parts = 4;
genAmat_h = @(M,N) genAmat(M,N,in);
in.P = P+1-floor(LM/N);
in.LM = LM; in.M = M; in.N = N;
PHI = create_PHI(in);

%% Overall streaming system matrix
T_length = 2*eta+P*N;
t_ind = 1:T_length;

sigt = sig(t_ind); % Signal under the LOT window at time t
x = PSI'*sigt; % Sparse LOT coefficients

if SM3_overlap == 0
    PSI_M = PSI(1:end-2*eta1,:);
else
    PSI_M = PSI(1:end-2*eta1,:);
end
A = PHI*PSI_M;
if SM3_overlap == 0
    y = PHI*sigt(1:end-2*eta1);
else
    y = PHI*sigt(1:end-2*eta1);
end

leny = length(y);
sigma = sqrt(norm(y)^2/10^(SNR/10)/leny);
e = randn(leny,1)*sigma;
y = y+e;

% parameter selection
% tau = sigma*sqrt(log(N));
tau = max(1e-2*max(abs(A'*y)),sigma*sqrt(log(P*N)));

maxiter = 2*P*N;
err_fun = @(z) (norm(x-z)/norm(x))^2;

%% Initialize by solving a rwt L1 problem
in = [];
in.tau = tau; W = tau;
in.W = W;
in.delx_mode = delx_mode;
in.debias = 0;
in.verbose = 0;
in.plots = 0;
in.record = 1;
in.err_fun = err_fun;
tic
for wt_itr = 1:5
    
    out = l1homotopy(A,y,in);
    xh = out.x_out;
    iter_bpdn = out.iter;
    time_bpdn = toc;
    gamma_bpdn = out.gamma;
    
    % Update weights
    xh_old = xh;
    
    alpha = 1; epsilon = 1;
    beta = M*(norm(xh,2)/norm(xh,1))^2;
    W = tau/alpha./(beta*abs(xh)+epsilon);
    
    W_old = W;
    yh = A*xh;
    Atr = A'*(A*xh-y);
    u =  -W.*sign(xh)-Atr;
    pk_old = Atr+u;
    
    in = out;
    in.xh_old = xh;
    in.pk_old = pk_old;
    in.u = u;
    in.W_old = W_old;
    in.W = W;
end
W = W_old;

sim = 1;
x_vec((sim-1)*N+1:sim*N,1) = x(1:N);
xh_vec((sim-1)*N+1:sim*N,1:3) = [xh(1:N) xh(1:N) xh(1:N)];

s_ind = 1:N+2*eta1;
sig_temp = Psi*xh(1:N);
s_offset = 0;

sig_vec(s_ind) = sigt(s_ind);
sigh_vec(s_ind,1:3) = sigh_vec(s_ind,1:3)+[sig_temp sig_temp sig_temp];

% average instantaneous estimates before committing to output...
estimate_buffer = repmat(xh(1:(P-1)*N,1),1,P-1)/(P-1);

xh_streamingRWT = xh;
x_sparsa = xh;
x_yall1 = xh;

done = 0;
while ~done
    
    
    %% Update the solution after updating the measurement matrix and/or the
    % sparse signal
    y_old = y; x_old = x;
    
    sigt_old = sigt; t_ind_old = t_ind;
    PHI_old = PHI;
    PSI_old = PSI;
    
    % change the LOT parameters for the incoming window?
    % eta or lp ??
    Psi_out = PSI(1:N+ETA(1)+ETA(2),1:N);
    ETA_old = ETA;
    ETA = [ETA(2:end); eta];
    Lp_old = Lp;
    Lp = [Lp(2:end); N];
    eta = ETA(1); eta1 = ETA(end);
    in.ETA = ETA; in.Lp = Lp;
    % PSI = create_PSI(in);
    
    % Shift the sampling window
    t_ind = t_ind+N;
    if t_ind(end) > length(sig)
        break;
    end
    sigt = sig(t_ind);
    
    % System matrix setup...
    % Shift up and left
    PHI(1:end-M,1:end-N) = PHI(M+1:end,N+1:end);
    % new measurement matrix
    Phi = genAmat_h(M,LM);
    PHI(end-M+1:end,end-LM+1:end) = Phi;
    
    A = PHI*PSI_M;
    if SM3_overlap == 0
        y = PHI*sigt(1:end-2*eta);
    else
        y = PHI*sigt(1:end-2*eta1);
    end
    % y(1:M) = y(1:M)-PHI(1:M,1:2*eta)*(Psi(end-2*eta+1:end,:)*xh(1:N));
    
    
    % Sparse LOT coeffs. irrespective of SM
    x = PSI'*sigt;
    
    % shift old measurements and add one new set of measurementts
    e(1:end-M) = e(M+1:end);
    e(end-M+1:end) = randn(M,1)*sigma;
    y= y+e;
    
    A0 = A; x0 = x; y0 = y;
    for solver = {'l1homotopy','sparsa','yall1'}
        solver = char(solver);
        switch solver
            case 'l1homotopy'
                xh = xh_streamingRWT;
            case 'sparsa'
                xh = x_sparsa;
            case 'yall1'
                xh = x_yall1;
        end
        y = y0; A = A0; x = x0;
        xh_old = xh;
        
        % REMOVE the part of outgoing LOT projection in the overlapping region
        % on left side of streaming window...
        % LOT projection of sigt = \Psi_p[O_p] \alpha_p + \Psi_{p-1}[O_p] \alpha_{p-1}
        % where latter is the unwanted part in sigt
        % sigt_proj0 = Psi(end-2*eta+1:end,:)*xh_old(1:N);
        % sigt_proj = sigt;
        % sigt_proj(1:2*eta) = sigt(1:2*eta)-sigt_proj0;
        % ys = y-PHI(:,1:2*eta)*sigt_proj0;
        if LeftProj_cancel
            y(1:M) = y(1:M)-PHI(1:M,1:2*eta)*(Psi_out(end-2*eta+1:end,:)*xh_old(1:N));
        end
        
        %% Update the signal estimate (for warm start)
        xh(1:end-N) = xh(N+1:end);
        % truncate small values in xh
        % xh(abs(xh)<tau/sqrt(log(P*N))) = 0;
        
        %--------------------------------------------%
        % Linear prediction for the incoming samples %
        %--------------------------------------------%
        
        xh1 = xh(end-N+1:end);
        
        if strcmpi(eType,'per')
            % periodic extension
            xh1 = xh(end-N+1:end);
        else
            sig_temp = PSI_M(1:end-N,1:(P-1)*N)*xh(1:end-N);
            sig_temp(1:2*eta) = sig_temp(1:2*eta)+Psi(end-2*eta+1:end,:)*xh_old(1:N);
            % if sim == 1
            %     fprintf('Extension type: %s \n',eType);
            % end
            switch eType
                case 'asym'
                    % anti-symmetric extension
                    sig_temp2 = [sig_temp; -sig_temp(end:-1:1); sig_temp];
                case 'sym'
                    % symmetrci extension
                    sig_temp2 = [sig_temp; sig_temp(end:-1:1); sig_temp];
                otherwise
                    disp('cant do it sire');
            end
            sig_temp2 = sig_temp2(1:P*N+2*eta);
            xh_temp = PSI'*sig_temp2;
            xh1 = xh_temp(end-N+1:end);
        end
        
        gamma1 = find(abs(xh1)>tau);
        [val ind] = sort(abs(xh1),'descend');
        gamma1 = ind(1:min(length(gamma1),ceil(M/log2(N))));
        
        yt = y(end-M+1:end)-A(end-M+1:end,1:end-N)*xh(1:end-N);
        At = A(end-M+1:end,end-N+1:end);
        gamma3 = gamma1;
        %
        % NOTE TODO: Need to replace prediction with rank-update...
        %
        % gamma_all = [find(xh(1:end-N)); gamma3+(P-1)*N];
        % xh_all = 0*xh;
        % xh_all(gamma_all) = A(:,gamma_all)\y; % inv(A(:,gamma_all)'*A(:,gamma_all))*(A(:,gamma_all)'*y-W(gamma_all).*sign(xh(gamma_all)));
        
        xh_t1 = 0*xh(end-N+1:end);
        xh_t1(gamma3) = pinv(At(:,gamma3)'*At(:,gamma3)+2*tau*eye(length(gamma3)))*(At(:,gamma3)'*yt);
        xh_t1(abs(xh_t1)<tau/sqrt(log(N))) = 0;
        xh(end-N+1:end) = xh_t1;
        
        % Filter previous instances?
        %
        % [a b] = lpc(sigh);
        % est_sig = filter([0 -a(2:end)],1,sigh);
        %
        % [val ind] = sort(abs(xh(end-N+1:end)),'descend');
        % xh_t(ind(M+1:end)) = 0;
        %
        %     if P>1
        %         % Learn the variations in the support ???
        %         xh2 = xh(end-2*N+1:end-N);
        %         gamma2 = find(abs(xh2)>0*tau);
        %         [val ind] = sort(abs(xh2),'descend');
        %         gamma2 = ind(1:min(length(gamma2),length(gamma1)));
        %     end
        
        % Truncate small coefficients... ?
        % xh(abs(xh)<tau/sqrt(log(P*N))) = 0;
        fig(111);
        plot([x xh]);
        
        % Remove the top-left edge of the system matrix
        if LeftEdge_trunc
            % fprintf('Consider oldest set of LOT coefficients to be fully known, and remove their contribution from the measurements... \n');
            alpha0h = xh(1:N);
            xh = xh(N+1:end);
            y = y-A(:,1:N)*alpha0h;
            
            A = A(:,N+1:end);
            alpha0 = x(1:N);
            x = x(N+1:end);
        end
        
        % Update weights
        epsilon = 1;
        beta = max(epsilon,M*(norm(xh,2)/norm(xh,1))^1);
        W = tau/alpha./(beta*abs(xh)+epsilon);
        
        W_old = W;
        
        
        if strcmpi(solver,'l1homotopy')
            
            homotopy_mode = 'dummy';
            switch homotopy_mode
                case 'dummy'
                    % create a dummy variable...
                    % use homotopy on the measurements...
                    % in principle, we can start with any xh with this formulation
                    % and any starting value of tau or W...
                    gamma = find(xh);
                    M_trunc = size(A,1); % P*(M-1);
                    % M_trunc = round(P*(N-1)/5);
                    if length(gamma) >= M_trunc
                        disp('length of gamma exceeded number of rows');
                        [xh_sort ind_sort] = sort(abs(xh),'descend');
                        xh(ind_sort(M_trunc+1:end)) = 0;
                        gamma = ind_sort(1:M_trunc);
                    end
                    Atr = A'*(A*xh-y);
                    u =  -W.*sign(xh)-Atr;
                    pk_old = Atr+u;
                otherwise
                    didp('Go back ... no escape');
            end
            
            
            in = out;
            gamma_old = gamma;
            in.gamma = gamma_old;
            switch delx_mode
                case 'mil';
                    % in.delx_mode = 'mil';
                    % The following gram matrix and its inverse can be used from the
                    % previous homotopy. Too lazy to include that right now...
                    % wt BPDN homotopy update
                    AtAgx = A(:,gamma_old)'*A(:,gamma_old);
                    iAtAgx = pinv(AtAgx);
                    in.iAtA = iAtAgx;
                case {'qr','chol'};
                    % in.delx_mode = 'qr';
                    [Q R] = qr(A(:,gamma_old),0);
                    in.Q = Q; in.R = R;
                case 'qrM'
                    % in.delx_mode = 'qrM';
                    [Q0 R0] = qr(A(:,gamma_old));
                    in.Q0 = Q0; in.R0 = R0;
            end
            
            in.xh_old = xh;
            in.pk_old = pk_old;
            in.u = u;
            in.W = W;
            in.delx_mode = delx_mode;
            in.debias = 0;
            in.verbose = 0;
            in.plots = 0;
            in.record = 1;
            in.err_fun = @(z) (norm(x-z)/norm(x))^2;
            tic
            out = l1homotopy(A,y,in);
            time_streamingRWT = toc;
            xh_streamingRWT = out.x_out;
            gamma_streamingRWT = out.gamma;
            iter_streamingRWT = out.iter;
            
            %     fig(333); clf; plot([xh_streamingRWT])
            %     hold on;
            %     plot(xh_old,':k')
            %     plot(setxor(gamma_old,gamma_streamingRWT),0,'r.')
            %     plot(x,'--m')
            
            % Check the solution using BPDN directly with l1homotopy
            %     in = []; x_old = x;
            %     in.W = W;
            %     in.delx_mode = delx_mode;
            %     in.debias = 0;
            %     in.verbose = 0;
            %     in.maxiter = maxiter;
            %     in.plots = 0;
            %     in.record = 1;
            %     in.err_fun = err_fun;
            %     tic
            %     out = l1homotopy(A,y,in);
            %     xh = out.x_out;
            %     iter_bpdn = out.iter;
            %     time_bpdn = toc;
            %     gamma_bpdn = out.gamma;
            %     err_bpdn = out.error_table;
            if LeftEdge_trunc
                xh_streamingRWT = [alpha0h; xh_streamingRWT];
                x = [alpha0; x];
            end
        elseif  strcmpi(solver,'sparsa')
            %% SpaRSA
            x_sparsa = xh; W_sparsa = W/tau; iter_sparsa = 0; time_sparsa = 0;
            if norm(y) > 1e-3
                psi_function = @(x,tau) soft(x,tau*W_sparsa);
                phi_function = @(x) sum(abs(W_sparsa.*x));
                tic;
                [x_sparsa,x_debias_SpaRSA_m,obj_SpaRSA_m_cont,...
                    times_SpaRSA_m_cont,debias_start_SpaRSA_m,mse_SpaRSA_m,taus_SpaRSA_m, numA, numAt]= ...
                    SpaRSA_adpW(y,A,tau,...
                    'Monotone',0,...
                    'adp_wt',0,...
                    'W_new',W_sparsa,...
                    'Debias',0,...
                    'Initialization',x_sparsa,...
                    'StopCriterion',2,...
                    'ToleranceA',1e-4,...
                    'psi',psi_function,...
                    'phi',phi_function,...
                    'Safeguard',1,...
                    'MaxiterA',maxiter,...
                    'Verbose',0,...
                    'True_x',x,...
                    'Continuation',1,...
                    'Continuationsteps',-1);
                
                time_sparsa = toc;
                iter_sparsa = (numA+numAt)/2;
                error_sparsa = norm(x-x_sparsa)/norm(x);
            end
            if LeftEdge_trunc
                x_sparsa = [alpha0h; x_sparsa];
                x = [alpha0; x];
            end
        elseif strcmpi(solver,'yall1')
            %% YALL1
            % set options
            digit = 4; if sigma > 0; digit = 4; end
            opts = [];
            opts.tol = 10^(-digit);
            opts.weights = W/tau;
            opts.print = 0;
            opts.maxit = maxiter;
            opts.x0 = xh;
            opts.nu = 0; opts.rho = tau;
            tic;
            [x_yall1,Out_yall1] = yall1(A,y,opts);
            % time_yall1 = [time_yall1 Out.cputime];
            time_yall1 = toc;
            iter_yall1 = (Out_yall1.cntA+Out_yall1.cntAt)/2;
            err_yall1 = norm(x-x_yall1)/norm(x);
            if LeftEdge_trunc
                x_yall1 = [alpha0h; x_yall1];
                x = [alpha0; x];
            end
        end
    end
    %% Record results
    sim = sim+1;
    SIM_stack{sim} = [sim, tau, ...
        norm(x-xh_streamingRWT)^2/norm(x)^2, sum(iter_streamingRWT,2), sum(time_streamingRWT,2), ...
        norm(x-x_sparsa)^2/norm(x)^2, sum(iter_sparsa,2), sum(time_sparsa,2), ...
        norm(x-x_yall1)^2/norm(x)^2, sum(iter_yall1,2), sum(time_yall1,2)];
    
    % print and plot
    if mod(sim-1,verbose) == 0 && verbose
        fprintf('streaming iter. %d. tau = %3.4g, (err,iter,time): streamingRWT homotopy-%3.4g,%3.4g,%3.4g, SpaRSA-%3.4g,%3.4g,%3.4g, YALL1-%3.4g,%3.4g,%3.4g. \n', ...
            SIM_stack{sim});
    end
    
    %% Plot LOT coeffs. on the window
    fig(1); subplot(211);
    plot([x xh_streamingRWT x_sparsa x_yall1]);
    title('Comparison betweeen the original and reconstructed signal')
    
    %% Reconstructed signal
    xh = xh_streamingRWT;
    x_vec((sim-1)*N+1:sim*N,1) = x(1:N);
    
    % remove the oldest estimate, shift the remaining up and left, and add the new estimate
    estimate_buffer = [[estimate_buffer(N+1:end,2:end); zeros(N,P-2)] xh(1:end-N)/(P-1)];
    if avg_output
        xh_est = xh(1:N);
        xh(1:N) = sum(estimate_buffer(1:N,:),2);
        % fig(123); plot([xh_est xh(1:N) x(1:N)])
        if sim == 2
            disp('output is averaged');
        end
    end
    xh_vec((sim-1)*N+1:sim*N,1) = xh(1:N);
    xh_vec((sim-1)*N+1:sim*N,2) = x_sparsa(1:N);
    xh_vec((sim-1)*N+1:sim*N,3) = x_yall1(1:N);
    
    if SM == 0
        % obsolete
        s_ind = t_ind(1:N+2*eta)-2*eta;
        sig_vec(t_ind(1:N)) = sigt(1:N);
        sigh_vec(s_ind,1) = sigh_vec(s_ind,1)+Psi*xh(1:N);
    else
        s_ind = t_ind(1:N+2*eta);
        sig_vec(s_ind) = sigt(1:N+2*eta);
        sigh_vec(s_ind,1) = sigh_vec(s_ind,1)+Psi*xh(1:N);
        sigh_vec(s_ind,2) = sigh_vec(s_ind,2)+Psi*x_sparsa(1:N);
        sigh_vec(s_ind,3) = sigh_vec(s_ind,3)+Psi*x_yall1(1:N);
    end
    
    
    %% plot recovered signals
    fig(1); subplot(212)
    if sim < P*N+eta*P
        plot([sig_vec(1:s_ind(end)) sigh_vec(1:s_ind(end),1)]);
    else
        plot([sig_vec(s_ind(end)-P*N-2*eta+1:s_ind(end)) sigh_vec(s_ind(end)-P*N-2*eta+1:s_ind(end),1)]);
    end
    
    %% spectrogram
    if mod(sim,25) == -1
        fig(3);
        subplot(211); spectrogram(sig_vec(1:s_ind(end)),N,0,N,length(sig),'yaxis');
        subplot(212); spectrogram(sigh_vec(1:s_ind(end),1),N,N-1,N,length(sig),'yaxis'); shg
    end
    drawnow;
end

mS =  sum(cell2mat(SIM_stack),1);
fprintf('Summed results: streaming_iter %d. tau = %3.4g, \n solver-(err,iter,time): \n streamingRWT homotopy-%3.4g,%3.4g,%3.4g; \n SpaRSA-%3.4g,%3.4g,%3.4g; \n YALL1-%3.4g,%3.4g,%3.4g. \n', streaming_iter, mS(2:end));
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: streaming_iter %d. tau = %3.4g, \n solver-(err,iter,time): \n streamingRWT homotopy-%3.4g,%3.4g,%3.4g; \n SpaRSA-%3.4g,%3.4g,%3.4g; \n YALL1-%3.4g,%3.4g,%3.4g. \n', streaming_iter, mS(2:end));

% l1homotopy-%3.4g,%3.4g,%3.4g;
err_l1homotopy = norm(sig_vec(2*eta+1:N*sim)-sigh_vec(2*eta+1:N*sim,1))^2/norm(sig_vec(2*eta+1:N*sim))^2;
err_sparsa = norm(sig_vec(2*eta+1:N*sim)-sigh_vec(2*eta+1:N*sim,2))^2/norm(sig_vec(2*eta+1:N*sim))^2;
err_yall1 = norm(sig_vec(2*eta+1:N*sim)-sigh_vec(2*eta+1:N*sim,3))^2/norm(sig_vec(2*eta+1:N*sim))^2;
fprintf('Signal MSE: l1homotopy-%3.4g, sparsa-%3.4g, yall1-%3.4g.\n',([err_l1homotopy,err_sparsa,err_yall1]));
fprintf('Signal SER (in dB): l1homotopy-%3.4g, sparsa-%3.4g, yall1-%3.4g.\n',-10*log10([err_l1homotopy,err_sparsa,err_yall1]));


%% plot signal and reconstruction error
x_len = min(length(x_vec),length(xh_vec))-(P-1)*N;
sig_len = min(length(sig_vec),length(sigh_vec))-(P-1)*N-2*eta;
x_vec1 = x_vec(1:x_len);
xh_vec1 = xh_vec(1:x_len,1);
sig_vec1 = sig_vec(1:sig_len);
sigh_vec1 = sigh_vec(1:sig_len,1);

fig(123);
subplot(221);
plot((1:sig_len)/N,sig_vec1, 'LineWidth',1);
axis tight;
title('original signal')
subplot(2,2,2)
plot((1:sig_len)/N,sigh_vec1-sig_vec1, 'LineWidth',1);
axis tight
title('reconstruction error')
subplot(2,2,3);
imagesc(reshape(x_vec1,N,x_len/N)); axis xy;
axis tight;
title('LOT coefficients');
colorbar
subplot(2,2,4);
imagesc(reshape(x_vec1-xh_vec1,N,x_len/N),[0 max(abs(x_vec1))/20]); axis xy
axis tight
title('reconstruction error (LOT)');
colorbar


fig(3);
subplot(211); spectrogram(sig_vec(1:s_ind(end)),N,N-1,N,length(sig),'yaxis');
subplot(212); spectrogram(sigh_vec(1:s_ind(end),1),N,N-1,N,length(sig),'yaxis'); shg