% demo_KalmanRWT
%
% Solves the following dynamic BPDN problem over a window t = t_1,...,t_L
% min_x \sum_t \|W_t x_t\|_1 + 1/2*||A_t x_t - y_t||_2^2 + 1/2||F_t x_t - x_t+1\|_2^2
%
% which updates the solution as the signal changes according to a linear
% dynamical system.
%
% for instance, y_t = A_t x_t + e_t
%               x_t+1 = F_t x_t + f_t 
%       where F_t is a partially known function that models prediction
%       between the consecutive x_t and f_t denotes the prediction error 
%       (e.g., a random drift)
%
% Applications:
%       streaming signal recovery using a dynamic model
% 
%   	track a signal as y, A, and/or x change... 
%       predict an estimate of the solution and
%       update weights according to the predicted solution
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: November 2012

clear
close all force

% Limit the number of computational threads (for profiling)
maxNumCompThreads(1);

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
mType = 'sign'; % {'randn','orth','rdct'};
mFixed = 0; % measurement system time-(in)variant
sType = 'pcwreg'; % {'heavisine', 'pcwreg', 'blocks','pcwPoly'}
SNR = 35;       % additive Gaussian noise

wt_pred = sqrt(0.5);

N = 256;   % signal length
R = 4; % compression rate
M = round(N/R);    % no. of measurements

LM = 1*N; % LM: length of measurement window
LS_Kalman = 'smooth'; % {'filter','smooth','inst'};

% streaming window
P = 3; % size of the working window is P*N

% signal length
sig_length = 2^15; % 128*128;

% signal dynamics
dType = 'crshift'; % type of dynamics 'crshift', or 'static'
cshift = -1;
rshift_max = 0.5;
rshift_h = @(z) (rand-0.5)*rshift_max*2;

% DWT parameters
% type of scaling function
% depth of scaling functions (number of wavelet levels)
% type of extension for DWT (0 - periodic extension, 3 - streaming)
wType = 'daub79'; sym = 1;
wType = 'daub8'; sym = 0;
J = log2(N)-3;

% rank-1 update mode
delx_mode = 'mil'; % mil or qr

% add snapshots of the signal in streaming window and average them before comitting to the output.
avg_output = 0; 

verbose = 0;


%% SM: Sampling modes
% % LM: length of measurement window
% LM = 2*N; % 'universal' sampling scheme (align before the overlapping regions of DWT windows that are measured)
if LM > N
    LeftEdge_trunc = 1;
else
    LeftEdge_trunc = 0;
end
LeftProj_cancel = 1;


%%
fprintf('CS-Kalman tracking a dynamical signal and reweighting..\n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,R) = %d, %d, %d, P = %d, LM = %d, LS_Kalman-%s \n wType-%s, J = %d, sym = %d, specified signal-length = %d, \n dType-%s, cshift = %d, rshift_max = %0.3g, wt_pred = %0.3g. ', mType, sType, SNR, N, round(N/R), R, P, LM, LS_Kalman, wType, J, sym, sig_length, dType, cshift, rshift_max, wt_pred);
disp(str0);

%% DWT setup
% DWT parameters
% Length of each window is L. (extend to adaptive/dyadic windows later?)
% wType = 'daub4'; % type of scaling function
% J = 3; % depth of scaling functions (number of wavelet levels)
% sym = 3; % type of extension for DWT (0 - periodic extension, 3 - streaming)
in_Psi = []; in_Psi.N = N; in_Psi.J = J; in_Psi.wType = wType; in_Psi.sym = sym;
Psi = create_DWT(in_Psi); % DWT synthesis matrix over a window
L = size(Psi,1);

%% Signal generation

% Setup dynamical model
% At every time instance, add to the original/previous signal
% an integer circular shift that is known
% a random drift that is unknown
if strcmpi(dType, 'crshift')
    % Generate a signal by circular shift and a random drift in a seed
    % signal
    in = []; in.type = sType; in.randgen = 0; in.take_fwt = 0;
    [x_init sig wave_struct] = genSignal(N,in);
    
    F_h = @(x,cshift,rshift) interp1(1:N,circshift(x,cshift),[1:N]+rshift,'linear','extrap')';
    
    F0 = zeros(N);
    for ii = 1:N;
        F0(:,ii) = F_h(circshift([1; zeros(N-1,1)],ii-1),cshift,0);
    end
    sigt = sig; sig = [];
    for ii = 1:round(sig_length/N);
        rshift = rshift_h(1);
        sigt = F_h(sigt, cshift, rshift);
        sig = [sig; sigt];
    end
else
    % Generate a predefined streaming signal
    in = []; in.type = sType; in.randgen = 0; in.take_fwt = 0;
    [x_init sig wave_struct] = genSignal(N,in);
    
    cshift = 0; rshift = 0;
    F_h = @(x,cshift,rshift) x;
    F0 = eye(N);
end
% sig = [zeros(L-N,1);sig];
sig_length = length(sig);

% view DWT coefficients...
alpha_vec = apply_DWT(sig,N,wType,J,sym);
figure(123);
subplot(211); imagesc(reshape(alpha_vec,N,length(alpha_vec)/N));
axis xy;
subplot(212); plot(alpha_vec);

% view innovations in the signal.. 
% dsig = []; for n = 0:N:length(sig)-N; dsig = [dsig; sig(n+1:n+N)-circshift(sig(n+N+1:n+2*N),1)]; figure(1); plot([sig(n+1:n+N) sig(n+1:n+N)-circshift(sig(n+N+1:n+2*N),1)]); pause; end

% Simulation parameters
streaming_iter = ceil(length(sig)/N);
SIM_stack = cell(streaming_iter,1);
SIM_memory = cell(streaming_iter,1);

x_vec = zeros(N*streaming_iter,1);
xh_vec = zeros(N*streaming_iter,3);
sig_vec = zeros(length(sig),1);
sigh_vec = zeros(length(sig),3);

%% Setup sensing matrices
in = []; in.type = mType;
if mFixed
    At = genAmat(M,LM,in);
    genAmat_h = @(m,n) At;
else
    genAmat_h = @(M,N) genAmat(M,N,in);
end
in.P = P-(LM-N)/N;
in.LM = LM; in.M = M; in.N = N;
PHI = create_PHI(in);

%% Dynamics matrix
F = zeros(P*N,(P+1)*N);
for p = 1:P
    F((p-1)*N+1:p*N,(p-1)*N+1:(p+1)*N) = [F0 -eye(N)];
end
F = wt_pred*F(:,N+1:end);

%% Create analysis/synthesis matrix explicitly and compute sparse coeffs.
in = [];
in.P = P; in.Psi = Psi;
% in.P = P; in.Jp = Jp; in.wType = wType; in.N = N; in.sym = sym;
PSI = create_PSI_DWT(in);

% Sparse coefficients...
T_length = size(PSI,1);
t_ind = 1:T_length;

sigt = sig(t_ind); % Signal under the LOT window at time t
if sym == 1 || sym == 2
    x = pinv(PSI'*PSI)*(PSI'*sigt); % Sparse LOT coefficients
else
    x = PSI'*sigt;
end

%% initialize with a predicted value of first x
% xt = x(1:N);
%
% At = genAmat_h(M,N);
% sigma = sqrt(norm(At*xt)^2/10^(SNR/10)/M);
% e = randn(M,1)*sigma;
% yt = At*xt+e;
%
% tau = max(1e-2*max(abs(At'*yt)),sigma*sqrt(log(N)));
%
% % rwt L1 with the first set of measurement...
% in = [];
% in.tau = tau; W = tau;
% in.delx_mode = delx_mode;
% for wt_itr = 1:5
%     W_old = W;
%
%     out = l1homotopy(At,yt,in);
%     xh = out.x_out;
%
%     % Update weights
%     xh_old = xh;
%     alpha = 1; epsilon = 1;
%     beta = M*(norm(xh_old,2)/norm(xh_old,1))^2;
%     W = tau/alpha./(beta*abs(xh_old)+epsilon);
%
%     yh = At*xh_old;
%     Atr = At'*(At*xh-yt);
%     u =  -W.*sign(xh)-Atr;
%     pk_old = Atr+u;
%
%     in = out;
%     in.xh_old = xh;
%     in.pk_old = pk_old;
%     in.u = u;
%     in.W_old = W_old;
%     in.W = W;
% end
% xh(abs(xh)<tau/sqrt(log(N))) = 0;


% Another way to initialize...
% Best M/2-sparse signal...
% [val_sort ind_sort] = sort(abs(x),'descend');
% xh = x;
% xh(ind_sort(P*N/2+1:end)) = 0;

% Oracle value for the initialization
xh = x; disp('oracle initialization');

% model for the outgoing window...
sim = 1;
st_ind = N;
t_ind = st_ind+t_ind;
s_ind = t_ind(1:L);

sig_out = PSI(st_ind+1:st_ind+N,:)*xh;
xh = xh(st_ind+1:end);

xh_out = xh(1:N);
x_vec((sim-1)*N+1:sim*N,1) = x(st_ind+1:st_ind+N);
xh_vec((sim-1)*N+1:sim*N,1:3) = [xh_out xh_out xh_out];

sig_temp = Psi*xh_out;
sig_temp = [sig_out; sig_temp(N+1:end)];
sig_vec(s_ind) = sigt(s_ind);
sigh_vec(s_ind,1:3) = sigh_vec(s_ind,1:3)+[sig_temp sig_temp sig_temp];


%% Generate complete measurement system
% Sparse coefficients...
t_ind = t_ind + N;
sigt = sig(t_ind); % Signal under the LOT window at time t
if sym == 1 || sym == 2
    x = pinv(PSI'*PSI)*(PSI'*sigt); % Sparse LOT coefficients
else
    x = PSI'*sigt;
end

y = PHI*sigt(1:end-(L-N));

leny = length(y);
sigma = sqrt(norm(y)^2/10^(SNR/10)/leny);
e = randn(leny,1)*sigma;
y = y+e;


PSI_M = PSI(1:end-(L-N),:);
A = [PHI; F]*PSI_M;


sig_out = sigh_vec(t_ind(1:N)-N,1);
y = [y; -wt_pred*F0*sig_out; zeros((P-1)*N,1)];

% REMOVE the part of outgoing DWT projection in the overlapping region
% on left side of streaming window...
if LeftProj_cancel
    y = y-[PHI(:,1:(L-N));F(:,1:(L-N))]*(Psi(end-(L-N)+1:end,:)*xh_out(1:N));
end

%% parameter selection 
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

sim = sim+1;
x_vec((sim-1)*N+1:sim*N,1) = x(1:N);
xh_vec((sim-1)*N+1:sim*N,1:3) = [xh(1:N) xh(1:N) xh(1:N)];

s_ind = t_ind(1:L);
sig_temp = Psi*xh(1:N);
sig_vec(s_ind) = sigt(1:L);
sigh_vec(s_ind,1:3) = sigh_vec(s_ind,1:3)+[sig_temp sig_temp sig_temp];

% average instantaneous estimates before committing to output...
estimate_buffer = repmat(xh(1:(P-1)*N,1),1,P-1)/(P-1);

xh_streamingRWT = xh;
x_sparsa = xh;
x_yall1 = xh;


%% Kalman initialization
if LM == N
    Pk_1 = eye(N)/(wt_pred)^2;
    sig_kalman = sig_vec(t_ind(1:N)-N,1);
    
    Ak = PHI(1:M,1:N);
    yk = y(1:M);
    x_k = F0*sig_kalman;
    P_k = F0*Pk_1*F0'+1/(wt_pred^2)*eye(N);
    PAt = P_k*Ak';
    Kk = PAt*(pinv(Ak*PAt+eye(M)));
    Pk_1 = P_k-Kk*PAt';
    sig_kalman = x_k + Kk*(yk-Ak*x_k);
    
    sig_temp = sigh_vec(t_ind(1:N)-N,3);
    y_kalman = y;
    y_kalman(P*M+1:P*M+N) = -wt_pred*(F0*sig_temp);
    
    switch LS_Kalman
        case 'inst'
            sig_P = [PHI;F]\y_kalman;
            sig_kalman = sig_P(1:N);
        case 'smooth'
            % solves for x_1 using the prediction covariance matrix
            % from all previous measurements and smoothing with P-1 future measurements
            % minimize 1/2 (x_1-x_1|0)'*P_1|0(x_1-x_1|0)
            % + \sum_{k = 1,...,P} 1/2||y_k-A_k x_k||_2^2 + lambda/2||F_k
            % x_k-x_k+1||_2^2
            
            iP_k = pinv(P_k);
            Pmat = PHI'*PHI + F'*F;
            Pmat(1:N,1:N) = Pmat(1:N,1:N)-wt_pred^2*eye(N)+iP_k;
            Pty = PHI'*y_kalman(1:M*P)+[iP_k*(F0*sig_temp); zeros((P-1)*N,1)];
            sig_P2 = pinv(Pmat)*Pty;
            sig_kalman = sig_P2(1:N);
        case 'filter'
            % no change...
    end
    
    sigh_vec(t_ind(1:N),3) = sig_kalman;
end


%% GO...

done = 0;
while ~done
    
    % Update the solution after updating the measurement matrix and/or the
    % sparse signal
    x_old = x;
    y_old = y; A_old = A;
    
    sigt_old = sigt; t_ind_old = t_ind;
    PHI_old = PHI;
    
    % Shift the sampling window
    t_ind = t_ind+N;
    if t_ind(end) > length(sig)
        break;
    end
    sigt = sig(t_ind); % Signal under the LOT window at time t
    if sym == 1 || sym == 2
        x = pinv(PSI'*PSI)*(PSI'*sigt); % Sparse LOT coefficients
    else
        x = PSI'*sigt;
    end
    
    % System matrix setup...
    % Shift up and left
    PHI(1:end-M,1:end-N) = PHI(M+1:end,N+1:end);
    % new measurement matrix
    Phi = genAmat_h(M,LM);
    PHI(end-M+1:end,end-LM+1:end) = Phi;
    
    % shift old measurements and add one new set of measurementts
    y = PHI*sigt(1:end-(L-N));
    e(1:end-M) = e(M+1:end);
    e(end-M+1:end) = randn(M,1)*sigma;
    y= y+e;
    
    A = [PHI; F]*PSI_M;
    
    A0 = A; x0 = x; y0 = y;
    for solver = {'l1homotopy','sparsa'}
        solver = char(solver);
        switch solver
            case 'l1homotopy'
                xh = xh_streamingRWT;
                sig_out = sigh_vec(t_ind(1:N)-N,1);
            case 'sparsa'
                xh = x_sparsa;
                sig_out = sigh_vec(t_ind(1:N)-N,2);
            case 'yall1'
                xh = x_yall1;
                sig_out = sigh_vec(t_ind(1:N)-N,3);
        end
        y = y0; A = A0; x = x0;
        xh_old = xh;
        y = [y; -wt_pred*F0*sig_out; zeros((P-1)*N,1)];
        
        % REMOVE the part of outgoing DWT projection in the overlapping region
        % on left side of streaming window...
        if LeftProj_cancel
            y = y-[PHI(:,1:(L-N));F(:,1:(L-N))]*(Psi(end-(L-N)+1:end,:)*xh_old(1:N));
        end
        
        % Update the signal estimate (for warm start)
        xh(1:end-N) = xh(N+1:end);
        sigh_old = PSI(1:end-L,:)*xh;
        % sigh_pred = [sigh_old; F_h(sigh_old(end-N+1:end),cshift,0); zeros(L-N,1)];
        sigh_pred = [sigh_old; F_h(sigh_old(end-N+1:end),cshift,0)];
        
        if sym == 3
            sigh_temp = F_h(sigh_pred(end-N+1:end),cshift,0);
            sigh_pred = [sigh_pred; sigh_temp(1:L-N)];
            % sigh_pred = [sigh_pred; linspace(sigh_pred(end),0,L-N)'];
        end
        if sym == 1 || sym == 2
            xh = pinv(PSI'*PSI)*(PSI'*sigh_pred); % Sparse LOT coefficients
        else
            xh = PSI'*sigh_pred;
        end
        xh(abs(xh)<tau/sqrt(log(P*N))) = 0;
        %         xh_temp = xh(end-N+1:end);
        %         xh_temp(abs(xh_temp)<tau/sqrt(log(P*N))) = 0;
        %         xh(end-N+1:end) = xh_temp;
        if sym == 3 % truncate coefficients for overlapping wavelets... 
            for p = 2.^(0:J)
                xh((P-1)*N+N/p) = 0;                
            end
        end
        
        fig(111);
        plot([x xh]);
        
        % Remove the top-left edge of the system matrix
        if LeftEdge_trunc
            % fprintf('Consider oldest set of LOT coefficients to be fully known, and remove their contribution from the measurements... \n');
            alpha0h = xh(1:N);
            xh = xh(N+1:end);
            y = y-A(:,1:N)*alpha0h;
            A = A(:,N+1:end);
            
            A_old = A; y_old = y;
            A(size(PHI,1)+1:size(PHI,1)+N,:) = [];
            y(size(PHI,1)+1:size(PHI,1)+N) = [];
            
            alpha0 = x(1:N);
            x = x(N+1:end);
        end
        
        % Update weights
        alpha = 1; epsilon = 1;
        beta = M*(norm(xh,2)/norm(xh,1))^1;
        W = tau/alpha./(beta*abs(xh)+epsilon);
        W_old = W;
        
        if strcmpi(solver,'l1homotopy');
            
            homotopy_mode = 'dummy';
            switch homotopy_mode
                case 'dummy'
                    % create a dummy variable...
                    % use homotopy on the measurements...
                    % in principle, we can start with any xh_old with this formulation
                    % and any starting value of tau or W...
                    gamma = find(xh);
                    M_trunc = size(A,1); % P*(M-1);
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
            % Reconstructed signal
            if LeftEdge_trunc
                x = [alpha0; x];
                xh_streamingRWT = [alpha0h; xh_streamingRWT];
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
            % Reconstructed signal
            if LeftEdge_trunc
                x = [alpha0; x];
                x_sparsa = [alpha0h; x_sparsa];
            end
        elseif strcmpi(solver,'yall1')
            
            %% YALL1
            disp('yall1 only works when A is underdetermined');
            % set options
            digit = 4; if sigma > 0; digit = 4; end
            opts = [];
            opts.tol = 10^(-digit);
            opts.weights = W/tau;
            opts.print = 0;
            opts.maxit = maxiter;
            opts.nonorth = 1;
            % opts.x0 = xh;
            opts.nu = 0; opts.rho = tau;
            tic;
            [x_yall1,Out_yall1] = yall1(A,y,opts);
            % time_yall1 = [time_yall1 Out.cputime];
            time_yall1 = toc;
            iter_yall1 = (Out_yall1.cntA+Out_yall1.cntAt)/2;
            err_yall1 = norm(x-x_yall1)/norm(x);
            % Reconstructed signal
            if LeftEdge_trunc
                x = [alpha0; x];
                x_yall1 = [alpha0h; x_yall1];
            end
            
            if max(abs(x_yall1)) > 500
                stp = 1;
            end
        end
    end
    
    %% Plot DWT coeffs. on the window
    fig(1); subplot(211);
    plot([x xh_streamingRWT x_sparsa]);
    title('Comparison betweeen the original and reconstructed signal')
    
    %% Reconstructed signal
    sim = sim+1;
    x_vec((sim-1)*N+1:sim*N,1) = x(1:N);
    xh = xh_streamingRWT;
    
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
    
    s_ind = t_ind(1:L);
    sig_vec(s_ind) = sigt(1:L);
    sigh_vec(s_ind,1) = sigh_vec(s_ind,1)+Psi*xh(1:N);
    sigh_vec(s_ind,2) = sigh_vec(s_ind,2)+Psi*x_sparsa(1:N);
    
    
    
    %% plot recovered signals
    fig(1); subplot(212)
    plot([sig_vec(1:s_ind(end)) sigh_vec(1:s_ind(end),1)]);
    
    drawnow;
    
    %% Kalman recursion
    if LM == N
        
        Ak = PHI(1:M,1:N);
        yk = y(1:M);
        x_k = F0*sig_kalman;
        P_k = F0*Pk_1*F0'+1/(wt_pred^2)*eye(N);
        PAt = P_k*Ak';
        Kk = PAt*(pinv(Ak*PAt+eye(M)));
        sig_kalman = x_k + Kk*(yk-Ak*x_k);
        Pk_1 = P_k-Kk*PAt';
        %         if mod(sim,50) == 0
        %             Pk_1 = eye(N)/wt_pred^2;
        %         end
        
        sig_temp = sigh_vec(t_ind(1:N)-N,3);
        y_kalman = y;
        y_kalman(P*M+1:P*M+N) = -wt_pred*(F0*sig_temp);
        
        switch LS_Kalman
            case 'inst'
                sig_P = [PHI;F]\y_kalman;
                sig_kalman = sig_P(1:N);
            case 'smooth'
                % solves for x_1 using the prediction covariance matrix
                % from all previous measurements and smoothing with P-1 future measurements
                % minimize 1/2 (x_1-x_1|0)'*P_1|0(x_1-x_1|0)
                % + \sum_{k = 1,...,P} 1/2||y_k-A_k x_k||_2^2 + lambda/2||F_k
                % x_k-x_k+1||_2^2
                
                iP_k = pinv(P_k);
                Pmat = PHI'*PHI + F'*F;
                Pmat(1:N,1:N) = Pmat(1:N,1:N)-wt_pred^2*eye(N)+iP_k;
                Pty = PHI'*y_kalman(1:M*P)+[iP_k*(F0*sig_temp); zeros((P-1)*N,1)];
                sig_P2 = pinv(Pmat)*Pty;
                sig_kalman = sig_P2(1:N);
            case 'filter'
                % do nothing
        end
        
        sigh_vec(t_ind(1:N),3) = sig_kalman;
        
        fig(33)
        plot([sig_vec(s_ind(1:N),1) sigh_vec(s_ind(1:N),1) sigh_vec(s_ind(1:N),3)]);
        title('comparison with LS-Kalman');
        err_l1 = norm(sig_vec(s_ind(1:N),1)-sigh_vec(s_ind(1:N),1))^2/norm(sig_vec(s_ind(1:N),1))^2;
        err_ls = norm(sig_vec(s_ind(1:N),1)-sigh_vec(s_ind(1:N),3))^2/norm(sig_vec(s_ind(1:N),1))^2;
        if mod(sim-1,verbose) == 0 && verbose
            fprintf('L1 vs LS Kalman: sim %d  --%3.4g : %3.4g--\n',sim, err_l1, err_ls);
        end
    end
    
    %% Record results
    SIM_stack{sim} = [sim, tau, ...
        norm(x-xh_streamingRWT)^2/norm(x)^2, sum(iter_streamingRWT,2), sum(time_streamingRWT,2), ...
        norm(x-x_sparsa)^2/norm(x)^2, sum(iter_sparsa,2), sum(time_sparsa,2), ...
        err_ls];
    
    % print and plot
    if mod(sim-1,verbose) == 0 && verbose
        fprintf('streaming iter. %d. tau = %3.4g, (err,iter,time): streamingRWT homotopy-%3.4g,%3.4g,%3.4g, SpaRSA-%3.4g,%3.4g,%3.4g, LS-Kalman-%3.4g\n', ...
            SIM_stack{sim});
    end
end
fprintf('\n');
mS =  sum(cell2mat(SIM_stack),1);
fprintf('Summed results: streaming_iter %d. tau = %3.4g, \n solver-(err,iter,time): \n streamingRWT homotopy-%3.4g,%3.4g,%3.4g; \n SpaRSA-%3.4g,%3.4g,%3.4g; \n LS-Kalman-%3.4g. \n', streaming_iter, mS(2:end));
% mS =  mean(cell2mat(SIM_stack),1);
% fprintf('Average results: streaming_iter %d. tau = %3.4g, \n solver-(err,iter,time): \n streamingRWT homotopy-%3.4g,%3.4g,%3.4g; \n SpaRSA-%3.4g,%3.4g,%3.4g; \n LS-Kalman-%3.4g. \n', streaming_iter, mS(2:end));

% l1homotopy-%3.4g,%3.4g,%3.4g;
st_ind = 2*N;
err_l1homotopy = norm(sig_vec(st_ind+1:N*sim)-sigh_vec(st_ind+1:N*sim,1))^2/norm(sig_vec(st_ind+1:N*sim))^2;
err_sparsa = norm(sig_vec(st_ind+1:N*sim)-sigh_vec(st_ind+1:N*sim,2))^2/norm(sig_vec(st_ind+1:N*sim))^2;
if LM == N
    err_kalman = norm(sig_vec(st_ind+1:N*sim)-sigh_vec(st_ind+1:N*sim,3))^2/norm(sig_vec(st_ind+1:N*sim))^2;
    fprintf('Signal MSE: l1homotopy-%3.4g, sparsa-%3.4g, LS-kalman-%3.4g.\n',([err_l1homotopy,err_sparsa,err_kalman]));
    fprintf('Signal SER (in dB): l1homotopy-%3.4g, sparsa-%3.4g, LS-kalman-%3.4g.\n',-10*log10([err_l1homotopy,err_sparsa,err_kalman]));

else
    fprintf('Signal MSE: l1homotopy-%3.4g, sparsa-%3.4g.\n',err_l1homotopy,err_sparsa);
    fprintf('Signal SER (in dB): l1homotopy-%3.4g, sparsa-%3.4g.\n',-10*log10([err_l1homotopy,err_sparsa]));
end

%% plot signal and reconstruction error
x_len = min(length(x_vec),length(xh_vec))-(P-1)*N;
sig_len = min(length(sig_vec),length(sigh_vec))-(P-1)*N-L-N;
x_vec1 = x_vec(1:x_len);
xh_vec1 = xh_vec(1:x_len,1);
sig_vec1 = sig_vec(1:sig_len);
sigh_vec1 = sigh_vec(1:sig_len,1);

fig(123);
subplot(221);
imagesc(reshape(sig_vec1,N,sig_len/N));
axis tight;
title('original signal')
subplot(2,2,2)
imagesc(reshape(sigh_vec1,N,sig_len/N));
% plot((1:sig_len)/N,sigh_vec1-sig_vec1, 'LineWidth',1);
axis tight
title('reconstruction error')
subplot(2,2,3);
imagesc(reshape(x_vec1,N,x_len/N)); axis xy;
axis tight;
title('DWT coefficients');
colorbar
subplot(2,2,4);
imagesc(reshape(x_vec1-xh_vec1,N,x_len/N),[0 max(abs(x_vec1))/20]); axis xy
axis tight
title('reconstruction error');
colorbar


fig(3);
% view DWT coefficients...
alpha_vec1 = apply_DWT(sig_vec1,N,wType,J,sym);
alphah_vec1 = apply_DWT(sigh_vec1,N,wType,J,sym);
subplot(221); imagesc(reshape(alpha_vec1,N,length(alpha_vec1)/N));
axis xy;
title('original')
subplot(222); imagesc(reshape(alphah_vec1,N,length(alpha_vec1)/N));
axis xy;
title('reconstructed');
subplot(212); plot([alpha_vec1 alphah_vec1]);
title('comparison');

%%
if LM == N
    st_ind = 2*N;
    fig(4); clf; hold on;
    s_len = length(sig_vec(st_ind+1:N*sim));
    rshp1 = @(x) reshape(x(st_ind+1:N*sim),N,s_len/N);
    err_fun = @(x) -10*log10(sum(rshp1(sig_vec-x).^2,1)./sum(rshp1(sig_vec).^2,1));
    plot(err_fun(sigh_vec(:,1)),'b');
    % plot(err_fun(sigh_vec(:,2)),'k');
    plot(err_fun(sigh_vec(:,3)),'r');
    title('error evolution');
    xlabel('iteration');
    ylabel('ser in db');
    legend('L1','LS');
    shg
end