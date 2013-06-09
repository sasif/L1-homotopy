% demo_streamingDWT
%
% Solves the following streaming problem
% min_a  \|a\|_1 + 1/2*||Phi Psi a + y||_2^2
%
% with Psi as an overlapping representation matrix and 
% by adding and removing one set of measurements at every time instant
%
% Applications:
%           streaming signal recovery using (overlapping) DWT
%           compare recovery using block-based DWT (sym=0) and
%           overlapping DWT (sym = 3) near line 80
% We can also add any other regularization operator in the reconstruction
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: March 2013
%
% TODO: this code needs some debugging (not stable)... 

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
sType = 'pcwreg'; % % menu('Choose function:', 'heavisine', 'pcwpoly (N>256,daub4)', 'pcwreg (N>128,daub4)', 'doppler', ...);
SNR = 35;       % additive Gaussian noise

% system parameters
N = 256;   % signal length
R = 4;      % compression rate
M = round(N/R);    % no. of measurements

LM = N; % LM: length of measurement window 

% streaming window
P = 5; % size of the working window is P*N

% signal length
sig_length = 2^15; % 128*128;
Np = N; % interval of kalman-type signal evolution

% signal extension for prediction of new coefficients in the streaming window.
eType = 'per';

% signal dynamics
dType = 'static'; % type of dynamics 'crshift', or 'static'
cshift = -1;
rshift_max = 0.5;
rshift_h = @(z) (rand-0.5)*rshift_max*2;

% DWT parameters
% type of scaling function
% depth of scaling functions (number of wavelet levels)
% type of extension for DWT (0 - periodic extension, 3 - streaming)
%
%-------------------------------------------------------------------------
% COMPARE recovery using block DWT (sym=0) and overlapping DWT (sym = 3)
%-------------------------------------------------------------------------
wType = 'daub4'; sym = 0;
% wType = 'daub4'; sym = 3;
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

fprintf('Streaming signal measureed and reconstructed with overlapping windows using DWT ... \n');
str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,R) = %d, %d, %d, P = %d, LM = %d, \n wType-%s, J = %d, sym = %d, specified signal-length = %d, \n dType-%s, cshift = %d, rshift_max = %3.4g, Np = %d.', mType, sType, SNR, N, round(N/R), R, P, LM, wType, J, sym, sig_length, dType, cshift, rshift_max, Np);
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
if strcmpi(dType,'crshift')    
    % Generate a random signal
    in = []; in.type = sType; in.randgen = 0; in.take_fwt = 0; 
    [x_init sig wave_struct] = genSignal(N,in);
    
    F_h = @(x,cshift,rshift) interp1(1:N,circshift(x,cshift),[1:N]+rshift,'linear','extrap')';
    
    F0 = zeros(N);
    for ii = 1:Np;
        F0(:,ii) = F_h(circshift([1; zeros(N-1,1)],ii-1),cshift,0);
    end
    maxsim = round(sig_length)/N;
    sigt = sig; sig = [];
    for ii = 1:maxsim;
        rshift = rshift_h(1);
        sigt = F_h(sigt, cshift, rshift);
        sig = [sig; sigt];
    end
else
    % Generate a random signal
    in = []; in.type = sType; in.randgen = 0; in.take_fwt = 0; in.Np = Np;
    [x_init sig wave_struct] = genSignal(sig_length,in);
end
sig = [zeros(L-N,1); sig]; 
sig_length = length(sig);
 
% view DWT coefficients...
alpha_vec = apply_DWT(sig,N,wType,J,sym);
figure(123); 
subplot(211); imagesc(reshape(alpha_vec,N,length(alpha_vec)/N)); 
axis xy;
subplot(212); plot(alpha_vec); 
 

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
% in.P = P; in.Jp = Jp; in.wType = wType; in.Np = Np; in.sym = sym;
PSI = create_PSI_DWT(in);

%% Setup sensing matrices
in = []; in.type = mType; in.parts = 4;
genAmat_h = @(M,N) genAmat(M,LM,in);
in.P = P-(LM-N)/N;
in.LM = LM; in.M = M; in.N = N;
PHI = create_PHI(in);

%% Overall streaming system matrix
T_length = size(PSI,1);
t_ind = 1:T_length;

sigt = sig(t_ind); % Signal under the LOT window at time t
if sym == 1 || sym == 2
    x = pinv(PSI'*PSI)*(PSI'*sigt); % Sparse LOT coefficients
else
    x = PSI'*sigt;
end

PSI_M = PSI(1:end-(L-N),:);
A = PHI*PSI_M;
y = PHI*sigt(1:end-(L-N));

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

s_ind = 1:L;
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
    y = PHI*sigt(1:end-(L-N));
    
    % Sparse LOT coeffs. irrespective of SM
    if sym == 1 || sym == 2
        x = pinv(PSI'*PSI)*(PSI'*sigt);
    else
        x = PSI'*sigt;
    end
    
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
        
          
        % REMOVE the part of outgoing DWT projection in the overlapping region
        % on left side of streaming window... 
        if LeftProj_cancel
            y(1:M) = y(1:M)-PHI(1:M,1:(L-N))*(Psi(end-(L-N)+1:end,:)*xh_old(1:N));
        end
        
        %% Update the signal estimate (for warm start)
        xh(1:end-N) = xh(N+1:end);
        % truncate small values in xh
        % xh(abs(xh)<tau/sqrt(log(P*N))) = 0;
        
        % --------------------------------------------%
        %  Linear prediction for the incoming samples %
        % --------------------------------------------%
        %
        %         xh1 = xh(end-N+1:end);
        %
        %         if strcmpi(eType,'per')
        %             % periodic extension
        %             xh1 = xh(end-N+1:end);
        %         else
        %             sig_temp = PSI_M(1:end-N,1:(P-1)*N)*xh(1:end-N);
        %             sig_temp(1:L-N) = sig_temp(1:L-N)+Psi(end-(L-N)+1:end,:)*xh_old(1:N);
        %             % if sim == 1
        %             %     fprintf('Extension type: %s \n',eType);
        %             % end
        %             switch eType
        %                 case 'asym'
        %                     % anti-symmetric extension
        %                     sig_temp2 = [sig_temp; -sig_temp(end:-1:1); sig_temp];
        %                 case 'sym'
        %                     % symmetrci extension
        %                     sig_temp2 = [sig_temp; sig_temp(end:-1:1); sig_temp];
        %                 otherwise
        %                     disp('cant do it sire');
        %             end
        %             sig_temp2 = sig_temp2(1:size(PSI,1));
        %             if sym == 1 || sym == 2
        %                 xh_temp = pinv(PSI'*PSI)*(PSI'*sig_temp2);
        %             else
        %                 xh_temp = PSI'*sig_temp2;
        %             end
        %             xh1 = xh_temp(end-N+1:end);
        %         end
        %
        %         gamma1 = find(abs(xh1)>tau);
        %         [val ind] = sort(abs(xh1),'descend');
        %         gamma1 = ind(1:min(length(gamma1),ceil(M/log2(N))));
        %
        %         yt = y(end-M+1:end)-A(end-M+1:end,1:end-N)*xh(1:end-N);
        %         At = A(end-M+1:end,end-N+1:end);
        %         gamma3 = gamma1;
        %         %
        %         % NOTE TODO: Need to replace prediction with rank-update...
        %         %
        %         % gamma_all = [find(xh(1:end-N)); gamma3+(P-1)*N];
        %         % xh_all = 0*xh;
        %         % xh_all(gamma_all) = A(:,gamma_all)\y; % pinv(A(:,gamma_all)'*A(:,gamma_all))*(A(:,gamma_all)'*y-W(gamma_all).*sign(xh(gamma_all)));
        %
        %         % if SM == 3
        %         %     Aty = At'*yt;
        %         %     [val ind] = sort(abs(Aty),'descend');
        %         %     gamma2 = ind(1:min(length(gamma1),ceil(M/4)));
        %         %     gamma3 = union(gamma1,gamma2);
        %         % end
        %         %
        %         xh_t1 = 0*xh(end-N+1:end);
        %         xh_t1(gamma3) = pinv(At(:,gamma3)'*At(:,gamma3)+2*tau*eye(length(gamma3)))*(At(:,gamma3)'*yt);
        %         xh_t1(abs(xh_t1)<tau/sqrt(log(P*N))) = 0;
        %         xh(end-N+1:end) = xh_t1;
        
        
        % Truncate small coefficients... ?
        xh(abs(xh)<tau/sqrt(log(P*N))) = 0;
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
    
	%% Plot DWT coeffs. on the window
    fig(1); subplot(211);
    plot([x xh_streamingRWT x_sparsa x_yall1]);
    title('Comparison betweeen the original and reconstructed signal')
    
	%% Reconstructed signal    
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
	xh_vec((sim-1)*N+1:sim*N,3) = x_yall1(1:N);
    
    
    s_ind = t_ind(1:L);
    sig_vec(s_ind) = sigt(1:L);
    sigh_vec(s_ind,1) = sigh_vec(s_ind,1)+Psi*xh(1:N);
    sigh_vec(s_ind,2) = sigh_vec(s_ind,2)+Psi*x_sparsa(1:N);
    sigh_vec(s_ind,3) = sigh_vec(s_ind,3)+Psi*x_yall1(1:N);

    
    
    %% plot recovered signals
    fig(1); subplot(212)
    plot([sig_vec(1:s_ind(end)) sigh_vec(1:s_ind(end),1)]);
 
    drawnow;
end

mS =  sum(cell2mat(SIM_stack),1);
fprintf('Summed results: streaming_iter %d. tau = %3.4g, \n solver-(err,iter,time): \n streamingRWT homotopy-%3.4g,%3.4g,%3.4g; \n SpaRSA-%3.4g,%3.4g,%3.4g; \n YALL1-%3.4g,%3.4g,%3.4g. \n', streaming_iter, mS(2:end));
mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: streaming_iter %d. tau = %3.4g, \n solver-(err,iter,time): \n streamingRWT homotopy-%3.4g,%3.4g,%3.4g; \n SpaRSA-%3.4g,%3.4g,%3.4g; \n YALL1-%3.4g,%3.4g,%3.4g. \n', streaming_iter, mS(2:end));

% l1homotopy-%3.4g,%3.4g,%3.4g;
st_ind = 2*N;
err_l1homotopy = norm(sig_vec(st_ind+1:N*sim)-sigh_vec(st_ind+1:N*sim,1))^2/norm(sig_vec(st_ind+1:N*sim))^2;
err_sparsa = norm(sig_vec(st_ind+1:N*sim)-sigh_vec(st_ind+1:N*sim,2))^2/norm(sig_vec(st_ind+1:N*sim))^2;
err_yall1 = norm(sig_vec(st_ind+1:N*sim)-sigh_vec(st_ind+1:N*sim,3))^2/norm(sig_vec(st_ind+1:N*sim))^2;
% err_l1homotopy = norm(sig_vec(L-N+1:N*sim)-sigh_vec(L-N+1:N*sim,1))^2/norm(sig_vec(L-N+1:N*sim))^2;
% err_sparsa = norm(sig_vec(L-N+1:N*sim)-sigh_vec(L-N+1:N*sim,2))^2/norm(sig_vec(L-N+1:N*sim))^2;
% err_yall1 = norm(sig_vec(L-N+1:N*sim)-sigh_vec(L-N+1:N*sim,3))^2/norm(sig_vec(L-N+1:N*sim))^2;
% fprintf('Signal error: l1homotopy-%3.4g, sparsa-%3.4g, yall1-%3.4g.\n',err_l1homotopy,err_sparsa,err_yall1); 
fprintf('Signal MSE: l1homotopy-%3.4g, sparsa-%3.4g, yall1-%3.4g.\n',([err_l1homotopy,err_sparsa,err_yall1]));
fprintf('Signal SER (in dB): l1homotopy-%3.4g, sparsa-%3.4g, yall1-%3.4g.\n',-10*log10([err_l1homotopy,err_sparsa,err_yall1]));



%% plot signal and reconstruction error
x_len = min(length(x_vec),length(xh_vec))-(P-1)*N;
sig_len = min(length(sig_vec),length(sigh_vec))-(P-1)*N-L-N;
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
title('DWT coefficients');
colorbar 
subplot(2,2,4); 
imagesc(reshape(x_vec1-xh_vec1,N,x_len/N),[0 max(abs(x_vec1))/20]); axis xy
axis tight
title('reconstruction error (DWT)');
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