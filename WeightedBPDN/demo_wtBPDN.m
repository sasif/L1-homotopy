% demo_wtBPDN
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% and dynamically update the weights w_i
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: March 16, 2011

% Comparison between IRW-H and ARW-H, 
% INIT: solution of standard BPDN problem (un-weighted L1)
% also included OMP, oracle-support LS, and oracle-weighted L1 

clear
close all force

%% parfor setup
% numCores = str2double(getenv('NUMBER_OF_PROCESSORS'));
% mpSize = numCores-1;
% if matlabpool('size') ~= mpSize
%     if matlabpool('size')~=0
%         matlabpool close;
%     else
%         matlabpool('open', mpSize);
%     end
% end

%% Setup path
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath ../Pursuits_Homotopy/
addpath ../utils/ 
addpath ../utils/utils_Wavelet/
addpath src/

disp(['--------------------',datestr(now),'-------------------------'])

% load RandomStates
%
% rseed = 0;
% rseed = sum(100*clock);
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));

rseed = 2012;
rand('state',rseed);
randn('state',rseed);

% simulation parameters
mType = 'randn'; % {'randn','orth','rdct'};
sType = 'blocks'; % {'randn','sign','highD', 'blocks','pcwPoly'}
SNR = 40;
rwt_mode = 5;   % mode for selecting parameters in iteratirve reweighting
lambda = 0;     % mode for selecting regularization parameter for BPDN

N = 512;   % signal length
R = 4;
M = round(N/R);    % no. of measurements
T = round(M/5);    % sparsity level

% reweighted setup
rwt = 5;        % number of reweighting iterations
rwt_adp = 0;    % number of reweighting iterations after adaptive reweighting

% rank-1 update mode 
delx_mode = 'mil'; % mil or qr

str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d, rwt_mode-%d, lambda%3.4g.', mType, sType, SNR, N, M, T, rwt_mode, lambda);
disp(str0);


%% setup waitbar
% h = waitbar(0,'1','Name','Adaptive wtBPDN...',...
%     'CreateCancelBtn',...
%     'setappdata(gcbf,''canceling'',1)');
% setappdata(h,'canceling',0)
% reverseStr = '';

maxsim = 10;
SIM_stack = cell(maxsim,1);
SIM_memory = cell(maxsim,1);
                                

%% Setting up figures for proper display... 
fig(1); 
Position = get(gcf,'Position');
ScreenSize = get(0,'ScreenSize');
Position(1) = 1;
Position(2) = ScreenSize(4)/3; 
Position(3:4) = ScreenSize(3:4)/2;
set(gcf,'Position',Position);
fig(2); 
Position = get(gcf,'Position');
ScreenSize = get(0,'ScreenSize');
Position(1) = ScreenSize(3)/2;
Position(2) = ScreenSize(4)/3; 
Position(3:4) = ScreenSize(3:4)/2;
set(gcf,'Position',Position);
for sim = 1:maxsim
    %     % Check for Cancel button press
    %     if getappdata(h,'canceling')
    %         break
    %     end
    %     % Report current estimate in the waitbar's message field
    %     waitbar(sim/maxsim,h,sprintf('sim %d of %d',sim,maxsim))
    
    % Generate a random signal
    in = []; in.type = sType; in.T = T; in.randgen = 1;
    switch sType
        case 'blocks'
            in.take_fwt = 1;
            in.wType = 'haar';
        case 'heavisine'
            in.take_fwt = 1;
            in.wType = 'daub4';
        case 'pcwpoly';
            in.take_fwt = 1;
            in.wType = 'daub8';
    end
            
    x = genSignal(N,in); 
    [val ind] = sort(abs(x),'descend');
    ind_pos = ind(val>1e-1);
    gamma_orig = ind_pos(1:min(length(ind_pos),M-1));
    

    % measurement matrix
    in = []; in.type = mType;
    A = genAmat(M,N,in);
    
    % measurements
    sigma = sqrt(norm(A*x)^2/10^(SNR/10)/M);
    %sigma = .05;
    e = randn(M,1)*sigma;
    y = A*x+e;
    
    %     % orthogonalize rows of A
    %     [Q, R] = qr(A',0);
    %     A = Q'; y = R' \ y;
    
    % parameter selection
    if lambda > 0
        tau = lambda*max(abs(A'*y));
    else
        % tau = sigma*sqrt(log(N));
        tau = max(1e-4*max(abs(A'*y)),sigma*sqrt(log(N)));
    end
    
    maxiter = 4*N;
    
    %% Adaptive support & weight selection (ARW-H in the paper)
    in = [];
    in.A = A; in.y = y; in.x = x;
    in.x_init = zeros(N,1); in.max_rwt = rwt_adp;
    in.tau = tau;
    in.rwt_mode = rwt_mode;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    if in.verbose
        fprintf('sType-%s, mType-%s, (N,M,T)-%d,%d,%d ',sType,mType,N,M,T);    
    end
    out = script_rwtBPDN_adaptive(in);
    xh_adp = out.x_out;
    iter_adp = out.iter;
    time_adp = out.time;
    gamma_adp = out.gamma;
    W_adp = out.W_new;
    err_adp = out.err;
    supp_diff_adp = out.supp_diff;
    
    
    %% Iterative reweighting (IRW-H in the paper)
    in = [];
    in.A = A; in.y = y; in.x = x;
    in.x_init = zeros(N,1); in.max_rwt = rwt;
    in.tau = tau;
    in.rwt_mode = rwt_mode;
    in.delx_mode = delx_mode;
    out = script_rwtBPDN_iterative(in);
    xh_rwt = out.x_out;
    xh_rwt_init = out.x_init;
    iter_rwt = out.iter;
    time_rwt = out.time;
    gamma_rwt = out.gamma;
    err_rwt = out.err;
    supp_diff_rwt = out.supp_diff;
    W_rwt = out.W_new;
    
    %% Oracle rwt BPDN
    alpha = 5; beta = 10;
    W_new = ones(N,1); W_new(gamma_orig) = 1/alpha./(beta*abs(x(gamma_orig)));    
    % epsilon = 0.1; W_new = 1/alpha./(beta*(abs(x))+epsilon);
    
    % To check the accuracy of rwtBPDN_adaptive set
    % W_new = W_adp/tau;
    
    in = [];
    in.tau = tau;
    in.maxiter = maxiter;
    in.x_orig = x;
    in.record = 1;
    in.delx_mode = delx_mode;
    %     in.Te = rwt_itr;
    AW = A*diag(1./W_new);
    out_new = BPDN_homotopy_function(AW, y, in); %BPDN
    xh_orac = out_new.x_out.*(1./W_new);
    gamma_orac = out_new.gamma;
    iter_orac = out_new.iter;
    time_orac = out_new.time;
    
    %% oracle LS
    x_LS = zeros(N,1);
    x_LS(gamma_orig) = A(:,gamma_orig)\y;
                                        
    %% OMP
    in = [];
    in.Te = round(1.2*T);
    out = OMP_function(y,A,in);
    x_omp = out.x_out;
    iter_omp = out.iter;
    gamma_omp = out.gamma;
    
    
    % fprintf('sim %d -- iter: %d; %d; %d; %d. Error: %3.4g; %3.4g; %3.4g; %3.4g; %3.4g. Diff-supp = %d; %d; %d; %d.\n', sim, iter_orac, iter_adp, iter_rwt, iter_omp, norm(x-x_LS)/norm(x), norm(x-xh_orac)/norm(x), norm(x-xh_adp)/norm(x), norm(x-xh_rwt)/norm(x), norm(x-x_omp)/norm(x), length(setdiff(gamma_orac,gamma_orig)), length(setdiff(gamma_adp,gamma_orig)),length(setdiff(gamma_rwt,gamma_orig)),length(setdiff(gamma_omp,gamma_orig)));
    % msg = sprintf('sim %d -- iter: %d; %d; %d; %d. Error: %3.4g; %3.4g; %3.4g; %3.4g; %3.4g. Diff-supp = %d; %d; %d; %d.\n', sim, iter_orac, iter_adp, iter_rwt, iter_omp, norm(x-x_LS)/norm(x), norm(x-xh_orac)/norm(x), norm(x-xh_adp)/norm(x), norm(x-xh_rwt)/norm(x), norm(x-x_omp)/norm(x), length(setdiff(gamma_orac,gamma_orig)), length(setdiff(gamma_adp,gamma_orig)),length(setdiff(gamma_rwt,gamma_orig)),length(setdiff(gamma_omp,gamma_orig)));
    % fprintf([reverseStr, msg]);
    % reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    % SIM_stack(sim,:) = [sim, iter_orac, iter_adp(1), sum(iter_adp,2), iter_rwt(1), sum(iter_rwt,2), iter_omp, norm(x-x_LS)/norm(x),  norm(x-xh_orac)/norm(x), err_adp(1), norm(x-xh_adp)/norm(x), err_rwt(1), norm(x-xh_rwt)/norm(x), norm(x-x_omp)/norm(x), length(setdiff(gamma_orac,gamma_orig)), length(setdiff(gamma_adp,gamma_orig)),length(setdiff(gamma_rwt,gamma_orig)),length(setdiff(gamma_omp,gamma_orig))];    
    % fprintf('sim %d -- iter: %d; %d; %d; %d. Error: %3.4g; %3.4g; %3.4g; %3.4g; %3.4g. Diff-supp = %d; %d; %d; %d.\n', sim, iter_orac, iter_adp, iter_rwt, iter_omp, norm(x-x_LS)/norm(x), norm(x-xh_orac)/norm(x), norm(x-xh_adp)/norm(x), norm(x-xh_rwt)/norm(x), norm(x-x_omp)/norm(x), length(setdiff(gamma_orac,gamma_orig)), length(setdiff(gamma_adp,gamma_orig)),length(setdiff(gamma_rwt,gamma_orig)),length(setdiff(gamma_omp,gamma_orig)));
    
    SIM_stack{sim} = [sim, tau, norm(x-x_LS)^2/norm(x)^2, norm(x-xh_orac)^2/norm(x)^2, iter_orac, ...
        norm(x-xh_adp)^2/norm(x)^2, sum(iter_adp,2), sum(time_adp,2), ...
        norm(x-xh_rwt)^2/norm(x)^2, sum(iter_rwt,2), sum(time_rwt,2), ...
        norm(x-x_omp)^2/norm(x)^2, iter_omp];  
    
    % adaptive rwt
    exp = 1;
    SIM_memory{sim}{exp,1} = 'ARW-H';
    SIM_memory{sim}{exp,2} = iter_adp;
    SIM_memory{sim}{exp,3} = err_adp;
    SIM_memory{sim}{exp,4} = supp_diff_adp;
    SIM_memory{sim}{exp,5} = time_adp;
    % iterative rwt
    exp = exp+1;
    SIM_memory{sim}{exp,1} = 'IRW-H';
    SIM_memory{sim}{exp,2} = iter_rwt;
    SIM_memory{sim}{exp,3} = err_rwt;
    SIM_memory{sim}{exp,4} = supp_diff_rwt;
    SIM_memory{sim}{exp,5} = time_rwt;

    mS =  SIM_stack{sim};
    fprintf('sim %d. tau = %3.4g, oracLS=%3.4g. (err,iter,time): oracWT-%3.4g,%3.4g; ARW-H-%3.4g,%3.4g,%3.4g; IRW-H-%3.4g,%3.4g,%3.4g; OMP-%3.4g,%3.4g. \n', sim, mS(2:end));

    %% plot recovered signals
    fig(1); clf
    subplot(221); plot([x xh_adp xh_rwt xh_rwt_init]); title('ORIG .. ARW-H .. IRW-H .. INIT estimates');
    subplot(222); plot([W_adp W_rwt]); title('Weights: ARW-H .. IRW-H'); 
    subplot(223); plot(xh_adp,x,'bd'); title('Original vs ARW-H estimate'); axis square
    subplot(224); plot(xh_rwt_init, x,'k.'); 
    hold on; plot(xh_rwt, x,'rs');    title('Original vs INIT and IRW-H estimate'); axis square;
    fig(2); clf
    hold on; plot(xh_rwt_init, x,'k.');    
    hold on; plot(xh_rwt, x,'rs');        
    plot(xh_adp,x,'bd');
    XLim = get(gca,'XLim');
    XLim(1) = sign(XLim(1))*max(abs(XLim(:)));
    XLim(2) = sign(XLim(2))*max(abs(XLim(:)));
    set(gca,'XLim',XLim);
    set(gca,'YLim',XLim);    
    axis square; shg    
    title('Comparison betweeen INIT .. IRW-H .. ARW-H and the original signal')
    legend('INIT','IRW-H','ARW-H','Location','NorthWest');
    
    if sim == 1
        fprintf('\n Hit any key to continue... resize figures (if you want) \n')    
        pause;   
    end
end

mS =  mean(cell2mat(SIM_stack),1);
fprintf('Average results: maxsim %d. tau = %3.4g, oracLS=%3.4g. (err,iter,time): oracWT-%3.4g,%3.4g; ARW-H-%3.4g,%3.4g,%3.4g; IRW-H-%3.4g,%3.4g,%3.4g; OMP-%3.4g,%3.4g. \n', maxsim, mS(2:end));

%% Plot results...
EXP = {'ARW-H','IRW-H'};
exp_plot = [1 2];
fig(101);
for exp_no = 1:length(exp_plot)
    exp = exp_plot(exp_no);
    sim_err = [];
    sim_itr = [];
    sim_time = [];
    sim_supp_diff = [];
    for ii = 1:maxsim
        sim_err = [sim_err; [SIM_memory{ii}{exp,3}]];
        sim_itr = [sim_itr; [SIM_memory{ii}{exp,2}]];
        sim_time = [sim_time; [SIM_memory{ii}{exp,5}]];
        sim_supp_diff = [sim_supp_diff; [SIM_memory{ii}{exp,4}]];
    end
    subplot(2,length(exp_plot),exp_no); plot(sim_err); title(sprintf('rel. error for %s',char(EXP{exp})));
    subplot(2,length(exp_plot),length(exp_plot)+exp_no); plot(sim_itr); title(sprintf('iter. count for %s',char(EXP{exp})));
    % subplot(3,length(exp_plot),2*length(exp_plot)+exp_no); plot(sim_supp_diff); title(sprintf('supp diff for %s',char(EXP{exp})));    
end

% delete(h)       
% DELETE the waitbar; don't try to CLOSE it.
% close all force % After Ctrl+C