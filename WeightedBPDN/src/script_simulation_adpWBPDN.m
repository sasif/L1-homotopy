% partial script for running different solvers for adaptive reweighting

% reweighted setup
adp_wt = 0;  % adaptive selection of the weights?
rwt = 3*(1-adp_wt)+1;        % number of reweighting iterations

I_ALL = [];
EXP = [0 1; 0 4; 1 1]; 
% 0 1 -- standard BPDN
% 0 x -- iterative reweighting with x reweighting iterations
% 1 1 -- adaptive reweigthing 
for exp = 1:3;
    % set up random generator
    rseed = 2012;
    rand('state',rseed);
    randn('state',rseed);

    % RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));
    % seed = rng('shuffle');
    % seed = rng(0);

    adp_wt = EXP(exp,1);
    rwt = EXP(exp,2);
    
    fprintf('adp_wt = %d, rwt = %d ',adp_wt,rwt);
    
    
    err_fun = @(z) -20*log10(z);
    
    %% Main loop
    SIM_stack = cell(maxsim,1);
    SIM_memory = cell(maxsim,1);
    
    Ir_stack = [];
    for sim = 1:maxsim
        %% Generate a random signal
        in = []; in.type = sType; in.T = T; in.randgen = 1; in.take_fwt = 1;
        [x img wave_struct] = genSignal(N,in);
        % gamma_orig = find(abs(x)>0);
        [val ind] = sort(abs(x),'descend');
        ind_pos = ind(val>0);
        gamma_orig = ind_pos(1:min(length(ind_pos),M-1));
        
        % measurement matrix
        in = []; in.type = mType;
        AA = genAmat(M,N,in);
        if largescale
            wType = wave_struct.wType;
            FilterBank = wave_struct.FilterBank;
            h0 = FilterBank(1,:);
            h1 = FilterBank(2,:);
            g0 = FilterBank(3,:);
            g1 = FilterBank(4,:);
            J = wave_struct.J;
            sym = wave_struct.sym;
            
            n = sqrt(N);
            rshp1 = @(z) reshape(z,n,n);
            rshp2 = @(z) z(:);
            
            % FT = @(x) rshp2(fwt2(rshp1(x),h0,h1,J,sym));
            % iFT = @(x) rshp2(ifwt2(rshp1(x),g0,g1,J,sym));
            FT = @(z) z; 
            iFT = @(z) z; 
            
            A = @(x) AA(iFT(x),1);
            AT = @(x) FT(AA(x,2));
            % measurements
            Ax = A(x);
            sigma = sqrt(norm(Ax)^2/10^(SNR/10)/M);
            %sigma = .05;
            e = randn(M,1)*sigma;
            
            y = Ax+e;
            Aty = AT(y);
        else
            A = AA;
            AT = A';
            % measurements
            Ax = A*x;
            sigma = sqrt(norm(Ax)^2/10^(SNR/10)/M);
            %sigma = .05;
            e = randn(M,1)*sigma;
            y = A*x+e;
            Aty = AT*y;
        end
        
        maxiter = 2*N;
        
        %% wSPGL1
        % call SPGL1
        iter_spgl1 = []; err_spgl1 = []; supp_diff_spgl1 = []; time_spgl1 = [];
        % if ~exist('spgSetParms','file'); error('Solver SPGL1 is not found.'); end
        
        W_new = ones(N,1);
        delta = norm(e);
        delta = max(sqrt(M)*sigma,1e-1);
        x_spgl1 = [];
        
        for rwt_itr = 1:rwt*0
            if adp_wt == 0
                spg_opts = spgSetParms('verbosity',0,'weights',W_new,'iterations',maxiter);
                [x_spgl1,r,g,info] = spgl1(AA,y,0,delta,x_spgl1,spg_opts);
                % [x_spgl1,r,g,info] = spgl1(AA,y,1e5,0,x_spgl1,spg_opts);
            else
                spg_opts = wspgSetParms('verbosity',0,'omega',0.3);
                [x_spgl1,r,g,info] = wspgl1(AA,y,0,delta,x_spgl1,spg_opts);
            end
            
            rerr = norm(x_spgl1-x)/norm(x);
            iter_spgl1 = [iter_spgl1  (info.nProdA+info.nProdAt)/2];
            err_spgl1 = [err_spgl1 rerr];
            supp_diff_sgpl1 = [supp_diff_spgl1 length(setxor(gamma_orig,find(x_spgl1)))];
            time_spgl1 = [time_spgl1 info.timeTotal];
            [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_spgl1,M);
            W_new = 1/alpha./(beta*abs(x_spgl1)+epsilon);
            % W_new = 1./(abs(x_spgl1)+0.1);
            r_spgl1 = A(x_spgl1)-y; Atr_spgl1 = AT(r_spgl1);        
        end
        
        %% SpaRSA
        tau = max(1e-3*max(abs(Aty)),sigma*sqrt(log(N)));
        % tau = max(abs(Atr_spgl1));
        
        iter_sparsa = [];
        err_sparsa = [];
        time_sparsa = [];
        x_sparsa = 0;
        supp_diff_sparsa = [];
        W_new = 1;
        
        % [alpha beta epsilon] = weight_param(5,rwt_itr,x,M);
        % W_new = 1/alpha./(beta*abs(x)+epsilon);
        % W_new = abs(Atr_spgl1)/tau;
        
        for rwt_itr = 1:rwt
            psi_function = @(x,tau) soft(x,tau*W_new);
            phi_function = @(x) sum(abs(W_new.*x));
            tic
            [x_SpaRSA,x_debias_SpaRSA_m,obj_SpaRSA_m_cont,...
                times_SpaRSA_m_cont,debias_start_SpaRSA_m,mse_SpaRSA_m,taus_SpaRSA_m, numA, numAt]= ...
                SpaRSA_adpW(y,A,tau,...
                'Monotone',0,...
                'AT',AT,...
                'adp_wt',adp_wt,...
                'Debias',0,...
                'Initialization',x_sparsa,...
                'StopCriterion',2,...
                'ToleranceA',1e-4,...
                'MaxiterA',maxiter,...
                'psi',psi_function,...
                'phi',phi_function,...
                'Verbose',0,...
                'W_new',W_new,...
                'True_x',x,...
                'Continuation',1,...
                'Continuationsteps',-1);
            x_sparsa = x_SpaRSA;
            time_sparsa = [time_sparsa toc];
            
            iter_sparsa = [iter_sparsa (numA+numAt)/2];
            rerr = norm(x-x_sparsa)/norm(x);
            err_sparsa = [err_sparsa rerr];
            
            supp_diff_sparsa = [supp_diff_sparsa length(setxor(gamma_orig,find(x_sparsa)))];
            
            [alpha beta epsilon] = weight_param(rwt_mode,rwt_itr,x_sparsa,M);
            W_new = 1/alpha./(beta*abs(x_sparsa)+epsilon);
            
            %         gamma = find(abs(x)>tau*100);
            %         W_new = ones(N,1);
            %         ewt_b = 2*length(y)*(norm(x,2)/norm(x,1))^2;
            %         W_new(gamma) = 1./abs(x(gamma))/ewt_b;
        end
        % r_sparsa = A(x_sparsa)-y; Atr_sparsa = AT(r_sparsa);
        
        %% NESTA
        %     opts = [];
        %     opts.Verbose = 0;
        %     opts.tolvar = 1e-5;
        %     opts.MaxIntIter = 5;
        %     delta = max(sigma*sqrt(M),1e-6);
        %     muf = 1e-3;
        %
        %     iter_nesta = [];
        %     err_nesta = [];
        %     supp_diff_nesta = [];
        %     x_nesta = [];
        %     tic;
        %     % constrained
        %     [x_nesta,niter,resid,outData,optsOut] = NESTA_adpW(A,AT,y,muf,delta,opts);
        %     % unconstrained
        %     % La = norm( A*A' );
        %     % [x_nesta,niter,resid,outData,optsOut] = NESTA_UP(A,[],y,tau,La,muf,opts);
        %     time_nesta = toc;
        %     iter_nesta = [iter_nesta niter*3];
        %     rerr = norm(x-x_nesta)/norm(x);
        %     err_nesta = [err_nesta rerr];
        %     supp_diff_nesta = [supp_diff_nesta length(setxor(gamma_orig,find(x_nesta)))];
        
        if largescale == 0
            %% Adaptive support & weight selection
            in = [];
            in.A = A; in.y = y; in.x = x;
            in.x_init = zeros(N,1); in.max_rwt = 0;
            in.tau = tau;
            in.rwt_mode = rwt_mode;
            in.delx_mode = delx_mode;
            in.debias = 0;
            in.verbose = 0;
            out = script_rwtBPDN_adaptive(in);
            xh_adp = out.x_out;
            iter_adp = out.iter;
            gamma_adp = out.gamma;
            W_new = out.W_new;
            xh_adp_init = out.x_init;
            err_adp = out.err;
            time_adp = out.time;
            supp_diff_adp = out.supp_diff;
            
            %% Oracle rwt BPDN
            alpha = 5; beta = 10;
            W_new = ones(N,1); W_new(gamma_orig) = 1/alpha./(beta*abs(x(gamma_orig)));
            % epsilon = 0.1; W_new = 1/alpha./(beta*(abs(x))+epsilon);
            
            % To check the accuracy of rwtBPDN_adaptive set
            % W_new = out.W_new/tau;
            
            in = [];
            in.tau = tau;
            in.maxiter = maxiter;
            in.x_orig = x;
            in.record = 1;
            in.delx_mode = 'mil';
            %     in.Te = rwt_itr;
            AW = A*diag(1./W_new);
            tic;
            out_new = BPDN_homotopy_function(AW, y, in); %BPDN
            % time_orac = out_new.time;
            time_orac = toc;
            xh_orac = out_new.x_out.*(1./W_new);
            gamma_orac = out_new.gamma;
            iter_orac = out_new.iter;
            
            
            %% oracle LS
            x_LS = zeros(N,1);
            x_LS(gamma_orig) = A(:,gamma_orig)\y;
        end
                
        %%
        if largescale
            fprintf('Image-%s-%d, wavelet-%s, adp_wt = %d @ rwt = %d ... ',sType,n,wType,adp_wt, rwt);
            
            I = ifwt2(reshape(x,n,n),g0,g1,J,sym); Ir_max = max(abs(I(:)));
            Ir_sparsa = ifwt2(reshape(x_sparsa,n,n),g0,g1,J,sym);
            % Ir_sparsa(Ir_sparsa<0) = 0;
            % Ir_sparsa(Ir_sparsa>255) = 255;
            
            mse_sparsa = norm(Ir_sparsa(:)-I(:))^2/N;
            psnr_sparsa = 10*log10(Ir_max*Ir_max/mse_sparsa);
            ser_sparsa = norm(I(:))^2/norm(Ir_sparsa(:)-I(:))^2;
            fprintf('(PSNR,MSE,SER): SpaRSA-%3.4g,%3.4g,%3.4g  ', psnr_sparsa, mse_sparsa, ser_sparsa);
            
            % Ir_spgl1 = ifwt2(reshape(x_spgl1,n,n),g0,g1,J,sym);
            % mse_spgl1 = norm(Ir_spgl1(:)-I(:))^2/N;
            % psnr_spgl1 = 10*log10(Ir_max*Ir_max/mse_spgl1);
            % ser_spgl1 = norm(I(:))^2/norm(Ir_spgl1(:)-I(:))^2;
            % fprintf('  SPGL1-%3.4g,%3.4g,%3.4g.',psnr_spgl1, mse_spgl1, ser_spgl1)
            
            fprintf('\n');
            
            Ir_stack = cat(3,Ir_stack,Ir_sparsa);
        end
        
        %% Save results...
        
        % fprintf('sim iter. %d: \n',sim);
        if largescale            
            SIM_stack{sim} = [sim, tau, ...
                psnr_sparsa, mse_sparsa, ser_sparsa, sum(iter_sparsa,2), sum(time_sparsa,2)];
            %   psnr_spgl1, mse_spgl1, ser_spgl1, sum(iter_spgl1,2), sum(time_spgl1,2)];
            % norm(x-x_nesta)/norm(x), sum(iter_nesta,2), sum(time_nesta,2)];
        else
            SIM_stack{sim} = [sim, tau, (norm(x-x_LS)/norm(x)), (norm(x-xh_orac)/norm(x)), iter_orac, ...
                (norm(x-xh_adp)/norm(x)), sum(iter_adp,2), sum(time_adp,2), ...
                (norm(x-x_sparsa)/norm(x)), sum(iter_sparsa,2), sum(time_sparsa,2)]; 
            %     (norm(x-x_spgl1)/norm(x)), sum(iter_spgl1,2), sum(time_spgl1,2)];
            % norm(x-x_nesta)/norm(x), sum(iter_nesta,2), sum(time_nesta,2)];
        end
        index = 0;
        if largescale == 0
            % adaptive rwt
            index = index+1;
            SIM_memory{sim}{index,1} = 'adp';
            SIM_memory{sim}{index,2} = iter_adp;
            SIM_memory{sim}{index,3} = err_adp;
            SIM_memory{sim}{index,4} = supp_diff_adp;
            SIM_memory{sim}{index,5} = time_adp;
        end
        
        % sparsa
        index = index+1;
        SIM_memory{sim}{index,1} = 'sparsa';
        SIM_memory{sim}{index,2} = iter_sparsa;
        SIM_memory{sim}{index,3} = err_sparsa;
        SIM_memory{sim}{index,4} = supp_diff_sparsa;
        SIM_memory{sim}{index,5} = time_sparsa;
        % spgl1
        index = index+1;
        SIM_memory{sim}{index,1} = 'spgl1';
        SIM_memory{sim}{index,2} = iter_spgl1;
        SIM_memory{sim}{index,3} = err_spgl1;
        SIM_memory{sim}{index,4} = supp_diff_spgl1;
        SIM_memory{sim}{index,5} = time_spgl1;
        %     % nesta
        %     index = index+1;
        %     SIM_memory{sim}{index,1} = 'nesta';
        %     SIM_memory{sim}{index,2} = iter_nesta;
        %     SIM_memory{sim}{index,3} = err_nesta;
        %     SIM_memory{sim}{index,4} = supp_diff_nesta;
        %     SIM_memory{sim}{index,5} = time_nesta;
        %     % salsa
        %     index = index+1;
        %     SIM_memory{sim}{6,1} = 'salsa';
        %     SIM_memory{sim}{6,2} = iter_salsa;
        %     SIM_memory{sim}{6,3} = err_salsa;
        %     SIM_memory{sim}{6,4} = supp_diff_salsa;
        
        % mS = mean(cell2mat(SIM_stack),1);
        % fprintf('maxsim %d. oracLS=%3.4g. (err,iter,time): oracWT-%3.4g,%3.4g; adp-%3.4g,%3.4g,%3.4g; rwt-%3.4g,%3.4g,%3.4g; yall1-%3.4g,%3.4g,%3.4g; sparsa-%3.4g,%3.4g,%3.4g; omp-%3.4g,%3.4g; fpc-%3.4g,%3.4g; spgl1-%3.4g,%3.4g. \n', maxsim, mS(2:end));
        
    end
    
    
    %%
    mS = mean(cell2mat(SIM_stack),1);
    if largescale == 0
        % str1 = 'maxsim %d. tau = %3.4g, oracLS=%3.4g. (err,iter,time): oracWT-%3.4g,%3.4g; adp-%3.4g,%3.4g,%3.4g; sparsa-%3.4g,%3.4g,%3.4g; spgl1-%3.4g,%3.4g,%3.4g; nesta-%3.4g,%3.4g,%3.4g;';
        str1 = 'maxsim %d. tau = %3.4g, oracLS=%3.4g. (err,iter,time): oracWT-%3.4g,%3.4g; adp-%3.4g,%3.4g,%3.4g; sparsa-%3.4g,%3.4g,%3.4g; spgl1-%3.4g,%3.4g,%3.4g';
    else
        % str1 = 'maxsim %d. tau = %3.4g, sparsa-%3.4g,%3.4g,%3.4g; spgl1-%3.4g,%3.4g,%3.4g; nesta-%3.4g,%3.4g,%3.4g;';
        str1 = 'maxsim %d. tau = %3.4g,\n~sparsa-%3.4g,%3.4g,%3.4g,%3.4g,%3.4g;'; %\n spgl1-%3.4g,%3.4g,%3.4g';
    end
    str2 = sprintf([str1,' \n'], maxsim, mS(2:end));
    fprintf(str2);       
    fprintf('(N,M,SNR)-%d,%d,%d, adp-%d,rwt-%d, PSNR: $%2.4g$dB --- ($%d$) \n',N,M,SNR,adp_wt,rwt,10*log10(Ir_max*Ir_max/mS(4)),round(mS(6)));
    
    Ir_avg = mean(Ir_sparsa,3);
    
    figure(100); imagesc(I,[0 255]); colormap gray; axis image off;        
    figure(101); imagesc(Ir_avg,[0 255]); colormap gray; axis image off;
    figure(102); imagesc(Ir_avg-I,[0 255/10]); colormap gray; axis image off;
    % figure(103); imagesc(Ir_spgl1,[0 255]); colormap gray; axis image off;
    % figure(104); imagesc(Ir_spgl1-I,[0 255/10]); colormap gray; axis image off;
    
    if SAVE_RESULTS
        figure(100); set(gcf, 'Color', 'w');
        eval(sprintf('export_fig %s.pdf',sType));
        file_name = sprintf('%s_N%d_M%d_SNR%d_adp%d_rwt%d_%s',sType,N,M,SNR,adp_wt,rwt,wType);
        figure(101); set(gcf, 'Color', 'w');
        eval(sprintf('export_fig %s.pdf',file_name));
        % file_name = sprintf('%s_N%d_M%d_SNR%d_adp%d_rwt%d_PSNR%3.4g_ITER%d_TIME%3.4g_spgl1',sType,N,M,SNR,adp_wt,rwt,psnr_spgl1,round(sum(iter_spgl1,2)),sum(time_spgl1,2));
        % figure(103); set(gcf, 'Color', 'w');
        % eval(sprintf('export_fig %s.pdf',file_name));
    end

    fprintf('\n');    
        
    I_ALL = [I_ALL Ir_avg];
    
    figure(200);
    % subaxis(2,3,exp,'SpacingVert',0.01,'MT',0.05,'MB',0.05, 'ML',0.01,'MR',0.01);
    subplot(2,3,exp)
    imagesc(Ir_sparsa, [0 255]); axis off image; colormap gray
    % subaxis(2,3,exp+3,'SpacingVert',0.01,'MT',0.05,'MB',0.05, 'ML',0.01,'MR',0.01);
    subplot(2,3,exp+3)    
    imagesc(abs(I-Ir_sparsa), [0 255/10]); axis off image; colormap gray
end

if SAVE_RESULTS
    figure(200); set(gcf, 'Color', 'w');
    file_name = sprintf('%s_N%d_M%d_SNR%d_%s_ALL',sType,N,M,SNR,wType);
    eval(sprintf('export_fig %s.pdf',file_name));    
end