%% Shrinkage parameters selection
% Selects target weights for the next homotopy step such that the weights
% shrink at a faster rate on the active set than the inactive set. 

switch shrinkage_mode
    case 'Tsteps'       
        % active weights are reduced to tau/ewt
        tau_vec = ones(N,1)*tau;
        tau_vec(gamma_xh) = tau/ewt;
    case 'frac'
        % acitve weights are reduced to a fraction of their present value
        tau_vec = ones(N,1)*min(epsilon/ewt,max(epsilon_vec));
        tau_vec(gamma_xh) = epsilon_vec(gamma_xh)/ewt;
    case 'Trwt'
        % active weights are modified based on current signal estimate
        [alpha beta epsilon] = weight_param(5,1,xk_1,M);    
        ewt_b = beta; % M*(norm(xk_1,2)/norm(xk_1,1))^2;
        tau_vec = ones(N,1)*tau;
        % tau_vec(gamma_xh) = min([tau_vec(gamma_xh) tau./abs(xk_1(gamma_xh))/ewt_b],[],2);  
        tau_vec(gamma_xh) = min([tau_vec(gamma_xh) tau./(beta*abs(xk_1(gamma_xh))+epsilon)],[],2);
        % alpha = 0.75;
        % tau_vec = alpha*tau_vec+(1-alpha)*epsilon_vec;
    case 'rwt'
        % hybrid scheme...
        %  first stage: reduce active weights by a fraction to same values
        %  second stage: reduce active weights according to signal estimate
        %
        
        % this two-step procedure is often slower, but more stable, than Trwt
        if epsilon <= tau*ewt;
            rwt_step2 = 1;
        end
        if rwt_step2
            [alpha beta epsilon] = weight_param(5,1,xk_1,M);    
            ewt_b = beta; % M*(norm(xk_1,2)/norm(xk_1,1))^2;
            
            tau_vec = ones(N,1)*tau;
            % tau_vec(gamma_xh) = min([tau_vec(gamma_xh) norm(xk_1,2)/sqrt(N*M)*tau./abs(xk_1(gamma_xh)*ewt)],[],2);
            % how about the following selection for which w_i = \tau/beta, where beta = sqrt(M)/\|x\|_1*\|x\|_2 which implies that sqrt(M)> beta > sqrt(M/K)
            % tau_vec(gamma_xh) = min([tau_vec(gamma_xh) norm(xk_1,1)/M*tau./abs(xk_1(gamma_xh)*ewt)],[],2);
            % tau_vec(gamma_xh) = min([tau_vec(gamma_xh) (norm(xk_1,1)/norm(xk_1,2)/sqrt(M))^1*tau./abs(xk_1(gamma_xh)*ewt)],[],2);
            tau_vec(gamma_xh) = min([tau_vec(gamma_xh) tau./abs(xk_1(gamma_xh))/ewt_b],[],2);
        else            
            tau_vec = ones(N,1)*min(epsilon/ewt,max(epsilon_vec));
            % tau_vec(gamma_xh) = ones(length(gamma_xh),1)*epsilon_old;
            tau_vec(gamma_xh) = min([epsilon_vec(gamma_xh) epsilon/ewt*ones(length(gamma_xh),1)],[],2);
            % tau_vec(gamma_xh) = epsilon/ewt*ones(length(gamma_xh),1);
        end
    case 'OLS'
        % set according to LS solution on the active support
        %   such as w_gamma = 1./abs((A_gamma'*A_gamma)^-1*A_gamma'*y);
        % xols = (A_gamma'*A_gamma)^-1*A_gamma'*y;
        xols = zeros(N,1);
        %         switch delx_mode
        %             case 'mil'
        %                 xols(gamma_xh) = iAtA*(A(:,gamma_xh)'*y);
        %             case 'qr'
        %                 xols(gamma_xh) = R\(R'\(A(:,gamma_xh)'*y));
        %         end
        xols(gamma_xh) = A(:,gamma_xh)\y;
        [alpha beta epsilon] = weight_param(5,1,xols,M);   
        % alpha =1; epsilon = 1; beta = M*(norm(xk_1,2)/norm(xk_1,1))^2;
        tau_vec = ones(N,1)*tau;
        % tau_vec(gamma_xh) = min([tau_vec(gamma_xh) tau./(beta*abs(xols(gamma_xh))+epsilon)],[],2); 
        tau_vec(gamma_xh) = max([tau_vec(gamma_xh)./(beta*abs(xols(gamma_xh))+epsilon) epsilon_vec(gamma_xh)./(beta*abs(xols(gamma_xh))+epsilon)],[],2);
        figure(107); plot([x_orig xols xk_1]); pause(1/60);
end