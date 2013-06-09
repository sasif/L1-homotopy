function [alpha beta epsilon] = weight_param(rwt_mode,itr,varargin)
% Compute alpha,beta, and epsilon to update weights as 
% w(i) = tau/alpha/(beta*abs(x(i)) + epsilon)
% Inputs:
% rwt_mode
% iteration number
% varargin - x : running signal estimate
%            M : number of measurements

switch rwt_mode
    case 1
        % Increase alpha and beta in such a way that beta > alpha > 1
        % This way weights decrease on all the indices, but they decrease
        % at a faster rate on the active set. 
        alpha = sqrt(itr); beta = z; epsilon = 1;
    case 2
        st_wt = 2; alpha = st_wt+itr; beta = 2.5*st_wt+itr; epsilon = 0.1;
    case 3
        st_wt = 5; alpha = st_wt+itr; beta = 2*st_wt+itr; epsilon = 0.1;
    case 4
        % Fixed alpha and beta (boring)
        alpha = 1; beta = 1; epsilon = 0.1;
    case 5
        % Increase weights if solution becomes dense (L1/L2 norm as a proxy
        % for the support size, note that \|x\|_2 \le \|x\|_1 \|
        % sqrt(K)\|x\|_2
        % Decrease the weights when M is large. 
        % Update beta at every reweighting iteration according to the
        % solution. 
        x = varargin{1};
        M = varargin{2};
        alpha = 1; beta = M*(norm(x,2)/norm(x,1))^2; epsilon = 1;
end