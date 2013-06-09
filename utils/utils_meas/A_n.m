function y = A_n(x,OMEGA,P)
% Takes noiselet measurements.
% Usage: y = A_n(x)
% x - N vector
% y - M vector 
% OMEGA - M-length vector denoting which noiselet coefficients to use
% P - Column permutation matrix 

N = length(x);
n = sqrt(N);

yt = realnoiselet(x(P))/n;

y = yt(OMEGA);