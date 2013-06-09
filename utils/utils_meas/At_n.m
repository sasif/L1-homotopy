function x = At_n(y, N, OMEGA,P) 
% Takes transpose of noiselet measurements.
% Usage: x = At_n(y)
% x - N vector
% y - M vector 
% OMEGA - M-length vector denoting which noiselet coefficients to use
% P - Column permutation matrix 

n = sqrt(N);

x_noise = zeros(N,1);
x_noise(OMEGA) = y;

x_inoise = realnoiselet(x_noise)/n;

x = zeros(N,1);
x(P) = x_inoise(:);