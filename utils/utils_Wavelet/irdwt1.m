% irdwt1.m
%
% Inverse redundant discrete wavelet transform.
% Inverts "rdwt1.m"
% Usage : x = irdwt1(w, g0, g1, nlev)
% w - rdwt coefficients, nlev+1xN
% g0,g1 - reconstruction filters
% nlev - number of levels in the decomposition
%
% Written by : Justin Romberg
% Created : 8/10/2001

function x = irdwt1(w, g0, g1, nlev)

N = size(w,2);

w = flipud(w);
w(nlev+1,:) = fliplr(cshift(w(nlev+1,:), 1, 'l'));
m0 = length(g0);
m1 = length(g1);
g0i = reshape([zeros(2^(nlev-1)-1, m0); g0], [1 m0*2^(nlev-1)]);
g1i = reshape([zeros(2^(nlev-1)-1, m1); g1], [1 m1*2^(nlev-1)]);
for ll = nlev:-1:1
  m0i = length(g0i);
  m1i = length(g1i);
  zi = floor((m1i-m0i)/4)+1;
  si = cshift(w(ll+1,:), zi, 'r');
  di = fliplr(cshift(w(ll,:), 1, 'l'));
  w(ll,:) = (cconv(si, g0i, floor(m1i/2)) + cconv(di, g1i, floor(m1i/2)))/2;
  g0i = g0i(2:2:end);
  g1i = g1i(2:2:end);
end
x = w(1,:);
