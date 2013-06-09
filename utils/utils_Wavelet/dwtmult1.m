% dwtmult1.m
%
% Performs multiple levels of the discrete wavelet transform (using
% dwtlevel1).  Lines up stuff for "tree" structure
% Usage : w = dwtmult1(x, h0, h1, L)
%
% Written by : Justin Romberg
% Created : 5/1/2001

function w = dwtmult1(x, h0, h1, L, sym)

if (nargin == 4), sym = 0; end

N = length(x);

w = x;
for ll = 1:L
  w(1:N*2^(-ll+1)) = dwtlevel1(w(1:N*2^(-ll+1)), h0, h1, sym);
end
