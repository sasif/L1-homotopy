% dwtmult2.m
%
% Performs multiple levels of the 2D discrete wavelet transform.
% Usage : w = dwtmult2(x, h0, h1, L, sym)
%
% Written by : Justin Romberg
% Created : 6/26/2001

function w = dwtmult2(x, h0, h1, L, sym)

if (nargin == 4), sym = 0; end

[Nr,Nc] = size(x);

w = x;
for ll = 1:L
  % rows
  for ii = 1:(Nr*2^(-ll+1))
    w(ii,1:Nc*2^(-ll+1)) = dwtlevel1(w(ii,1:Nc*2^(-ll+1)), h0, h1, sym);
  end
  % columns
  for jj = 1:(Nc*2^(-ll+1))
    w(1:Nr*2^(-ll+1),jj) = dwtlevel1(w(1:Nr*2^(-ll+1),jj)', h0, h1, sym)';
  end
end

