% idwtmult1.m
%
% Inverts dwtmult1.
%
% Written by : Justin Romberg
% Created : 5/1/2001

function x = idwtmult1(w, g0, g1, L, sym)

if (nargin == 4), sym=0; end

N = length(w);

for ll = L:-1:1
  w(1:N*2^(-ll+1)) = idwtlevel1(w(1:N*2^(-ll+1)), g0, g1, sym);
end
x = w;
  