% cdfspline.m
%
% Filter coefficients for Cohen, Daubechies, Faveau biorth spline filters
% Usage : [h0,g0] = cdfspline(N,M)
% h0 - bspline poly of order N
% g0 - poly using rest of factors from daubpoly((N+M)/2)
% h0 and g0 are scaling filters
%
% Written by : Justin Romberg
% Created : 3/23/2004

function [h0, g0] = cdfspline(N, M)


if (mod(N+M,2) ~= 0)
  error('N+M must be even.')
end


[P, Q] = daubpoly((N+M)/2);
qr = roots(Q);

h0r = -ones(N,1); 
g0r = [-ones(M,1); qr];

h0p = poly(h0r);
g0p = poly(g0r);
h0 = [zeros(1,M-1) sqrt(2)*h0p/sum(h0p) zeros(1,M-1)];
g0 = sqrt(2)*g0p/sum(g0p);


