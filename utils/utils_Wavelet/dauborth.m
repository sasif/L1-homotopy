% dauborth.m
%
% Filter coefficients for Daubechies' compactly supported wavelets.
% Usage: [h0,h1,g0,g1] = dauborth(M)
% h0 - analysis scaling filter
% h1 - analysis wavelet filter, M/2 vanishing moments
% g0 - synthesis scaling filter
% g1 - synthesis wavelet filter
%
% Written by: Justin Romberg
% Created: 3/23/2004

function [h0,h1,g0,g1] = dauborth(M)

if (mod(M,2) ~= 0)
  error('M must be even.')
end

[P,Q] = daubpoly(M/2);

qr = sort(roots(Q));
hr = [-ones(M/2,1); qr(1:M/2-1)];
gr = [-ones(M/2,1); qr(M/2:M-2)];

h0p = poly(hr);
g0p = poly(gr);
h0 = sqrt(2)*h0p/sum(h0p);
g0 = sqrt(2)*g0p/sum(g0p);
h1 = (-1).^(1:M).*g0;
g1 = (-1).^(0:M-1).*h0;

% flip h0 and g0 to produce functions in the text
% h0 = fliplr(h0); 
% g0 = fliplr(g0);
% h1 = (-1).^(1:M).*g0;
% g1 = (-1).^(0:M-1).*h0;

% another source: http://wavelets.pybytes.com/wavelet/db4/
