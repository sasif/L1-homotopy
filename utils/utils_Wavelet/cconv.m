% cconv.m
%
% Circular convolution.
% Usage : y = cconv(x,h,k)
% x - input signal
% h - filter
% k - output hole in filter (k=0 if the filter is causal)
%
% Written by : Justin Romberg
% Created : 2/16/98, Revised : 5/1/2001

function y = cconv(x,h,k)

if (nargin == 2), k = 0; end

N = length(x);
M = length(h);
inds = mod([1:N 1:M-1]-1, N) + 1;
y = cshift(conv2(x(inds), h, 'valid'),M-k-1,'r');

% old code
%
%inds = mod([N-M+2+k:N 1:N 1:k]-1, N) + 1;
%y = conv2(x(inds), h, 'valid');
