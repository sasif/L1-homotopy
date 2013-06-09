% fwt.m
%
% Fast wavelet transform.
%
% Usage: w = fwt(x, h0, h1, J, sym);
% x - nx1 input vector
% h0 - lowpass decomposition filter
% h1 - highpass decomposition filter
%      h0 and h1 should be zero padded appropriately so that they have the
%      same length, and this length is even.
% J - number of levels in the filter bank.
%     Right now, it must be chose so that n*2^(-J) >= length(h0)
% sym - How to extend the input.
%     sym=0: periodization
%     sym=1: type-I symmetric extension ( [... x(2) x(1) x(2) x(3) ...])
%            The wavelet filters must have type-I even symmetry 
%            (e.g. daub79)
%     sym=2: type-II symmetric extension ( [... x(2) x(1) x(1) x(2) x(3) ...])
%            The lowpass filter must have type-II even symmetry, 
%            The highpass filter must have type-II odd symmetry.
%            (e.g. daub1018)
%
% A friendly way to understand/implement symmetric wavelets and adjoints 
% is to treat symmetric extension and filtering as two separate operations. 
%
% Written by: Justin Romberg
% Created: April 2007

