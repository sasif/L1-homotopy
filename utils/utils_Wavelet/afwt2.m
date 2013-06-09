% afwt2.m
%
% Adjoint of fwt2
%
% Usage: x = afwt2(w, h0, h1, J, sym);
% x - mxn output image (min(m,n) must be divisible by 2^J) 
% w - mxn input image of wavelet coeffients
% h0 - lowpass decomposition filter
% h1 - highpass decomposition filter
%      h0 and h1 should be zero padded appropriately so that they have the
%      same length, and this length is even.
% J - number of levels in the filter bank.
%     Right now, it must be chosen so that n*2^(-J) >= length(h0)
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

%
% Written by: Justin Romberg
% Created: April 2007