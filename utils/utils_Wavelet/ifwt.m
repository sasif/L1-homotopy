% ifwt.m
%
% Inverse fast wavelet transform.
%
% Usage: x = ifwt(w, g0, g1, J, sym);
% w - nx1 input vector of wavelet coeffients
% g0 - lowpass reconstruction filter
% g1 - highpass reconstruction filter
%      g0 and g1 should be zero padded appropriately so that they have the
%      same length, and this length is even.
% J - number of levels in the filter bank.
%     Right now, it must be chose so that n*2^(-J) >= length(g0)
% sym - How the input was extended when the forward transform was taken.
%     sym=0: periodization
%     sym=1: type-I symmetric extension ( [... x(2) x(1) x(2) x(3) ...])
%            The wavelet filters must have type-I even symmetry 
%            (e.g. daub79)
%     sym=2: type-II symmetric extension ( [... x(2) x(1) x(1) x(2) x(3) ...])
%            The lowpass filter must have type-II even symmetry, 
%            The highpass filter must have type-II odd symmetry.
%            (e.g. daub1018)
%
% Written by: Justin Romberg
% Created: April 2007