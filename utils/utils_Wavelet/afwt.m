% afwt.m
%
% Adjoint of Forward Fast wavelet transform.
%
% Usage: x = afwt(w, h0, h1, J, sym);
% w - nx1 wavelet transoform vector
% h0 - lowpass decomposition filter
% h1 - highpass decomposition filter
%      h0 and h1 should be zero padded appropriately so that they have the
%      same length, and this length is even.
% J - number of levels in the filter bank.
% 
% sym - How to extend the input. (currently using only sym=1)
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
% Written by: M. Salman Asif
% Created: December 2007

