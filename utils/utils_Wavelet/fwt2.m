% fwt2.m
%
% Fast 2D wavelet transform.
%
% Usage: w = fwt(x, h0, h1, J, sym);
% x - mxn input image (min(m,n) must be divisible by 2^J) 
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

%    2D Example's  output and explanation:
%
%       The coefficients in w are arranged as follows.
%
%              .------------------.
%              |         |        |
%              |    0    |   1    |
%              |         |        |
%              |   L,L   |   H,L  |
%              |         |        |
%              --------------------
%              |         |        |
%              |    2    |   3    |
%              |         |        |
%              |   L,H   |  H,H   |
%              |         |        |
%              `------------------'
%       
%       where 
%            0 : Low pass vertically and Low pass horizontally 
%                (scaling coefficients)
%            1 : Low pass vertically and high pass horizontally
%            2 : High pass vertically and low  pass horizontally
%            3 : High pass vertically and high pass horizontally
%            
%
% Written by: Justin Romberg
% Created: April 2007

% Modified by: Salman Asif
%              December 2011
% Added the ability to work with non-square images. 
% The only requirement is that min(ROW,COL) should be dyadic upto the
% chosen scale J