% symext.m
%
% Symmetrically extends the input vector
% Usage : xe = symext(x, ln, rn, par, sym)
% x - vector to be extended, 1xN
% ln - number of values to tack onto the left of x
% rn - number of values to tack onto the right of x
% par - parity of the extension, 
%       'ee' or 'e' - repeating ext., 
%       'oo' or 'o' non-repeating,
%       'eo' - (2,1) extension (even on right, odd on left)
%       'oe' - (1,2) extension (odd on right, even on left)
%       
% sym - if par = 'e', sym gives whether we have a symmetric extension
% (sym='s') or asymmetric (sym='a')
%
% Written by : Justin Romberg
% Created : 9/4/2001

function xe = symext(x, ln, rn, par, sym)

if (nargin < 5)
  sym = 's';
end

N = length(x);
t = [fliplr(1-(1:ln)) 1:N N+1:N+rn];
if ((par == 'o') | (par == 'oo'))
  yl = 1;
  yh = N;
elseif ((par == 'e') | (par == 'ee'))
  yl = 0.5;
  yh = N+0.5;
elseif (par == 'eo')
  yl = 0.5;
  yh = N;
elseif (par == 'oe')
  yl = 1;
  yh = N+0.5;
end

if (abs(yl-yh)<1)
  tr = ones(1,ln+rn+1);
else
  tr = zeros(size(t));
  p = mod(floor((t-yl)/(yh-yl)),2);
  i0 = find(p == 0);
  tr(i0) = mod(t(i0)-yl, yh-yl) + yl;
  i1 = find(p == 1);
  tr(i1) = yh - mod(t(i1)-yl, yh-yl);
end

tr = round(tr);
if (sym == 'a')
  xe = (1-2*p).*x(tr);
elseif (sym == 's')
  xe = x(tr);
end
