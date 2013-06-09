% dwtlevel1.m
%
% One level of the discrete wavelet transform.  Periodic extension (for now)
% and the wavelets are lined up the best way for "tree structure".
% Usage : w = dwtlevel1(x, h0, h1, sym)
%
% Written by : Justin Romberg
% Created : 5/1/2001

function w = dwtlevel1(x, h0, h1, sym)

if (nargin == 3), sym=0; end

N = length(x);
m0 = length(h0);
m1 = length(h1);

if (m0 ~= m1)
  error('Use biorfilt to create filters');
end

% NEED TO FIX THIS SO IT WORKS LIKE BIORFILT()
% if wavelet filters have unequal length, zero pad the shorter one
% the difference in length has to be even
%if (m0 > m1)
%  h1 = [zeros(1,(m0-m1)/2) h1 zeros(1,(m0-m1)/2)];
%  m1 = length(h1);
%elseif (m1 > m0)
%  h0 = [zeros(1,(m1-m0)/2) h0 zeros(1,(m1-m0)/2)];
%  m0 = length(h0);
%end

% center on the middle of the filter
% odd length - hole at (m1-1)/2 
%    (bottom left node of tree centers on sample 1.5)
% even length - hole at m1/2
%    (bottom left node of tree centers on sample 1)
k = floor(m1/2);
if (sym == 0)
  % circular filtering
  c0 = cconv(x, h0, k);
  c1 = cconv(x, h1, k);
  % downsample
  w = [c0(1:2:N) c1(1:2:N)];
elseif (sym == 1)
  k = floor(m1/2)-1;
  % symmetrically extend
  xe = symext(x, k, m0-k-1, 'o', 's');
  % filter
  c0 = conv2(xe, h0, 'valid');
  c1 = conv2(xe, h1, 'valid');
  % downsample
  w = [c0(1:2:N) c1(1:2:N)];
elseif (sym == 2)
  % symmetrically extend
  xe = symext(x, m0-1-k, k, 'e', 's');
  % filter
  c0 = conv2(xe, h0, 'valid');
  c1 = conv2(xe, h1, 'valid');
  % downsample
  w = [c0(1:2:N) c1(1:2:N)];
elseif (sym == 3)
  % linear filtering
  c0 = conv(x, h0);
  c1 = conv(x, h1);
  % downsample
  w = [c0(1:2:end); c1(1:2:end)];
end

