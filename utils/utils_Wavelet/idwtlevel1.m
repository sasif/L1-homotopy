% idwtlevel1.m
%
% Inverse of one level of the discrete wavelet transform (dwtlevel1).
%
% Written by : Justin Romberg
% Created : 5/1/2001

function x = idwtlevel1(w, g0, g1, sym)

if (nargin == 3), sym = 0; end

N = length(w);
m0 = length(g0);
m1 = length(g1);

if (m0 ~= m1)
  error('Use biorfilt to create filters');
end


% NEED TO FIX THIS SO IT WORKS LIKE BIORFILT()
% if the filters are uneven lengths, zero pad them
%if (m0 > m1)
%  g1 = [zeros(1,(m0-m1)/2) g1 zeros(1,(m0-m1)/2)];
%  m1 = length(g1);
%elseif (m1 > m0)
%  g0 = [zeros(1,(m1-m0)/2) g0 zeros(1,(m1-m0)/2)];
%  m0 = length(g0);
%end

% center on the middle of the wavelet filter
% odd length - hole at (m1-1)/2
% even length - hole at m1/2
k = floor(m1/2)-1;
if (sym == 0)
  % upsample
  c0 = reshape([w(1:N/2); zeros(1,N/2)], [1 N]);
  c1 = reshape([w(N/2+1:N); zeros(1,N/2)], [1 N]);
  x = cconv(c0, g0, k) + cconv(c1, g1, k);
elseif (sym == 1)
  k = floor(m1/2);
  % symmetrically extend
  c0e = symext(w(1:N/2), k, m1-1-k, 'oe', 's');
  c1e = symext(w(N/2+1:N), k, m1-1-k, 'eo', 's');
  % upsample
  c0 = reshape([c0e; zeros(1,N/2+m1-1)], [1 N+2*(m1-1)]);
  c1 = reshape([c1e; zeros(1,N/2+m1-1)], [1 N+2*(m1-1)]);
  % filter
  x = conv2(c0(k+1:N+k+m0-1), g0, 'valid') + ...
      conv2(c1(k+1:N+k+m0-1), g1, 'valid');
elseif (sym == 2) 
  % symmetrically extend
  c0e = symext(w(1:N/2), m1-1-k, k, 'e', 's');
  c1e = symext(w(N/2+1:N), m1-1-k, k, 'e', 'a');
  % upsample
  c0 = reshape([c0e; zeros(1,N/2+m1-1)], [1 N+2*(m1-1)]);
  c1 = reshape([c1e; zeros(1,N/2+m1-1)], [1 N+2*(m1-1)]);
  % filter
  x = conv2(c0(m0-k:2*(m1-1)-k+N), g0, 'valid') + ...
      conv2(c1(m0-k:2*(m1-1)-k+N), g1, 'valid');
elseif (sym == 3)
  % linear filtering (no extension)
  c0 = reshape([w(1:length(w)/2); zeros(1,length(w)/2)], [1 length(w)]);
  c1 = reshape([w(length(w)/2+1:end); zeros(1,length(w)/2)], [1 length(w)]);
  x = conv(c0, g0) + conv(c1, g1);
end

