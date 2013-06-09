% cshift.m
%
% Circular shift of a row vector.
% Usage : y = cshift(x, t, dir)
% x - input vector
% t - number of spots to shift
% dir - either 'r' or 'l' (default is 'r')
%
% Written by : Justin Romberg

function y = cshift(x, t, dir)

if (nargin == 2), dir = 'r'; end

N = size(x,2);
t = mod(t,N);
if (dir == 'r')
  y = [x(:,N-t+1:N) x(:,1:N-t)];
else
  y = [x(:,1+t:N) x(:,1:t)];
end



%N = length(x);
%t = mod(t, N);
%if (dir == 'r')
%  y = [x(N-t+1:N) x(1:N-t)];
%else
%  y = [x(1+t:N) x(1:t)];
%end
