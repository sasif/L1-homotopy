% daubpoly.m
%
% Daubechies polynomial for compactly supported wavelets.
% Factor the answer in different ways to get differnet filter coeffs.
% Usage : [P, Q] = daubpoly(p)
% Q - polynomial for roots at locations other than -1 (for PR)
% P - polynomial conv(q,(1+z)^(2*p))
% Ex:  
%    let qr = roots(Q), pr = roots(P)
%    set h0 = poly([-ones(1,p); qr(1:p-1)])
%        g0 = poly([-ones(1,p); qr(p:2*p-2)])
%    for Daub orthonormal basis
%
% Written by : Justin Romberg
% Created : 3/23/2004

function [P, Q] = daubpoly(p)

B = binom(2*p,0:2*p);
qm = zeros(2*p-1);
B2 = [zeros(1,2*p-3) B zeros(1,2*p-3)];
for kk = 1:2*p-1
  qm(kk,:) = fliplr(B2(2*(kk-1)+1:2*(kk-1)+2*p-1));
end
Q = inv(qm)*[zeros(1,p-1) 1 zeros(1,p-1)]';
P = conv(Q,B);
