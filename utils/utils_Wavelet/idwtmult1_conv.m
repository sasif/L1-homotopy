% idwtmult1_conv.m
%
% adjoint of dwtmult1_conv.

function x = idwtmult1_conv(w, g0, g1, J)

sym = 3;

w = w(:)';

L = length(w);
x = [];
xl = w(1:L*2^(-J+1));
for j = J:-1:1
  xh = idwtlevel1(xl, g0, g1, 3);
  xh = xh(1:end-1);
  if j == 1
      x = xh; 
      return;
  end
  xp = w(L*2^(-j+1)+1:L*2^(-j+2));
  npad = length(xh)-length(xp);
  
  % center 
  % xl = [xh zeros(1,ceil(npad/2)) xp zeros(1,floor(npad/2))]; 
  
  % left 
  xl = [xh xp zeros(1,npad)];   
end
