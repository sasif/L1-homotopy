function a = binom(n,k)
%
% a = binom(n,k)
% BINOMIAL COEFFICIENTS
%
% allowable inputs:
%       n : integer, k : integer
%       n : integer vector, k : integer
%       n : integer, k : integer vector
%       n : integer vector, k : integer vector (of equal dimension)
%

% Ivan Selesnick
% selesi@taco.poly.edu
% Polytechnic University

nv = n;
kv = k;
if (length(nv) == 1) & (length(kv) > 1)
        nv = nv * ones(size(kv));
elseif (length(nv) > 1) & (length(kv) == 1)
        kv = kv * ones(size(nv));
end
a = nv;
for i = 1:length(nv)
   n = nv(i);
   k = kv(i);
   if n >= 0
      if k >= 0
         if n >= k
            c = prod(1:n)/(prod(1:k)*prod(1:n-k));
         else
            c = 0;
        end
     else
        c = 0;
     end
   else
      if k >= 0
         c = (-1)^k * prod(1:k-n-1)/(prod(1:k)*prod(1:-n-1));
      else
         if n >= k
            c = (-1)^(n-k)*prod(1:-k-1)/(prod(1:n-k)*prod(1:-n-1));
         else
            c = 0;
         end
      end
   end
   a(i) = c;
end


