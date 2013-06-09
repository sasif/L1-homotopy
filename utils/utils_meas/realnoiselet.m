%    Real-valued Noiselet Transform
% 
%    Usage: w = realnoiselet(x)
%    x must be a REAL VALUED COLUMN VECTOR or MATRIX
%    m = size(x,1) must be a POWER OF TWO
% 
%    Notes:
%    1) This implementation uses exactly m*log2(m) additions/subtractions.
%    2) This is a real-valued variation of "dragon" noiselets.
%    3) This is symmetric and orthogonal. To invert, apply again and
%       divide by m.
% 
%    Written by: Justin Romberg, Georgia Tech
%                Peter Stobbe, Caltech
%    Email: jrom@ece.gatech.edu, stobbe@acm.caltech.edu
%    Created: October 2006
%    Last Modified: December 2006
