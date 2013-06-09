% function lhs = mvprod(A,rhs,Gamma,flag)
%   Usage: lhs = mvprod(A,rhs,Gamma,flag);
%   flag = 0: lhs = A_Gamma * rhs
%   otherwise: lhs = A_Gamma' * rhs 
%
%  Matrix vector multiplication using indices of the matrix.
%  Written by: Salman Asif, Georgia Tech
