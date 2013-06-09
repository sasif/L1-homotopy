% function iAtB_mod = update_inverse(AtB, Atb, atB, atb);

% This function uses matrix inversion lemma to update the inverse of 
% square matrix after addition or removal of a row-column pair.
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu

function iAtB = update_inverse(AtB, iAtB_old,flag);

n = size(AtB,1);

%A11 = AtB(1:n-1,1:n-1);
A12 = AtB(1:n-1,n);
A21 = AtB(n,1:n-1);
A22 = AtB(n,n);

% add columns
if flag == 1
    iA11 = iAtB_old;
    iA11A12 = iA11*A12;
    A21iA11 = A21*iA11;
    S = A22-A21*iA11A12;
    Q11_right = iA11A12*(A21iA11/S); 
%     Q11 = iA11+ Q11_right;
%     Q12 = -iA11A12/S;
%     Q21 = -A21iA11/S;
%     Q22 = 1/S;

    iAtB = zeros(n);
    %iAtB = [Q11 Q12; Q21 Q22]; 
    iAtB(1:n-1,1:n-1) = iA11+ Q11_right;
    iAtB(1:n-1,n) = -iA11A12/S; 
    iAtB(n,1:n-1) =  -A21iA11/S;
    iAtB(n,n) = 1/S;
%delete columns
else if flag == 2
        Q11 = iAtB_old(1:n-1,1:n-1);
        Q12 = iAtB_old(1:n-1,n);
        Q21 = iAtB_old(n,1:n-1);
        Q22 = iAtB_old(n,n);
        Q12Q21_Q22 = Q12*(Q21/Q22);
        iAtB = Q11 - Q12Q21_Q22;
        %iAtB = iA11;
    end
end
