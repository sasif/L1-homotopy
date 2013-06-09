function out = OMP_function(y,A,in)

[M N] = size(A);
Te = M;
if isfield(in,'Te')
    Te = in.Te;
end
thresh = 1e-5;
if isfield(in,'thresh')
    thresh = in.thresh;
end

[val_max gamma_omp] = max(abs(A'*y));

for iter = 1:Te
    x_omp = zeros(N,1);
    x_omp(gamma_omp) = A(:,gamma_omp)\y;
    r_omp = y-A*x_omp;
    if norm(r_omp)< thresh
        break;
    end
    p_omp = A'*r_omp;
    gamma_ompC = setdiff([1:N],gamma_omp);
    [val_omp, ind_omp] = max(abs(p_omp(gamma_ompC)));
    gamma_omp = [gamma_omp; gamma_ompC(ind_omp)];
end
out = [];
out.x_out = x_omp;
out.gamma = gamma_omp;
out.iter = iter;