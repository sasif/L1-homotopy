function PHI = create_PHI(in);

% Inputs
%   P -- desired number of sampling windowsLOT window
%   N -- length of each sampling window
%   M -- number of measuremnets per sampling window
%   mType -- type of measurements (randn, sign)
%   SM -- sampling mode (starting and ending position of each row)
%       1 : Project sampling waveforms on the LOT basis 
%           Sample each window separately
%           Resultant system will be block diagonal with LOT coeffs
%           measureed with mType waveforms
%       2 : Project sampling waveform on the LOT basis
%           Sample multiple windows in each measurement
%           Resultant system is banded 
%       3 : General sampling system
%           Sample signal directly with mType samples
%           Resultant system will be banded with LOT coeffs. measured with
%           mType*PSI waveforms... 
%       4 : Overlapping samples
%           Applicable with compressive interference cancellation
%       5 : A windowed measurements and reconstructed with DCT-I basis at 
%           right-side border

% For SM = 1, 
% In the LOT coefficients domain, the final system matrix will be block 
% diganal, where each block measures only one set of LOT coefficints. 
% In the signal domain, the measurement waveforms are created by first 
% projecting random (Gaussian or Bernoulli) vectors onto LOT bases... 
% The sensing matrix rows will overlap and align with the boundary of each 
% LOT window boundaries 
% LM = N+2*eta*N; length of sensing waveform
% Note: Each set of LOT coefficients can be estimated separately.

% For SM = 2, (multiple LOT coefficients measured together--streaming setup)
% In the LOT coeffs domain, overlapping measurements of LOT coefficients 
% in d windows 
% LM = d*N+2*eta*N; length of sensing waveform in the signal domain
% Note: Need streaming setup here
        
% For SM = 3; 
% 'universal' sampling scheme 
% It might be better to align the measurements such that they overlap with
% a few LOT windows before the last overlapping interval
% LM = d*N; length of sensing waveforms overlapping d LOT windows

% Sensing matrix for each
genAmat_h = @(M,N) genAmat(M,N,in);

P = in.P; M = in.M; N = in.N; LM = in.LM;
PHI = zeros(P*M,(P-1)*N+LM);
for p = 1:P
    PHI((p-1)*M+1:p*M,(p-1)*N+1:(p-1)*N+LM) = genAmat_h(M,LM);
end
