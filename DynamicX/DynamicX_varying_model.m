function sig_mod = DynamicX_varying_model(Name,n);

% I have 'perturbed' the following function from WaveLab.
%
% MakeSignal -- Make artificial signal
%  Usage
%    sig = MakeSignal(Name,n)
%  Inputs
%    Name   string: 'HeaviSine', 'Bumps', 'Blocks',
%            'Doppler', 'Ramp', 'Cusp', 'Sing', 'HiSine',
%            'LoSine', 'LinChirp', 'TwoChirp', 'QuadChirp',
%            'MishMash', 'WernerSorrows' (Heisenberg),
%            'Leopold' (Kronecker), 'Piece-Regular' (Piece-Wise Smooth),
%	     'Riemann','HypChirps','LinChirps', 'Chirps', 'Gabor'
%	     'sineoneoverx','Cusp2','SmoothCusp','Gaussian'
%	     'Piece-Polynomial' (Piece-Wise 3rd degree polynomial)
%    n      desired signal length
%  Outputs
%    sig    1-d signal
%
%  References
%    Various articles of D.L. Donoho and I.M. Johnstone

% N = 2048;
%
% HeaviSine = MakeSignal('HeaviSine',N);
% Blocks    = MakeSignal('Blocks',N);
% Doppler   = MakeSignal('Doppler',N);
% Ramp      = MakeSignal('Ramp',N);
% Piece_Polynomial = MakeSignal('Piece-Polynomial',N);
%
% name = 1;
% switch name
%     case 1
%         Name = 'HeaviSine';
%     case 2
%         Name = 'Blocks';
%     case 3
%         Name = 'Doppler';
%     case 4
%         Name = 'Ramp';
%     case 5
%         Name = 'Piece-Polynomial';
%     otherwise
%         disp('hey ... ');
% end

N = n;
sig = MakeSignal(Name,N);

sig_old = sig;

% modify the signal randomly while preserving the structure.
t = (1:N)./N;
if strcmp(Name,'HeaviSine'),
    signal = (1+randn*.05)*4.*sin(4*pi.*t*(1+randn*.02));
    signal = signal - sign(t - .3*(1+randn*.01)) - sign(.72*(1+randn*.01) - t);

elseif strcmp(Name,'Doppler'),
    sig = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
    signal = (1+randn*.02)*sqrt(t.*(1-t)).*sin((2*pi*1.05*(1+randn*.02)) ./(t+.05*(1+randn*.02)));

elseif strcmp(Name,'Ramp')
    sig = t - (t >= .37);
    signal = (1+randn*.05)*t - (t >= .37*(1+randn*.02));
elseif strcmp (Name,'Blocks')
    pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
    sig = zeros(size(t));
    for j=1:length(pos)
        sig = sig + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
    end
    signal = sig;
    for j=1:length(pos)-1
        rand_s = randn;
        wt = (1+rand_s*.1);
        wt = 1+randsrc*rand*.2;
        signal(ceil(pos(j)*N):floor(pos(j+1)*N)) =  signal(ceil(pos(j)*N):floor(pos(j+1)*N))*wt;
    end

elseif strcmp(Name,'Piece-Polynomial'),
    n = N;
    t = (1:fix(n/5)) ./fix(n/5);
    sig1=20*(t.^3+t.^2+4)*(1+(randn*.1));
    sig3=40*(2.*t.^3+t)*(1+(randn*.1)) + 100*(1+(randn*.05));
    sig2=10.*t.^3*(1+(randn*.05)) + 45*(1+(randn*.1));
    sig4=16*t.^2+8.*t*(1+randn*.1)+16*(1+(randn*.05));
    sig5=20*(t+4)*(1+(randn*.1));
    sig6(1:fix(n/10))=ones(1,fix(n/10))*(1+(randn*.1));
    sig6=sig6*20*(1+(randn*.05));
    sig(1:fix(n/5))=sig1;
    sig(2*fix(n/5):-1:(fix(n/5)+1))=sig2;
    sig((2*fix(n/5)+1):3*fix(n/5))=sig3;
    sig((3*fix(n/5)+1):4*fix(n/5))=sig4;
    sig((4*fix(n/5)+1):5*fix(n/5))=sig5(fix(n/5):-1:1);
    diff=n-5*fix(n/5);
    sig(5*fix(n/5)+1:n)=sig(diff:-1:1);
    %sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=-ones(1,fix(n/10))*20;
    sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=ones(1,fix(n/10))*10;
    sig((n-fix(n/10)+1):(n+fix(n/20)-fix(n/10)))=ones(1,fix(n/20))*150;
    % zero-mean
    bias=sum(sig)/n;
    sig=sig-bias;
    signal = sig;
end

sig_mod = signal;

% Originally made by David L. Donoho.
% Function has been enhanced.

%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:39 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu