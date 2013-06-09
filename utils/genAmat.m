function A = genAmat(M,N,in);

switch in.type
    case 'randn'
        % Gaussian measurements
        A = (randn(M,N))/sqrt(M);
        % A = orth(A')';   % with orthogonal rows
    case 'hadamard'
        % Hadamard
        H = hadamard(N);
        q = randperm(N);
        A = H(q(1:M),:)/sqrt(M);
    case 'sign'
        % Bernoulli
        A = sign(randn(M,N))/sqrt(M);
    case 'orth'
        % Random Projection
        A = (randn(M,N))/sqrt(M);
        A = orth(A')';
    case 'rdct'
        % Randomly subsampled DCT
        D = dct(eye(N));
        q = randperm(N);
        A = D(q(1:M),:);
    case 'inc'
    case {'noiselets','rsample'}
        q = randperm(N);
        OMEGA = q(1:M);
        P = 1:N; % randperm(N);
        P = P(:); OMEGA = OMEGA(:);
        A = @(x,mode) Afun(x,mode,N, OMEGA, P,in.type);
    case 'subsample'
        q = randperm(N);
        OMEGA = q(1:M);
        A = eye(N);
        A = A(OMEGA,:);
    case 'streaming'        
        parts = in.parts;
        Mp = round(M/parts);
        Np = round(N/parts);
        A = zeros(M,N);
        for ii = 1:parts
            A((ii-1)*Mp+1:ii*Mp,(ii-1)*Np+1:ii*Np) = sign(randn(Mp,Np))/sqrt(Mp);
        end
    case 'rnd-demod'
    case 'rnd-conv'
end

function y = Afun(x,mode, N, OMEGA,P,mtype)
if mode == 1
    switch mtype
        case 'noiselets'
            y = A_n(x,OMEGA,P);
        case 'rsample'
            y = x(OMEGA);
    end
end
if mode == 2
    switch mtype
        case 'noiselets'
            y = At_n(x, N, OMEGA,P);
        case 'rsample'
            y = zeros(N,1);
            y(OMEGA) = x;
    end
end