function [x varargout] = genSignal(N,in)

%  Inputs
%   N -- length of the signal
%   in (input structure)
%       type    -- signal type (details below)
%       randgen -- randomize the signal (1) or reproduce the signal from
%                   Wavelab (0)
%       wType   -- type of wavelet transform
%
%   sType: 
%
%   T-sparse signals: randn, sign, ones, highD
%
%   The following functions adopted from Wavelab toolbox
%   Name string: 'HeaviSine', 'Blocks', 'Bumps'
%            'Doppler', 'Ramp', 'Cusp', 'Sing',            
%            'pcwPoly' (Piece-Wise 3rd degree polynomial)
%            'pcwreg' (Piece-Regular, Piece-Wise Smooth),
%            'LinChirp', 'TwoChirp', 'QuadChirp',
%            'MishMash', 
%            'Cusp2', 'HypChirps','LinChirps', 'Chirps',
%
%   Sounds: 
%       MATLAB options: 'chirp','gong','handel'
%       Wavelab options: 'Greasy','Tweet'
%
%   Images: 
%
%   References
%    Various articles of D.L. Donoho and I.M. Johnstone
%
% Originally made by David L. Donoho.
% Function has been enhanced.
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:39 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m

t = (1:N)./N;
if isfield(in,'randgen'); randgen = in.randgen; else randgen = 0; end
if isfield(in,'take_fwt'); take_fwt = in.take_fwt; else take_fwt = 0; end

wave_struct = [];
sig = [];

sType = in.type;

switch lower(sType)
    case 'randn'
        T = in.T;
        q = randperm(N);
        x = zeros(N,1);
        x(q(1:T)) = randn(T,1);
        sig = x; 
        
    case 'sign'
        T = in.T;
        q = randperm(N);
        x = zeros(N,1);
        x(q(1:T)) = sign(randn(T,1));
        sig = x; 
        
    case 'ones'
        T = in.T;
        q = randperm(N);
        x = zeros(N,1);
        x(q(1:T)) = ones(T,1);
        sig = x; 
        
    case 'highd'
        T = in.T;
        q = randperm(N);
        x = zeros(N,1);
        x(q(1:ceil(T/2))) = sign(randn(ceil(T/2),1));
        x(q(ceil(T/2)+1:T)) = sign(randn(T-ceil(T/2),1))./1e2;
        sig = x; 
        
    %-----------------%
    % Wavelab signals %
    %-----------------%
    case 'blocks'
        if randgen == 0
            pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
            hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
        else
            pos = sort(rand(1,11));
            % hgt = sign(randn(1,10)).*randint(1,10,[1 5]);
            hgt = sign(randn(1,11)).*randi([1 5],1,11);
        end
        sig = zeros(size(t));
        for j=1:length(pos)
            sig = sig + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
        end
        sig = sig(:); 
    case 'bumps'
        if randgen == 0
            pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
            hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
            wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
        else
            pos = sort(rand(1,11));
            hgt = 4+sign(randn(1,11)).*rand(1,11);
            wth = rand(1,11)/100;
        end
        sig = zeros(size(t));
        for j =1:length(pos)
            sig = sig + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
        end
        sig = sig(:); 
        
    case 'pcwpoly' % Piecewise polynomial
        n = N;
        if randgen == 0
            t = (1:fix(n/5)) ./fix(n/5);
            sig1=20*(t.^3+t.^2+4);
            sig3=40*(2.*t.^3+t) + 100;
            sig2=10.*t.^3 + 45;
            sig4=16*t.^2+8.*t+16;
            sig5=20*(t+4);
            sig6(1:fix(n/10))=ones(1,fix(n/10));
            sig6=sig6*20;
            sig(1:fix(n/5))=sig1;
            sig(2*fix(n/5):-1:(fix(n/5)+1))=sig2;
            sig((2*fix(n/5)+1):3*fix(n/5))=sig3;
            sig((3*fix(n/5)+1):4*fix(n/5))=sig4;
            sig((4*fix(n/5)+1):5*fix(n/5))=sig5(fix(n/5):-1:1);
            dif=n-5*fix(n/5);
            sig(5*fix(n/5)+1:n)=sig(dif:-1:1);
            %sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=-ones(1,fix(n/10))*20;
            sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=ones(1,fix(n/10))*10;
            sig((n-fix(n/10)+1):(n+fix(n/20)-fix(n/10)))=ones(1,fix(n/20))*150;
            % zero-mean
            bias=sum(sig)/n;
            sig=sig-bias;
        else
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
            dif=n-5*fix(n/5);
            sig(5*fix(n/5)+1:n)=sig(dif:-1:1);
            %sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=-ones(1,fix(n/10))*20;
            sig((fix(n/20)+1):(fix(n/20)+fix(n/10)))=ones(1,fix(n/10))*10;
            sig((n-fix(n/10)+1):(n+fix(n/20)-fix(n/10)))=ones(1,fix(n/20))*150;
            % zero-mean
            bias=sum(sig)/n;
            sig=sig-bias;
        end
        sig = sig(:);
         
    case 'heavisine'
        if randgen == 0
            sig = 4.*sin(4*pi.*t);
            sig = sig - sign(t - .3) - sign(.72 - t);
        else
            sig = (1+rand*.5)*4.*sin(4*pi.*t*(1+rand*.25));
            sig = sig - sign(t - .3*(1+randn*.1)) - sign(.72*(1+randn*.1) - t)*(1+randn);
        end
        sig = sig(:);
         
    case 'doppler'
        if randgen == 0
            sig = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
        else
            sig = (1+randn*.02)*sqrt(t.*(1-t)).*sin((2*pi*1.05*(1+randn*.02)) ./(t+.05*(1+randn*.02)));
        end
        sig = sig(:);
         
    case 'cusp'
        if randgen == 0
            sig = sqrt(abs(t - .37));
        else
            sig = sqrt(abs(t - (1+0.1*rand)*.37));
        end
        
        sig = sig(:); 
        
    case 'ramp'
        if randgen == 0
            sig = t - (t >= .37);
        else
            sig = (1+randn*.05)*t - (t >= .37*(1+randn*.02));
        end
        sig = sig(:);
              
    case 'pcwreg' % Piecewise regular       
        sig = MakeSignal('Piece-Regular',N);
        
        sig = sig(:);    
        
    case 'cusp2'
        n = N/8;
        i1 = (1:n)./n;
        x = (1-sqrt(i1)) + i1/2 -.5;
        sig = zeros(1,N);
        sig(N-1.5*n+1:N-.5*n) = x;
        sig(N-2.5*n+2:N-1.5*n+1) = fliplr(x);
        sig(3*n+1:3*n + n) = .5*ones(1,n);
        sig = sig(:); 
        
    case 'sing'
        sig = MakeSignal('Sing',N); 

        sig = sig(:);        
        
        
    %--------------%
    % Kalman model %
    %--------------%
    case lower('kalman-heavisine')        
        Np = in.Np;
        t = (1:Np)./Np;
        N = Np*ceil(N/Np);
        sig = zeros(N,1);
 
        randgen = 0;
        for p = 1:ceil(N/Np);
            sigt = zeros(Np,1);
            if randgen == 0
                sigt = 4.*sin(4*pi.*t);
                sigt = sigt - sign(t - .3) - sign(.72 - t);
            else
                sigt = (1+rand*.01)*4.*sin(4*pi.*t*(1+rand*.025));
                sigt = sigt - sign(t - .3*(1+randn*.01)) - sign(.72*(1+randn*.01) - t)*(1+randn*.02);
            end
            randgen = 1;
            sig((p-1)*Np+1:p*Np) = sigt;
        end
        sig = sig(:);
    case 'kalman-pcwpoly' % Piecewise polynomial
        Np = in.Np;
        n = Np;
        N = Np*ceil(N/Np);
        sig = zeros(N,1);
 
        randgen = 0;
        for p = 1:ceil(N/Np);
            sigt = zeros(Np,1);
            
            if randgen == 0                
                t = (1:fix(n/5)) ./fix(n/5);
                sig1=20*(t.^3+t.^2+4);
                sig3=40*(2.*t.^3+t) + 100;
                sig2=10.*t.^3 + 45;
                sig4=16*t.^2+8.*t+16;
                sig5=20*(t+4);
                sig6(1:fix(n/10))=ones(1,fix(n/10));
                sig6=sig6*20;
                sigt(1:fix(n/5))=sig1;
                sigt(2*fix(n/5):-1:(fix(n/5)+1))=sig2;
                sigt((2*fix(n/5)+1):3*fix(n/5))=sig3;
                sigt((3*fix(n/5)+1):4*fix(n/5))=sig4;
                sigt((4*fix(n/5)+1):5*fix(n/5))=sig5(fix(n/5):-1:1);
                dif=n-5*fix(n/5);
                sigt(5*fix(n/5)+1:n)=sigt(dif:-1:1);
                %sigt((fix(n/20)+1):(fix(n/20)+fix(n/10)))=-ones(1,fix(n/10))*20;
                sigt((fix(n/20)+1):(fix(n/20)+fix(n/10)))=ones(1,fix(n/10))*10;
                sigt((n-fix(n/10)+1):(n+fix(n/20)-fix(n/10)))=ones(1,fix(n/20))*150;
                % zero-mean
                bias=sum(sigt)/n;
                sigt=sigt-bias;
            else                
                t = (1:fix(n/5)) ./fix(n/5);
                sig1=20*(t.^3+t.^2+4)*(1+(randn*.1));
                sig3=40*(2.*t.^3+t)*(1+(randn*.1)) + 100*(1+(randn*.05));
                sig2=10.*t.^3*(1+(randn*.05)) + 45*(1+(randn*.1));
                sig4=16*t.^2+8.*t*(1+randn*.1)+16*(1+(randn*.05));
                sig5=20*(t+4)*(1+(randn*.1));
                sig6(1:fix(n/10))=ones(1,fix(n/10))*(1+(randn*.1));
                sig6=sig6*20*(1+(randn*.05));
                sigt(1:fix(n/5))=sig1;
                sigt(2*fix(n/5):-1:(fix(n/5)+1))=sig2;
                sigt((2*fix(n/5)+1):3*fix(n/5))=sig3;
                sigt((3*fix(n/5)+1):4*fix(n/5))=sig4;
                sigt((4*fix(n/5)+1):5*fix(n/5))=sig5(fix(n/5):-1:1);
                dif=n-5*fix(n/5);
                sigt(5*fix(n/5)+1:n)=sigt(dif:-1:1);
                %sigt((fix(n/20)+1):(fix(n/20)+fix(n/10)))=-ones(1,fix(n/10))*20;
                sigt((fix(n/20)+1):(fix(n/20)+fix(n/10)))=ones(1,fix(n/10))*10;
                sigt((n-fix(n/10)+1):(n+fix(n/20)-fix(n/10)))=ones(1,fix(n/20))*150;
                % zero-mean
                bias=sum(sigt)/n;
                sigt=sigt-bias;
            end
            randgen = 0;
            sig((p-1)*Np+1:p*Np) = sigt;
        end
        sig = sig(:);
        
    %--------%
    % images %
    %--------%
    case lower('image-barbara');
        I = double(imread('barbara.ras'));
        [ROW COL] = size(I);
        I = imresize(I, [ROW N]);        
        % I = I-repmat(mean(I,2),1,N);
        vec = @(z) z(:);
        % sig = double(vec(I(zigzag_order(ROW,COL))));
        sig = vec(I'); 
    case lower('image-baboon');
        I = double(imread('baboon.ras'));
        [ROW COL] = size(I);
        I = imresize(I, [ROW N]);        
        % I = I-repmat(mean(I,2),1,N);
        vec = @(z) z(:);
        % sig = double(vec(I(zigzag_order(ROW,COL))));
        sig = vec(I');     
    case lower('image-cameraman');
        I = double(imread('cameraman.tif'));
        [ROW COL] = size(I);
        I = imresize(I, [ROW N]);        
        % I = I-repmat(mean(I,2),1,N);
        vec = @(z) z(:);
        % sig = double(vec(I(zigzag_order(ROW,COL))));
        sig = vec(I');
    case lower('image-house');
        I = double(rgb2gray(imread('house.tif')));
        [ROW COL] = size(I);
        I = imresize(I, [ROW N/2]);        
        % I = I-repmat(mean(I,2),1,N);
        vec = @(z) z(:);
        % sig = double(vec(I(zigzag_order(ROW,COL))));        
        I = [I(1:2:end,:) fliplr(I(2:2:end,:))];
        figure(123);
        subplot(121); imagesc(I');
        subplot(122); imagesc(diff(I',[],2));

        sig = vec(I');
    case lower('image-boats');
        load('boats');
        I = double(sig);
        [ROW COL] = size(I);
        I = imresize(I, [ROW N/2]);        
        % I = I-repmat(mean(I,2),1,N);
        vec = @(z) z(:);
        % sig = double(vec(I(zigzag_order(ROW,COL))));
        % I = [I(1:2:end,:) fliplr(I(2:2:end,:))];        
        figure(123); 
        subplot(121); imagesc(I');
        subplot(122); imagesc(diff(I',[],2));
        
        sig = vec(I');
    case lower('image-saturn')
        I = double(rgb2gray(imread('saturn4.jpg')));
        [ROW COL] = size(I);
        I = imresize(I, [ROW N/2]);        
        % I = I-repmat(mean(I,2),1,N);
        vec = @(z) z(:);
        % sig = double(vec(I(zigzag_order(ROW,COL))));
        I = [I(1:2:end,:) fliplr(I(2:2:end,:))];

        figure(123); 
        subplot(121); imagesc(I');
        subplot(122); imagesc(diff(I',[],2));

        sig = vec(I');     
    case lower('image-chirp')
        I = double(rgb2gray(imread('pattern.jpg')));                 
        [ROW COL] = size(I);
        I = imresize(I, [ROW N]);        
        % I = I-repmat(mean(I,2),1,N);
        vec = @(z) z(:);
        % sig = double(vec(I(zigzag_order(ROW,COL))));
        sig = vec(I');
        % I = double(rgb2gray(imread('ohare1.jpg'))); I = I-mean(I(:)); [M N]=size(I);
        % for ii=1:2:M-1; figure(1); imagesc(I); figure(2); plot(I(ii:ii+1,:)'); figure(3); plot(fftshift(abs(fft(I(ii:ii+1,:)')))); title(ii); drawnow; pause(1/10); end
        
	%--------%
    % Chirps %
    %--------%
    case lower('LinChirp')
        sig = MakeSignal('LinChirp',N);
        sig = sig(:);        
    
    case  lower('TwoChirp')
        sig = MakeSignal('TwoChirp',N);
        sig = sig(:);
        
    case lower('QuadChirp')
        sig = MakeSignal('QuadChirp',N);
        sig = sig(:);
        
    case lower('MishMash')
        sig = MakeSignal('MishMash',N);
        sig = sig(:);
        
    case lower('LinChirps')
        sig = MakeSignal('LinChirps',N);
        sig = sig(:); 
        
    case lower('Chirps')
        sig = MakeSignal('Chirps',N);
        sig = sig(:);       
     
    case lower('HypChirps')
        % strcmp(Name,'HypChirps'), % Hyperbolic Chirps of Mallat's book
        N = round(N*1.75);
        sig = MakeSignal('HypChirps',N);
        sig = sig(N/10+1:round(.65*N)+N/10);      
        sig = sig(1:2^(floor(log2(length(sig)))));
        sig = sig(:);

    case lower('chirp-MATLAB')% % Chirp model 2
        % Fs = 200;
        % t = -1:1/Fs:1-1/Fs;
        % fo = Fs/3; f1 = 2*Fs/3;
        % X = chirp(t,fo,1,f1,'q',[],'convex')';
        % % X = chirp(t,fo,1,f1,'li')';
        % figure(501); spectrogram(X,hanning(256),128,256,Fs,'yaxis'); shg
        % len_stream = length(t);
        % K = 0;    

        
    %--------%
    % Sounds %
    %--------%
    case {'chirp','gong','handel'}
        eval(sprintf('load %s.mat',sType));
        sig = y; % interp1(1:length(y),y,1:1/2:length(y))';  
        x = sig(:);
    case {'greasy','tweet'}
        eval('sig = ReadSignal(sType);');
        x = sig; 
    %--------%
    % Images %
    %--------%
    case {'cameraman','barbara','lena','peppers','airplane','baboon','sailboat','tiffany','boats','shapes','pirate','house'}
        n = sqrt(N);
        switch sType
            case {'shapes','boats'}
                load(sType);
            case {'pirate','cameraman'}
                sig = double(imread(sType,'tif'));
            otherwise
                sig = double(imread(sType,'ras'));
        end
        sig = imresize(sig,n/size(sig,1));
        [h0 h1 g0 g1] = daub79; sym = 1; wType = 'daub79';
        % [h0 h1 g0 g1] = dauborth(4); sym = 0; wType = 'daub4';
        % [h0 h1 g0 g1] = dauborth(8); sym = 0; wType = 'daub8';
        
        fprintf('Image: %s, size: %dx%d, wavelet: %s \n',sType, n,n, wType);
        J = floor(log2(n)-4);
        x = fwt2(sig,h0,h1,J,sym);
        x = x(:);
        sig = x;
        N = length(x);
        
        wave_struct = {};
        wave_struct.wType = wType;
        wave_struct.J = J;
        wave_struct.FilterBank = [h0; h1; g0; g1];
        wave_struct.sym = sym;
        
        take_fwt = 0;
         
    case 'phantom'
        n = sqrt(N);
        sig = phantom('Modified Shepp-Logan',n);
        [h0 h1 g0 g1] = dauborth(2); sym = 0; wType = 'haar';
        
        fprintf('Image: %s, size: %dx%d, wavelet: %s \n',sType, n,n, wType);
        J = floor(log2(n)-3);
        x = fwt2(sig,h0,h1,J,sym);
        x = x(:);
        N = length(x);
        
        wave_struct = {};
        wave_struct.wType = wType;
        wave_struct.J = J;
        wave_struct.FilterBank = [h0; h1; g0; g1];
        wave_struct.sym = sym;
         
        take_fwt = 0;
    otherwise
        disp('NAO');
        
end

if take_fwt == 1
    if isfield(in,'wType')
        wType = in.wType;
    else
        wType = 'daub4';
    end
    if exist('wType','var')
        switch lower(wType)
            case 'haar'
                [h0 h1 g0 g1] = dauborth(2); sym = 0; wType = 'haar';
            case 'daub4'
                [h0 h1 g0 g1] = dauborth(4); sym = 0; wType = 'daub4';
            case 'daub8'
                [h0 h1 g0 g1] = dauborth(8); sym = 0; wType = 'daub8';
            case 'daub79'
                [h0 h1 g0 g1] = daub79; sym = 1; wType = 'daub79';
            case 'daub1018'
                [h0 h1 g0 g1] = daub1018; sym = 2; wType = 'daub1018';
            otherwise
                disp('NAO');
        end
    else
        [h0 h1 g0 g1] = dauborth(4); sym = 0; wType = 'daub4';
    end
    if isfield(in,'J');
        J = in.J;
    else 
        J = floor(log2(N));
    end
    x = fwt(sig,h0,h1,J,sym);
    W_h = @(x) fwt(x,h0,h1,J,sym);
    iW_h = @(x) ifwt(x,g0,g1,J,sym);
    
    wave_struct = {};
    wave_struct.W_h = W_h;
    wave_struct.iW_h = iW_h;
else
    x = sig;
end

if isfield(in,'plot')
    if in.plot
        figure(1);
        subplot(211); set(gca,'FontSize',16);
        plot(sig,'LineWidth',2); set(gca,'XLim',[1 N]);
        axis tight; YLim = get(gca,'YLim');
        dY = YLim(2)-YLim(1);
        YLim = [YLim(1)-dY/10 YLim(2)+dY/10];
        set(gca,'YLim',YLim);
        set(gca,'XtickLabel','')
        ylabel(sType);
        subplot(212); set(gca,'FontSize',16);
        plot(x,'LineWidth',2); set(gca,'XLim',[1 N]);
        axis tight; YLim = get(gca,'YLim');
        dY = YLim(2)-YLim(1);
        YLim = [YLim(1)-dY/10 YLim(2)+dY/10];
        set(gca,'YLim',YLim);
        set(gca,'XtickLabel','')
        ylabel(wType);
    end
end

nout = max(nargout,1) - 1;
if nout > 0
    varargout{1} = sig;
end
if nout > 1
    varargout{2} = wave_struct;
end

