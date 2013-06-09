% Creates a DWT representation matrix without signal extension

function Psi = create_DWT(in)

% A wavelet basis function using orthogonal wavelets
%
% inputs: 
%   J -- finest scale for wavelets
%   wType -- type of wavelets 
%   sym -- type of extension
%   N -- length of signal interval

N = in.N;  
wType = in.wType; 
J = in.J;
sym = in.sym;

if strcmpi(wType(1:4),'daub')
    switch wType
        case 'daub79'
            [h0 h1 g0 g1] = daub79; sym = 1;
        case 'daub1018'
            [h0 h1 g0 g1] = daub1018; sym = 2;
            %         case 'daub4'
            %             filter_length = 4;
            %             [h0 h1 g0 g1] = daub4;
            %         case 'daub8'
            %             filter_length = 8;
            %             [h0 h1 g0 g1] = daub8;
        otherwise
            filter_length = str2num(wType(5:end));
            [h0 h1 g0 g1] = dauborth(filter_length);
    end
else
    error('Not implemented yet');
end

% switch lower(wType)
%     case 'daub2'
%         [h0 h1 g0 g1] = dauborth(2);
%         filter_length = 2;
%     case 'daub4'
%         [h0 h1 g0 g1] = dauborth(4);
%         filter_length = 4;
%     case 'daub8'
%         [h0 h1 g0 g1] = dauborth(8);
%         filter_length = 8;
%     case 'daub12'
%         [h0 h1 g0 g1] = dauborth(12);
%         filter_length = 12;
%     otherwise
%         disp('NAO');
% end

% compute the lengths of functions at every scale
if sym == 3
    w_len = filter_length;
    % for j = 1:J-1; w_len(end+1) = w_len(end)+2^j*(filter_length-1); end
    for j = 1:J-1; w_len(end+1) = 2*w_len(end)-1+(filter_length-1); end
    w_len(end+1) = w_len(end); 
else
    w_len = N;
end
% 
if w_len(end) > N && sym == 3;
    error('reduce the scale J or increase length N - length of scaling function is larger than N');
end

% compute synthesis matrix... 
switch sym
    case 3
        iW_h = @(x) idwtmult1_conv(x,g0,g1,J);
        % iW_h = @(x) dwtmult1_conv(x,h0,h1,J);
    otherwise
        iW_h = @(x) ifwt(x,g0,g1,J,sym);
end
Psi = [];
for ii = 1:N;   
    Psi(:,end+1) = iW_h(circshift([1; zeros(N-1,1)],ii-1));    
end

% %-----------------% %
% % debug 
% %-----------------% %
% PSI = [[Psi; zeros(N)] [zeros(N); Psi]];
% vec = @(z) z(:); 
% max(abs(vec(PSI'*PSI-eye(2*N))))
% %%
% in = []; in.J = 4; in.N = 256; in.wType = 'daub10'; in.sym = 3;
% Psi = create_DWT(in);