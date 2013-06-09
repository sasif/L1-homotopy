function Compile_FWT

MEX_OK = 0;

% for file = {'fwt','ifwt','fwt2','ifwt2', 'fwt2_CWT', 'afwt2_CWT',  ...
%         'ifwt2_CWT', 'afwt','afwt2'}
%     
%     file = char(file);
%     if exist(file)~=3,
%         MEX_OK = 0;
%         break;
%     end
% end
disp('Compiling Wavelet mex files');

if ~MEX_OK
    Friend = computer;
    isPC = 0;
    if strcmp(Friend(1:2),'PC')
        isPC = 1;
    end
    if isPC
        mex -c fwt_level.c
        mex fwt.c fwt_level.obj -output fwt
        mex ifwt.c fwt_level.obj -output ifwt
        mex fwt2.c fwt_level.obj -output fwt2
        mex ifwt2.c fwt_level.obj -output ifwt2
        
    else
        mex -c fwt_level.c
        mex fwt.c fwt_level.o -output fwt
        mex ifwt.c fwt_level.o -output ifwt
        mex fwt2.c fwt_level.o -output fwt2
        mex ifwt2.c fwt_level.o -output ifwt2
    end
end
