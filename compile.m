%% May need to compile the following mex files.
% mvprod (matrix-vector product on restricted support)
% realnoiselet (noiselet transform for random sampling of images)
% Wavelet transform (fwt, fwt2, ifwt, ifwt2, ...)

mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

cd utils
disp('Compiling mvprod.c');
mex mvprod.c

cd utils_Wavelet
Compile_FWT
cd ..

cd utils_meas
disp('Compiling realnoiselet.c');
mex realnoiselet.c

disp('Done compiling!');

cd(mdir);