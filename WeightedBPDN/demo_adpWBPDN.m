% demo_adpWBPDN.m
%
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \Sum \w_i |x_i| + 1/2*||y-Ax||_2^2
%
% while adaptively selecting the weights w_i
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@gatech.edu
% Created: June 16, 2011

clear
close all

%% parfor setup
% numCores = str2double(getenv('NUMBER_OF_PROCESSORS'));
% mpSize = numCores-1;
% if matlabpool('size') ~= mpSize
%     if matlabpool('size')~=0
%         matlabpool close;
%     else
%         matlabpool('open', mpSize);
%     end
% end

%% Setup path
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath ../Pursuits_Homotopy/
addpath ../utils/
addpath ../utils/utils_Wavelet
addpath ../utils/utils_meas
addpath ../solvers/
addpath src/

disp(['--------------------',datestr(now),'-------------------------'])

% simulation parameters
rwt_mode = 5;
lambda = 0;


SAVE_RESULTS = false;

% SAVE_RESULTS = true;  diary(sprintf('%s-reproduce.txt',mname));

largescale = 1;
SNR = 40;
R = 2; 

IMG_LIST = {'barbara','boats', 'cameraman','house','peppers','shapes','lena','airplane','baboon','sailboat','tiffany'};
sType = IMG_LIST{3};

if largescale
    mType = 'noiselets';
    % sType = 'cameraman'; N = (256)^2;
    % sType = 'boats'; N = (256)^2;
    % sType = 'peppers'; N = (256)^2;
    % sType = 'house'; N = (256)^2;
    % sType = 'shapes'; N = (256)^2;
    % sType = 'pirate'; N = (512)^2;
    N = (256)^2;
    M = round(N/R);
    %M = 30000;
    T = N;
    str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M) = %d, %d, rwt_mode-%d, lambda%3.4g.', mType, sType, SNR, N, M, rwt_mode, lambda);
    disp(str0);
else
    mType = 'randn';
    sType = 'HeaviSine';
    N = 512;   % signal length
    M = round(N/2);    % no. of measurements
    T = round(M/3);    % sparsity level
    str0 = sprintf('mType-%s, sType-%s, SNR = %d, (N,M,T) = %d, %d, %d, rwt_mode-%d, lambda%3.4g.', mType, sType, SNR, N, M, T, rwt_mode, lambda);
    disp(str0);
end

% rank-1 update mode
delx_mode = 'mil'; % mil or qr

%% Simulation
maxsim = 1;
script_simulation_adpWBPDN
