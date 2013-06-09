% Run this file to set up MATLAB paths to all the sub-directories.
%
% Alternatively you can directly add L1_homotopy folder along with all the 
% subfolders, using MATLAB menu as
% (File --> Set Path --> Add with subfolders -->...)

% file_name = mfilename('fullpath');
% dbs = dbstack;
% 
% addpath(genpath(file_name(1:length(file_name)-length(dbs.name))));

mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath(genpath(mdir))
