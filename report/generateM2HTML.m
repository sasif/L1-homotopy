clear

% Setup path
m2html_path = 'C:\Users\Asif\Dropbox\Homotopy toolbox\m2html';
addpath(m2html_path);

mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

cd ..  
cd ..
m2html('mfiles','L1_homotopy_v2.0', 'htmldir','L1_homotopy_v2.0\report\doc','recursive','on', 'global','on','graph','on','template','frame','index','menu');