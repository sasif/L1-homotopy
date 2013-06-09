% dwtmult1_conv.m
%
% Performs multiple levels of the discrete wavelet transform with linear 
% filtering (using dwtlevel1 with sym == 3)
% Usage : w = dwtmult1(x, h0, h1, L)
%
% Modified from dwtmult1.m by Justin Romberg

% NOT TESTED YET... 

function w_conv = dwtmult1_conv(x, h0, h1, L)

sym = 3; 

N = length(x);
w_conv = [];
wl = x(:)';
for ll = 1:L
  w = dwtlevel1(wl, h0, h1, sym);
  wl = w(1,:);
  wh = w(2,:);
  w_conv = [wh w_conv];
end
w_conv = [wl w_conv];
