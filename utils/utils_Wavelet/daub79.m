% daub79.m
%
% Returns the filter coefficients (for lowpass and highpass) for the
% Daubechies 7,9 biorthogonal wavelet set
% Usage : [h0,h1,g0,g1] = daub79
% h0 - lowpass analysis
% h1 - highpass analysis
% g0 - lowpass synthesis
% g1 - highpass synthesis

function [h0, h1, g0, g1] = daub79()

b = sqrt(2)*[0.6029490182363579 0.2668641184428723 -0.07822326652898785 ...
  -0.01686411844287495 0.02674875741080976];
c = sqrt(2)/2*[1.115087052456994 -0.5912717631142470 -0.05754352622849957 ...
      0.09127176311424948]; 


h0 = [0 fliplr(b) b(2:5)];
h1 = [0 -fliplr(c) -c(2:4) 0 0];
g0 = ((-1).^(1:10)).*h1;
g1 = ((-1).^(0:9)).*h0;

% NOTE: This function is equivalent to the following using 
%  the MATLAB wavelet toolbox:
% [rf,df] = biorwavf('bior4.4');
% [h0,h1,g0,g1] = biorfilt(df,rf);

