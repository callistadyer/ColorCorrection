function LMSCalFormat = RGB2LMSCalFormat(RGBCalFormat,Disp)
% Function takes in gamma corrected RGB in cal format, and outputs LMS image in cal format 
%
% This expects a gamma corrected (ready to display) RGB image, and undoes
% the gamma correction (aka linearizes) the image before applying the
% transformation to LMS.
%
% Syntax:
%   LMSCalFormat = RGB2LMSCalFormat(RGBCalFormat,Disp)
%
% Description:
%
% Inputs:
%   RGBImage              - RGB image, gamma corrected.
%   Disp contains
%       Disp.d                     - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%       Disp.T_cones               - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor             - [nWlx3]. Display primaries
%
% Outputs:
%   lmsImage              - LMS image in image format
%
% Optional key/value pairs:
%   None
%


% Reverse the gamma correction
% gammaTable = displayGet(Disp.d,'gammatable');
% rgbLinImage = dac2rgb(RGBImage, gammaTable)*(2^displayGet(Disp.d,'dacsize')-1);
rgbLinCalFormat = RGB2rgbLin(RGBCalFormat, Disp);

% Convert to LMS
LMSCalFormat = rgbLin2LMSCalFormat(rgbLinCalFormat,Disp);

end