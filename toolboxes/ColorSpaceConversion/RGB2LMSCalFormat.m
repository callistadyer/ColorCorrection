function LMSCalFormat = RGB2LMSCalFormat(RGBImage,Disp)
% Function takes in RGB image in image format, and outputs LMS image in image format 
% also calls function rgb2LMSimg.m
%
% This expects a gamma corrected (ready to display) RGB image, and undoes
% the gamma correction (aka linearizes) the image before applying the
% transformation to LMS.
%
% Syntax:
%   LMSCalFormat = RGB2LMSCalFormat(RGBImage,Disp)
%
% Description:
%
% Inputs:
%   RGBImage              - RGB image, gamma corrected.
%   Disp contains
%       Disp.d                     - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%       Disp.T_cones               - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor             - [nWlx3]. Display primaries
%       Disp.m                     - Scalar.  Row dimension of image     
%       Disp.n                     - Scalar.  Column dimension of image  
%
% Outputs:
%   lmsImage              - LMS image in image format
%
% Optional key/value pairs:
%   None
%


% Reverse the gamma correction
gammaTable = displayGet(Disp.d,'gammatable');
rgbLinImage = dac2rgb(RGBImage, gammaTable)*(2^displayGet(Disp.d,'dacsize')-1);

rgbLinCalFormat = ImageToCalFormat(rgbLinImage,Disp.m,Disp.n);

% Convert to LMS
LMSCalFormat = rgbLin2LMSCalFormat(rgbLinCalFormat,Disp);

end