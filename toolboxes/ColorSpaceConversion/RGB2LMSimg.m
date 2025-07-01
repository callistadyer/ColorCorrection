function lmsImage = RGB2LMSimg(RGBImage,Disp,scaleFactor,bScale)
% Function takes in RGB image in image format, and outputs LMS image in image format 
% also calls function rgb2LMSimg.m
%
% This expects a gamma corrected (ready to display) RGB image, and undoes
% the gamma correction (aka linearizes) the image before applying the
% transformation to LMS.
%
% Syntax:
%   lmsImage = RGB2LMSimg(RGBImage,d,T_cones,P_monitor,scaleFactor,m,n,bScale)
%
% Description:
%
% Inputs:
%   RGBImage              - RGB image, gamma corrected.
%   d                     - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%   T_cones               - [3xnWl]. Cone spectral sensitivities
%   P_monitor             - [nWlx3]. Display primaries
%   scaleFactor           - 
%   m                     - Scalar.  Row dimension of image     
%   n                     - Scalar.  Column dimension of image     
%   bScale                - Boolean. Scale or not 
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

% Undo scaling and convert to LMS
lmsImage = rgbLin2LMSimg(rgbLinImage,Disp.T_cones,Disp.P_monitor,scaleFactor,Disp.m,Disp.n,bScale);

end