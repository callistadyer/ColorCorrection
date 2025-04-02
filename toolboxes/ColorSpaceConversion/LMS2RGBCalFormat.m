function [RGBCalFormat,rgbLinCalFormat] = LMS2RGBCalFormat(lmsImageCalFormat,Disp)

% Converts LMS cone images to RGB images. Outputs both linear and gamma
% corrected rgb/RGB values
%
% Syntax:
%   [RGBImageCalFormat,rgbLinImageCalFormat] = LMS2RGBCalFormat(lmsImageCalFormat,Disp)
%
% Description:
%
% Inputs:
%   lmsImageCalFormat     - [3xnPix] matrix. Cal formatted LMS image
%   Disp contains:
%       Disp.d                     - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%       Disp.T_cones               - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor             - [nWlx3]. Display primaries
%       Disp.m                     - Scalar.  Row dimension of image     
%       Disp.n                     - Scalar.  Column dimension of image     
%
% Outputs:
%   RGBImageCalFormat     - RGB image in cal format
%   rgbLinImageCalFormat  - linear rgb image in cal format
%
% Optional key/value pairs:
%   None
%

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
[rgbImageCalFormat] = LMS2rgbLinCalFormat(lmsImageCalFormat,Disp);
rgbImage = CalFormatToImage(rgbImageCalFormat,Disp.m,Disp.n);

% Linear rgb values make sure not below 0
rgbLinImage = rgbImage; 

rgbLinCalFormat = ImageToCalFormat(rgbLinImage);

% Gamma correct
iGtable = displayGet(Disp.d,'inversegamma');
RGBImage = rgb2dac(rgbLinImage,iGtable)/(2^displayGet(Disp.d,'dacsize')-1);

% Transform to cal format
RGBCalFormat = ImageToCalFormat(RGBImage);

end