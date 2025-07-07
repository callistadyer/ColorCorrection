function [RGBCalFormat,rgbLinCalFormat] = LMS2RGBCalFormat(lmsCalFormat,Disp,imgParams)

% Converts LMS cone images to RGB images. Outputs both linear and gamma
% corrected rgb/RGB values
%
% Syntax:
%   [RGBCalFormat,rgbLinCalFormat] = LMS2RGBCalFormat(lmsCalFormat,Disp,imgParams)
%
% Description:
%
% Inputs:
%   lmsCalFormat     - [3xnPix] matrix. Cal formatted LMS image
%   Disp contains:
%       Disp.d                     - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%       Disp.T_cones               - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor             - [nWlx3]. Display primaries
%   imgParams contains:
%       imgParams.m
%       imgParams.n
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
[rgbLinCalFormat] = LMS2rgbLinCalFormat(lmsCalFormat,Disp);

% Gamma correct
% iGtable  = displayGet(Disp.d,'inversegamma');
% RGBImage = rgb2dac(rgbLinImage,iGtable)/(2^displayGet(Disp.d,'dacsize')-1);
RGBCalFormat = rgbLin2RGB(rgbLinCalFormat, Disp);

end