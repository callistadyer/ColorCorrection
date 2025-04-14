function [rgbImageCalFormat] = LMS2rgbLinCalFormat(lmsImageCalFormat,Disp)

% function takes an LMS image in cal format and outputs and rgb image in cal format
% also outputs the scale factor used to normalize the rgbImage 
%
% Syntax:
%   [rgbImageCalFormat] = LMS2rgbLinimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n)
%
% Description:
%
% Inputs:
%   lmsImageCalFormat     - [3xnPix] matrix. Cal formatted LMS image
%   Disp contains
%       Disp.d                     - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%       Disp.T_cones               - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor             - [nWlx3]. Display primaries
%       Disp.m                     - Scalar.  Row dimension of image     
%       Disp.n                     - Scalar.  Column dimension of image     
%
% Outputs:
%   rgbLinImageCalFormat  - linear rgb image in cal format
%
% Optional key/value pairs:
%   None
%

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
rgbImageCalFormat = M_cones2rgb*lmsImageCalFormat;

end