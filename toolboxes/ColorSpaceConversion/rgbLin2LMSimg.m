function lmsImage = rgbLin2LMSimg(rgbImage,Disp,scaleFactor,bScale)

% function takes rgb image in image format and outputs LMS image in image format 
%
% Syntax:
%   lmsImage = rgbLin2LMSimg(rgbImage,T_cones,P_monitor,scaleFactor,m,n,bScale)
%
% Description:
%
% Inputs:
%   rgbImage              - linear rgb image
%   Disp includes:
%       Disp.T_cones               - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor             - [nWlx3]. Display primaries
%       Disp.m                     - Scalar.  Row dimension of image     
%       Disp.n                     - Scalar.  Column dimension of image  
%   scaleFactor           - scale factor used to scale img values   
%   bScale                - Boolean. Scale or not 
%
% Outputs:
%   lmsImage              - LMS image
%
% Optional key/value pairs:
%   None
%

if bScale == 1
% Undo the scaling 
rgbImage = rgbImage * scaleFactor;
end

% Cal format
rgbImageCalFormat = ImageToCalFormat(rgbImage);

% LMS image
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
lmsImageCalFormat = M_rgb2cones * rgbImageCalFormat;

% Turn back into image format from cal format
lmsImage = CalFormatToImage(lmsImageCalFormat,Disp.m,Disp.n);
end