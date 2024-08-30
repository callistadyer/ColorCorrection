function [RGBImageCalFormat,rgbLinImageCalFormat,scaleFactor] = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n,bScale)

% Converts LMS cone images to RGB images. Outputs both linear and gamma
% corrected rgb/RGB values
%
% Syntax:
%   [RGBImageCalFormat,rgbImageCalFormat,scaleFactor] = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n,bScale)
%
% Description:
%
% Inputs:
%   lmsImageCalFormat     - [3xnPix] matrix. Cal formatted LMS image
%   d                     - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%   T_cones               - [3xnWl]. Cone spectral sensitivities
%   P_monitor             - [nWlx3]. Display primaries
%   m                     - Scalar.  Row dimension of image     
%   n                     - Scalar.  Column dimension of image     
%   bScale                - Boolean. Scale or not 
%
% Outputs:
%   None
%
% Optional key/value pairs:
%   None
%

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
% THIS SHOULD CALL LMS2rgbLinimg().
rgbLinImageCalFormat = M_cones2rgb*lmsImageCalFormat;
rgbImage = CalFormatToImage(rgbLinImageCalFormat,m,n);

if bScale == 1
    % For right now, normalize so that maximum value in rgb is 1
    scaleFactor = max(rgbImage(:)); % save scale factor for later
    rgbImage = rgbImage/scaleFactor;
else
    scaleFactor = 1;
end

% Truncated version for gamma correction
rgbImageTruncate = rgbImage;
rgbImageTruncate(rgbImageTruncate < 0) = 0;

% Linear rgb values make sure not below 0
rgbLinImageCalFormat = rgbImageTruncate; 

% Gamma correct
iGtable = displayGet(d,'inversegamma');
RGBImage = rgb2dac(rgbImageTruncate,iGtable)/(2^displayGet(d,'dacsize')-1);

% Transform to cal format
RGBImageCalFormat = ImageToCalFormat(RGBImage);

end