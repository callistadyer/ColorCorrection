function [RGBImageCalFormat,rgbLinImageCalFormat,scaleFactor] = LMS2RGBCalFormat(lmsImageCalFormat,Disp,bScale)

% Converts LMS cone images to RGB images. Outputs both linear and gamma
% corrected rgb/RGB values
%
% Syntax:
%   [RGBImageCalFormat,rgbLinImageCalFormat,scaleFactor] = LMS2RGBCalFormat(lmsImageCalFormat,Disp,bScale)
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
%   bScale                - Boolean. Scale or not 
%
% Outputs:
%   RGBImageCalFormat     - RGB image in cal format
%   rgbLinImageCalFormat  - linear rgb image in cal format
%   scaleFactor           - scale factor used to scale img values
%
% Optional key/value pairs:
%   None
%

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
[rgbImageCalFormat,scaleFactor] = LMS2rgbLinCalFormat(lmsImageCalFormat,Disp,bScale);
rgbImage = CalFormatToImage(rgbImageCalFormat,Disp.m,Disp.n);

if bScale == 1
    % For right now, normalize so that maximum value in rgb is 1
    scaleFactor = max(rgbImage(:)); % save scale factor for later
    rgbImage = rgbImage/scaleFactor;
else
    scaleFactor = 1;
end

% Truncated version for gamma correction
rgbImageTruncate = rgbImage;
% rgbImageTruncate(rgbImageTruncate < 0) = 0;

% Linear rgb values make sure not below 0
rgbLinImage = rgbImageTruncate; 
rgbLinImageCalFormat = ImageToCalFormat(rgbLinImage);
% rgbImageTruncate(rgbImageTruncate > 1) = 1;

% Values of rgbImageTruncate should be between 0 and 1... if not, there's
% gonna be an error in rgb2dac
if any(rgbImageTruncate(:) > 1 | rgbImageTruncate(:) < 0)
    error(['LMS2RGBCalFormat: WARNING! rgb values are out of gamut... rgbImageTruncate values outside of the range [0 1]']);
end

% Gamma correct
iGtable = displayGet(Disp.d,'inversegamma');
RGBImage = rgb2dac(rgbImageTruncate,iGtable)/(2^displayGet(Disp.d,'dacsize')-1);

% Transform to cal format
RGBImageCalFormat = ImageToCalFormat(RGBImage);

end