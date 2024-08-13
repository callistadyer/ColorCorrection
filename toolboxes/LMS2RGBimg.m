function [RGBImageCalFormat,scaleFactor] = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n)

% function takes an LMS image in cal format and outputs and RGB image in cal format
% also outputs the scale factor used to normalize the rgbImage 
%

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
rgbImageCalFormat = M_cones2rgb*lmsImageCalFormat;
rgbImage = CalFormatToImage(rgbImageCalFormat,m,n);

% For right now, normalize so that maximum value in rgb is 1
scaleFactor = max(rgbImage(:)); % save scale factor for later 
rgbImage = rgbImage/scaleFactor;

% Truncated version for gamma correction
rgbImageTruncate = rgbImage;
rgbImageTruncate(rgbImageTruncate < 0) = 0;

% Gamma correct
iGtable = displayGet(d,'inversegamma');
RGBImage = rgb2dac(rgbImageTruncate,iGtable)/(2^displayGet(d,'dacsize')-1);

% Transform to cal format
RGBImageCalFormat = ImageToCalFormat(RGBImage);

end