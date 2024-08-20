function [rgbImageCalFormat,scaleFactor] = LMS2rgbLinimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n,bScale)

% function takes an LMS image in cal format and outputs and rgb image in cal format
% also outputs the scale factor used to normalize the rgbImage 
%

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
rgbImageCalFormat = M_cones2rgb*lmsImageCalFormat;
rgbImage = CalFormatToImage(rgbImageCalFormat,m,n);

if bScale == 1
    % For right now, normalize so that maximum value in rgb is 1
    scaleFactor = max(rgbImage(:)); % save scale factor for later
    rgbImage = rgbImage/scaleFactor;
else
    scaleFactor = 1;
end

% Transform to cal format
rgbImageCalFormat = ImageToCalFormat(rgbImage);

rgbImageTruncate = rgbImageCalFormat;
rgbImageTruncate(rgbImageTruncate < 0) = 0;

rgbImageCalFormat = rgbImageTruncate;
end