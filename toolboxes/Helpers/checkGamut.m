function [inGamut, rgbImageCalFormat] = checkGamut(LMSimageCalFormat,Disp,bScale)

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
[rgbImageCalFormat,scaleFactor] = LMS2rgbLinCalFormat(LMSimageCalFormat,Disp,bScale);
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

% Values of rgbImageTruncate should be between 0 and 1... if not, there's
% gonna be an error in rgb2dac
if any(rgbImageTruncate(:) > 1 | rgbImageTruncate(:) < 0)
    inGamut = 0;
    % min(rgbImage(:))
    % max(rgbImage(:))
else 
    inGamut = 1;
end
