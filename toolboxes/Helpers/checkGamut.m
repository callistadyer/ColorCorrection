function [inGamut, rgbImageCalFormat] = checkGamut(LMSimageCalFormat,Disp)

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
[rgbImageCalFormat] = LMS2rgbLinCalFormat(LMSimageCalFormat,Disp);
rgbImage = CalFormatToImage(rgbImageCalFormat,Disp.m,Disp.n);

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
