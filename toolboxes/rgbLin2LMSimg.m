function lmsImage = rgbLin2LMSimg(rgbImage,T_cones,P_monitor,scaleFactor,m,n)

% function takes rgb image in image format and outputs LMS image in image format 
% 
% 

% Undo the scaling 
rgbImage = rgbImage * scaleFactor;

% Cal format
rgbImageCalFormat = ImageToCalFormat(rgbImage);

% LMS image
M_rgb2cones = T_cones*P_monitor;
lmsImageCalFormat = M_rgb2cones * rgbImageCalFormat;

% Turn back into image format from cal format
lmsImage = CalFormatToImage(lmsImageCalFormat,m,n);
end