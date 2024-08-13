function lmsImage = RGB2LMSimg(RGBImage,d,T_cones,P_monitor,scaleFactor,m,n)

% function takes in RGB image in image format, and outputs LMS image in image format 
% also calls function rgb2LMSimg.m
% 

% Reverse the gamma correction
gammaTable = displayGet(d,'gammatable');
rgbImage = dac2rgb(RGBImage, gammaTable)*(2^displayGet(d,'dacsize')-1);

% Undo scaling and convert to LMS
lmsImage = rgbLin2LMSimg(rgbImage,T_cones,P_monitor,scaleFactor,m,n);

end