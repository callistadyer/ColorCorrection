function rgbLinCalFormat = RGB2rgbLin(RGBCalFormat, Disp, imgParams)
% RGB2rgbLin  Convert RGB to linear rgb
%   rgbLin = RGB2rgbLin(RGBImage, Disp, imgParams)
%
%   Inputs:
%     RGBCalFormat - RGB image in normalized 0-1 range 
%     Disp         
%   Outputs:
%     rgbLin   - linear RGB image (MxNx3), range [0,1]

% Image format for reverse gamma correction
RGBImage = CalFormatToImage(RGBCalFormat,imgParams.m,imgParams.n);

% Reverse the gamma correction
iGtable      = displayGet(Disp.d,'inversegamma');
rgbLinImage  = rgb2dac(RGBImage,iGtable)/(2^displayGet(Disp.d,'dacsize')-1);

% Cal format for output
rgbLinCalFormat = ImageToCalFormat(rgbLinImage);

end
