function rgbLinCalFormat = RGB2rgbLin(RGBCalFormat, Disp, imgParams)
% RGB2rgbLin  Convert RGB to linear rgb
%   rgbLin = RGB2rgbLin(RGBImage, Disp, imgParams)
%
%   Inputs:
%     RGBCalFormat      - gamma corrected RGB image
%     Disp              - Display structure
%   Outputs:
%     rgbLinCalFormat   - linear RGB image (MxNx3), range [0,1]

% Image format for reverse gamma correction
RGBImage = CalFormatToImage(RGBCalFormat,imgParams.m,imgParams.n);

% Reverse the gamma correction
iGtable      = displayGet(Disp.d,'inversegamma');
rgbLinImage  = rgb2dac(RGBImage,iGtable)/(2^displayGet(Disp.d,'dacsize')-1);

% Cal format for output
rgbLinCalFormat = ImageToCalFormat(rgbLinImage);

end
