function rgbLinCalFormat = RGB2rgbLin(RGBCalFormat, Disp, imgParams)
% RGB2rgbLin  Convert RGB to linear rgb
%   rgbLinCalFormat = RGB2rgbLin(RGBCalFormat, Disp, imgParams)
%
%   Inputs:
%     RGBCalFormat      - gamma corrected RGB image
%     Disp              - Display structure
%   Outputs:
%     rgbLinCalFormat   - linear RGB image (MxNx3), range [0,1]

% Image format for reverse gamma correction
RGBImage = CalFormatToImage(RGBCalFormat,imgParams.m,imgParams.n);
RGB_norm = double(RGBImage) / 255;

% dac is gamma corrected RGB values
% Use gamma table to do inverse gamma correction
% Use inverse gamma table to do gamma correction
% Reverse the gamma correction
gammaTable      = displayGet(Disp.d,'gamma table');
rgbLinImage = dac2rgb(RGB_norm, gammaTable);  % leave as [0,1] linear
% rgbLinImage = dac2rgb(RGBImage, gammaTable);
% rgbLinImage     = dac2rgb(RGB_norm,gammaTable)/(2^displayGet(Disp.d,'dacsize')-1);

% Cal format for output
rgbLinCalFormat = ImageToCalFormat(rgbLinImage);

end
