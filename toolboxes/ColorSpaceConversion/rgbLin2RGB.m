function RGBCalFormat = rgbLin2RGB(rgbLinCalFormat, Disp, imgParams)
% rgbLin2RGB  Convert linear rgb image to gamma corrected rgb
%   RGBCalFormat = rgbLin2RGB(rgbLinCalFormat, Disp, imgParams)
%
%   Inputs:
%     rgbLinCalFormat  - linear rgb in cal format 
%     Disp             - Struct containing:
%     imgParams
%
%   Outputs:
%     RGBCalFormat    - gamma-corrected RGB cal format 

% Image format for gamma correction
rgbLinImage = CalFormatToImage(rgbLinCalFormat, imgParams.m,imgParams.n);

% Gamma correction
invgammaTable = displayGet(Disp.d,'inverse gamma');
% RGBImage   = rgb2dac(rgbLinImage, invgammaTable)*(2^displayGet(Disp.d,'dacsize')-1);

RGBImage = rgb2dac(rgbLinImage, invgammaTable); 
% Put in range [0 1]
% RGBImage = double(RGBImage) / (2^displayGet(Disp.d,'dacsize') - 1);
RGBImage = double(RGBImage - 1) / (2^displayGet(Disp.d,'dacsize') - 1);

% Cal format for output
RGBCalFormat = ImageToCalFormat(RGBImage);

end
