function rgbLinCalFormat = RGB2rgbLin(RGBImage, Disp,imgParams)
% RGB2rgbLin  Convert sRGB DAC image to linear RGB using display gamma table.
%   rgbLin = RGB2rgbLin(RGBImage, Disp)
%
%   Inputs:
%     RGBImage - sRGB image in normalized 0-1 range (MxNx3)
%     Disp     - Struct containing:
%         Disp.d - display object with gamma info
%
%   Outputs:
%     rgbLin   - linear RGB image (MxNx3), range [0,1]

% Get display parameters
% Reverse the gamma correction
gammaTable = displayGet(Disp.d,'gammatable');
rgbLinImage = dac2rgb(RGBImage, gammaTable)*(2^displayGet(Disp.d,'dacsize')-1);

rgbLinCalFormat = ImageToCalFormat(rgbLinImage);

end
