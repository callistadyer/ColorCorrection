function RGBCalFormat = rgbLin2RGB(rgbLinCalFormat, Disp, imgParams)
% rgbLin2RGB  Convert linear rgb image to gamma corrected rgb
%   RGBImage = rgbLin2RGB(rgbLinImage, Disp)
%
%   Inputs:
%     rgbLinImage - linear RGB image (MxNx3), range [0,1]
%     Disp        - Struct containing:
%         Disp.d - display object with inverse gamma info
%
%   Outputs:
%     RGBImage    - gamma-corrected RGB image in normalized 0-1 range (MxNx3)

% Image format for gamma correction
rgbLinImage = CalFormatToImage(rgbLinCalFormat, imgParams.m,imgParams.n);

% Gamma correction
gammaTable = displayGet(Disp.d,'gammatable');
RGBImage = dac2rgb(rgbLinImage, gammaTable)*(2^displayGet(Disp.d,'dacsize')-1);

% Cal format for output
RGBCalFormat = ImageToCalFormat(RGBImage);

end
