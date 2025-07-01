function RGBImage = rgbLin2RGB(rgbLinImage, Disp)
% rgbLin2RGB  Convert linear RGB image to sRGB DAC image using display inverse gamma.
%   RGBImage = rgbLin2RGB(rgbLinImage, Disp)
%
%   Inputs:
%     rgbLinImage - linear RGB image (MxNx3), range [0,1]
%     Disp        - Struct containing:
%         Disp.d - display object with inverse gamma info
%
%   Outputs:
%     RGBImage    - gamma-corrected sRGB image in normalized 0-1 range (MxNx3)


iGtable  = displayGet(Disp.d,'inversegamma');
RGBImage = rgb2dac(rgbLinImage,iGtable)/(2^displayGet(Disp.d,'dacsize')-1);

end
