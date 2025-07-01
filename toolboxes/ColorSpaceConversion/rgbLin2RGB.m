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

% Get display parameters
inverseGtable = displayGet(Disp.d, 'inversegamma');
dacSize       = displayGet(Disp.d, 'dacsize');

% Convert linear RGB → DAC using inverse gamma table
dacValues = rgb2dac(rgbLinImage, inverseGtable);

% Normalize DAC to 0-1 sRGB
RGBImage = dacValues / (2^dacSize - 1);
end
