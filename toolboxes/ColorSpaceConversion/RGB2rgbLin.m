function rgbLin = RGB2rgbLin(RGBImage, Disp)
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
gammaTable = displayGet(Disp.d, 'gammatable');
dacSize    = displayGet(Disp.d, 'dacsize');

% Convert normalized 0-1 → DAC
dacValues = RGBImage * (2^dacSize - 1);

% Convert DAC → linear RGB using gamma table
rgbLin = dac2rgb(dacValues, gammaTable);
end
