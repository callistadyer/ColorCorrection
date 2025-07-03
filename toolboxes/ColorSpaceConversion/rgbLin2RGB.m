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
%
% I’m converting a linear RGB image into a gamma-corrected 
% RGB image by mapping each linear RGB value to its corresponding gamma-corrected 
% output using an inverse gamma table. The inverse gamma table has 3 columns 
% (for R, G, B) and nLevels rows (e.g., 1024 levels), which represent the 
% possible DAC output values for each channel. For each pixel’s R, G, and B 
% component—assumed to be within the [0, 1] range, I calculate an index into
% the inverse gamma table by scaling the input value to the number of table entries..
% index = linspace(0,1,1024) like that, rounding to the nearest integer, 
% and then clipping it to ensure it stays within  bounds (1 to nLevels).
% For each channel, I replace the linear input with the corresponding value
% from the inverse gamma table


% Image format for gamma correction
rgbLinImage = CalFormatToImage(rgbLinCalFormat, imgParams.m,imgParams.n);

% Gamma correction
invGammaTable = displayGet(Disp.d,'inverse gamma');
dacSize = displayGet(Disp.d, 'dacsize');
maxDAC  = 2^dacSize - 1;
invGammaTable = invGammaTable / maxDAC;

% Number of levels in the inverse gamma table
nLevels = size(invGammaTable, 1);
xvals = linspace(0,1,nLevels); 

% Compute lookup indices 
indices = round(rgbLinCalFormat * (nLevels - 1)) + 1;

% Make sure indices are all good and no numerical junk...
indices(indices < 1)       = 1;
indices(indices > nLevels) = nLevels;

% Initialize 
RGBCalFormat = zeros(size(rgbLinCalFormat));

% Loop over RGB 
for ch = 1:3
    RGBCalFormat(ch, :) = invGammaTable(indices(ch, :), ch)';
end


%% Alternatively, use rgb2dac
% RGBImage = rgb2dac(rgbLinImage, invgammaTable); 
% dacsize = displayGet(Disp.d, 'dacsize');
% RGBImage = (double(RGBImage) - 1) / (2^dacsize - 1);

% Cal format for output
% RGBCalFormat = ImageToCalFormat(RGBImage);

end
