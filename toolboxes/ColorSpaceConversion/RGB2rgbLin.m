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

% Get gamma table from display (should be normalized to [0,1] already)
gammaTable = displayGet(Disp.d, 'gamma table');
nLevels = size(gammaTable, 1);
xvals = linspace(0,1,nLevels); 

% Compute indices for mapping
indices = round(RGBCalFormat * (nLevels - 1)) + 1;

% Make sure indices are all good and no numerical junk...
indices(indices < 1)       = 1;
indices(indices > nLevels) = nLevels;

% Initialize
rgbLinCalFormat = zeros(size(RGBCalFormat));

% Loop over RGB
for ch = 1:3
    rgbLinCalFormat(ch, :) = gammaTable(indices(ch, :), ch)';
end

%% Alternatively, use the dac2rgb. Less sure about how this works.
% Reverse the gamma correction
% gammaTable      = displayGet(Disp.d,'gamma table');
% rgbLinImage     = dac2rgb(RGBImage,gammaTable); %/(2^displayGet(Disp.d,'dacsize')-1);

% Cal format for output
% rgbLinCalFormat = ImageToCalFormat(rgbLinImage);

end
