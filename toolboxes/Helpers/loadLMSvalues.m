function [triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = loadLMSvalues(img,renderType,Disp,imgParams,varargin)
% loadLMSvalues  Loads or generates an image and converts to LMS for dichromat simulation
%
% Syntax:
%   [triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = loadLMSvalues(img,renderType,Disp,imgParams,varargin)
%
% Inputs:
%   img:              Either 'ishihara' or a filename ('.png', '.jpg') or a hyperspectral identifier
%   renderType:       Type of dichromacy to simulate
%                         'Deuteranopia', 'Protanopia', or 'Tritanopia'
%   imgParams
%   Disp:             Structure containing display calibration, cone sensitivities, and image dimensions
%
% Outputs:
%   triLMSCalFormat:  LMS image in calibration format for trichromatic viewer
%   triRGBCalFormat:  RGB image 
%   Disp:             Updated display struct with possibly new dimensions
%
% Description:
%   Depending on the input image, this function performs one of several operations:
%   (1) For 'ishihara', it generates a synthetic Ishihara plate with embedded digit
%   (2) For file inputs (.png, .jpg), it loads and converts the image to LMS
%   (3) For hyperspectral inputs, it calls a renderer to generate simulated LMS scenes
%   In all cases, it produces LMS versions for both trichromat and dichromat viewers.
%
%
% History:
%   04/14/2025  cmd  Initial draft and structure
%
% Examples:
%{
Disp = loadDisplay();
imgParams = buildSetParameters('ishihara',1,128,128);
[triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = loadLMSvalues('ishihara','Deuteranopia',Disp,imgParams,'clearTestImages', false);

% Check that behavior has not changed since we declared it good.
% if (abs(sum(testLMS(:)) - 525.8024)/525.8024 > 1e-4)
%     error('No longer get same LMS image returned by loadLMSValues');
% end
%}

% key-value pairs
p = inputParser;
defaultClearTestImages = false;              % Set default value for 'clearTestImages' key
addParameter(p, 'clearTestImages', defaultClearTestImages, @(x) islogical(x) || isnumeric(x)); % Define optional key
parse(p, varargin{:});                       
clearTestImages = p.Results.clearTestImages; 

[didLoad, triLMSCalFormat, trirgbLinCalFormat, triRGBImage, ...
 diLMSCalFormat, dirgbLinCalFormat, diRGBImage, pathName, outputDir, outputSubdir] = handleLMSFileLoad(img, renderType, imgParams, Disp, clearTestImages);

if didLoad
    fprintf('Found precomputed LMS data for %s\n', img);
    
    % Show image pair of original and dichromat simulation
    figure();
    subplot(1,2,1)
    imagesc(triRGBImage);
    axis square
    title(sprintf('%s | %s | %dx%d', img, renderType, imgParams.m, imgParams.n));
    subtitle('Trichromat RGB image');

    subplot(1,2,2)
    imagesc(diRGBImage);
    axis square
    title(sprintf('%s | %s | %dx%d', img, renderType, imgParams.m, imgParams.n));
    subtitle('Dichromat RGB image');

    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(img,'ishihara')

    [insideColors, outsideColors] = chooseIshiharaColors(renderType,imgParams.plateType,Disp);
        % 1 -> gray with missing cone mod
        % 2 -> background random inside with missing cone mod
        % 3 -> LS background, M inside
        % 4 -> like 2 but constrained between .3 and .7 colors so more room for
        %      modulation

    % Generate plate now that you have the correct colors
    ishiharaRGB = generateIshiharaPlate('74', insideColors, outsideColors,imgParams.m);

    % Gamma corrected image (for visualization)
    triRGBImage = im2double(ishiharaRGB);

    triRGBCalFormat = ImageToCalFormat(triRGBImage);
     
    % Get linear rgb from gamma corrected RGB
    trirgbLinCalFormat = RGB2rgbLin(triRGBCalFormat,Disp);

    % Put modified image into LMS 
    triLMSCalFormat      = Disp.M_rgb2cones * trirgbLinCalFormat;
    
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)

    % Load in the image
    imgFolder = fullfile(outputDir, 'ColorCorrectionImages');
    imgFullPath = fullfile(imgFolder, img);

    imgRGB = im2double(imread(imgFullPath));
    % imgRGB = im2double(imread(img));   

    % Clip ends to ensure in gamut before we process the image
    imgRGB(imgRGB>1) = 1;
    imgRGB(imgRGB<0) = 0;
    
    % Resize the image according to m and n (see buildSetParameters.m)
    % Gamma corrected image (for visualization)
    triRGBImage = imresize(imgRGB, [imgParams.m, imgParams.n]);       

    triRGBCalFormat = ImageToCalFormat(triRGBImage);

    % Clip ends to ensure in gamut before we process the image
    triRGBCalFormat(triRGBCalFormat>1) = 1;
    triRGBCalFormat(triRGBCalFormat<0) = 0;
    
    % Get linear rgb from gamma corrected RGB
    trirgbLinCalFormat = RGB2rgbLin(triRGBCalFormat,Disp);

    % Convert to LMS 
    triLMSCalFormat = Disp.M_rgb2cones * trirgbLinCalFormat;    

else

    %   modType:       Type of cone modulation for hyperspectral cases
    %                  'L', 'M', 'S', or 'rand'
    switch (renderType)
        case 'Deuteranopia'
            modType = 'M';
        case 'Protanopia'
            modType = 'L';
        case 'Tritanopia'
            modType = 'S';
    end

    % Generate gray image with squares
    [~,triLMSCalFormat_plate] = generateGrayImage(imgParams.nSquares,modType,Disp,imgParams);

    % Just grab the version with the isochromatic plate
    triLMSCalFormat = triLMSCalFormat_plate; % do this when you just want to see the isochromatic plate square version (other is just gray)

    % Convert from LMS to linear rgb
    trirgbLinCalFormat = Disp.M_cones2rgb * triLMSCalFormat;
    triRGBCalFormat    = rgbLin2RGB(trirgbLinCalFormat,Disp);

    % Gamma corrected image (for visualization)
    triRGBImage        = CalFormatToImage(triRGBCalFormat,imgParams.m,imgParams.n);

end

 
% Create dichromat simulation
[diLMSCalFormat,  dirgbLinCalFormat] = DichromSimulateLinear(triLMSCalFormat, renderType, Disp);

% What if you use im2double on the gamma corrected version? this actually
% seems to create a gray image, unlike when we perform it on the
% uncorrected linear version
% [diLMSCalFormat,  dirgbLinCalFormat] = DichromSimulateLinear(Disp.M_rgb2cones * im2double(triRGBCalFormat), renderType, Disp);
% figure(); imagesc(CalFormatToImage(dirgbLinCalFormat,imgParams.m,imgParams.n))

% Get dichromat gamma corrected image (for visualization)
diRGBCalFormat = rgbLin2RGB(dirgbLinCalFormat,Disp);
diRGBImage = CalFormatToImage(diRGBCalFormat,imgParams.m,imgParams.n);


% saveLMSData(outputSubdir, img, triLMSCalFormat, trirgbLinCalFormat, triRGBImage, diLMSCalFormat, dirgbLinCalFormat, diRGBImage, Disp, imgParams);
saveLMSData(outputSubdir, triLMSCalFormat, trirgbLinCalFormat, triRGBImage, diLMSCalFormat, dirgbLinCalFormat, diRGBImage, Disp, imgParams)

fprintf('Generated and saved LMS data for %s\n', img);

% Show image pair of original and dichromat simulation
figure();
subplot(1,2,1)
imagesc(triRGBImage);
axis square
title(sprintf('%s | %s | %dx%d', img, renderType, imgParams.m, imgParams.n));
subtitle('Trichromat RGB image');

subplot(1,2,2)
imagesc(diRGBImage);
axis square
title(sprintf('%s | %s | %dx%d', img, renderType, imgParams.m, imgParams.n));
subtitle('Dichromat RGB image')

end