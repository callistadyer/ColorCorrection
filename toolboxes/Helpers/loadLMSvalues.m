function [triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = loadLMSvalues(img,renderType,Disp,imgParams)
% loadLMSvalues  Loads or generates an image and converts to LMS for dichromat simulation
%
% Syntax:
%   [triLMSCalFormat,diLMSCalFormat] = loadLMSvalues(img,renderType,modType,nSquares,constraintWL,plateType,Disp)
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
[triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = loadLMSvalues('ishihara','Deuteranopia',Disp,imgParams);

% Check that behavior has not changed since we declared it good.
if (abs(sum(testLMS(:)) - 525.8024)/525.8024 > 1e-4)
    error('No longer get same LMS image returned by loadLMSValues');
end
%}



%%%%%%%%%%%%%%%%%%%%%% Dealing with save directory %%%%%%%%%%%%%%%%%%%%%%
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');


% Make this a key value pair
clearTestImages = true;  
testImagesDir = fullfile(outputDir, 'testImages');
% Clear 'testImages' folder if requested
if exist('clearTestImages', 'var') && clearTestImages
    if exist(testImagesDir, 'dir')
        fprintf('Clearing entire folder: %s\n', testImagesDir);
        rmdir(testImagesDir, 's');  
    end
    mkdir(testImagesDir);  % recreate empty testImages folder
end


% Check to see if the image already exists
% Determine output subdirectory
if endsWith(img, {'.png', '.jpg'}, 'IgnoreCase', true)
    outputSubdir = fullfile(outputDir, 'testImages', renderType, img);
else
    outputSubdir = fullfile(outputDir, 'testImages', renderType, img, num2str(imgParams.setType));
end
if ~exist(outputSubdir, "dir")
    mkdir(outputSubdir);
end

% Save path after testImages. Can use this later when creating parallel
% file structure as the original images... for the transformed images
idx = strfind(outputSubdir, ['testImages' filesep]);
pathName = outputSubdir(idx + length('testImages') + 1 : end);

% Determine expected output file paths
triLMSPath    = fullfile(outputSubdir, 'triLMSCalFormat.mat');
trirgbLinPath = fullfile(outputSubdir, 'trirgbLinCalFormat.mat');
triRGBPath    = fullfile(outputSubdir, 'triRGBCalFormat.mat');

diLMSPath    = fullfile(outputSubdir, 'diLMSCalFormat.mat');
dirgbLinPath = fullfile(outputSubdir, 'dirgbLinCalFormat.mat');
diRGBPath    = fullfile(outputSubdir, 'diRGBCalFormat.mat');

dispPath      = fullfile(outputSubdir, 'Disp.mat');
imgParamsPath = fullfile(outputSubdir, 'imgParams.mat');
if endsWith(img, {'.png', '.jpg'}, 'IgnoreCase', true)
    imageBaseName = img;
else
    imageBaseName = [img, '.png'];
end

imageOutputPath = fullfile(outputSubdir, imageBaseName);

% If all outputs exist, load them and return
if exist(triLMSPath, 'file') && exist(trirgbLinPath, 'file') && exist(triRGBPath, 'file')...
        && exist(diLMSPath, 'file') && exist(dirgbLinPath, 'file') && exist(diRGBPath, 'file')...
        && exist(dispPath, 'file') && exist(imgParamsPath, 'file')
    fprintf('Found precomputed LMS data for %s', img);
    load(triLMSPath,    'triLMSCalFormat');
    load(trirgbLinPath, 'trirgbLinCalFormat');
    load(triRGBPath,    'triRGBImage');

    load(diLMSPath,     'diLMSCalFormat');
    load(dirgbLinPath,  'dirgbLinCalFormat');
    load(diRGBPath,     'diRGBImage');

    load(dispPath,      'Disp');
    load(imgParamsPath, 'imgParams');

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


    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% If the image doesn't exist, continue %%%%%%%
%   modType:       Type of cone modulation for hyperspectral cases
%                  'L', 'M', 'S', or 'rand'
switch (renderType)
    case 'Deuteranopia'
        modType = 'M';
        constraintWL = 585;
    case 'Protanopia'
        modType = 'L';
        error('ERROR: you need to set up constraint wavelength for Protanopia case')
    case 'Tritanopia'
        modType = 'S';
        error('ERROR: you need to set up constraint wavelength for Tritanopia case')
end

if strcmp(img,'ishihara')

    [insideColors, outsideColors] = chooseIshiharaColors(renderType,imgParams.plateType,Disp);
        % 1 -> gray with missing cone mod
        % 2 -> background random inside with missing cone mod
        % 3 -> LS background, M inside
        % 4 -> like 2 but constrained between .3 and .7 colors so more room for
        %      modulation

    % Generate plate now that you have the correct colors
    triRGBImage = generateIshiharaPlate('74', insideColors, outsideColors,imgParams.m);

    % Gamma corrected image (for visualization)
    % triRGBImage = im2double(ishiharaRGB);

    triRGBCalFormat = ImageToCalFormat(triRGBImage);
     
    % Get linear rgb from gamma corrected RGB
    trirgbLinCalFormat = RGB2rgbLin(triRGBCalFormat,Disp,imgParams);

    % Put modified image into LMS 
    triLMSCalFormat      = Disp.M_rgb2cones * trirgbLinCalFormat;
    
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)

    % Load in the image 
    imgRGB = im2double(imread(img));   

    % Resize the image according to m and n (see buildSetParameters.m)
    % Gamma corrected image (for visualization)
    triRGBImage = imresize(imgRGB, [imgParams.m, imgParams.n]);       

    % Get linear rgb from gamma corrected RGB
    imgrgbLin = RGB2rgbLin(triRGBImage,Disp);

    % Cal format linear rgb 
    trirgbLinCalFormat = ImageToCalFormat(imgrgbLin); 

    % Clip ends to ensure in gamut before we process the image
    trirgbLinCalFormat(trirgbLinCalFormat>1) = 1;
    trirgbLinCalFormat(trirgbLinCalFormat<0) = 0;

    % Convert to LMS 
    triLMSCalFormat = Disp.M_rgb2cones * trirgbLinCalFormat;    

else
    
    % Generate gray image with squares
    [triLMSCalFormat,triLMSCalFormat_plate] = generateGrayImage(img,setParams.nSquares,modType,Disp);

    % Just grab the version with the isochromatic plate
    triLMSCalFormat = triLMSCalFormat_plate; % do this when you just want to see the isochromatic plate square version (other is just gray)

    % Convert from LMS to linear rgb
    trirgbLinCalFormat = Disp.M_cones2rgb * triLMSCalFormat;
    triRGBCalFormat    = rgbLin2RGB(trirgbLinCalFormat,Disp,imgParams);

    % Gamma corrected image (for visualization)
    triRGBImage        = CalFormatToImage(triRGBCalFormat,imgParams.m,imgParams.n);

end

% Create dichromat simulation
[diLMSCalFormat,  dirgbLinCalFormat] = DichromSimulateLinear(triLMSCalFormat, renderType, Disp);

% Get dichromat gamma corrected image (for visualization)
diRGBCalFormat = rgbLin2RGB(dirgbLinCalFormat,Disp,imgParams);
diRGBImage = CalFormatToImage(diRGBCalFormat,imgParams.m,imgParams.n);

% Convert and save image
imwrite(triRGBImage, imageOutputPath);

% Save trichromat values
save(triLMSPath, 'triLMSCalFormat');
save(trirgbLinPath, 'trirgbLinCalFormat');
save(triRGBPath, 'triRGBImage') % should this be cal format or image?

% Save dichromat values
save(diLMSPath, 'diLMSCalFormat');
save(dirgbLinPath, 'dirgbLinCalFormat');
save(diRGBPath, 'diRGBImage') % should this be cal format or image?

% Save Display parameters
save(dispPath,  'Disp');

% Save image parameters
save(imgParamsPath,  'imgParams');

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