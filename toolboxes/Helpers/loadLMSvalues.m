function [triLMSCalFormat,trirgbLinCalFormat,pathName] = loadLMSvalues(img,renderType,Disp,imgParams)
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
imgParams = buildSetParameters('ishihara',1,128,128)
testLMS = loadLMSvalues('ishihara','Deuteranopia',Disp,imgParams);

% Check that behavior has not changed since we declared it good.
if (abs(sum(testLMS(:)) - 525.8024)/525.8024 > 1e-4)
    error('No longer get same LMS image returned by loadLMSValues');
end
%}



%%%%%%%%%%%%%%%%%%%%%% Dealing with save directory %%%%%%%%%%%%%%%%%%%%%%
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');

% Check to see if the image already exists
% Determine output subdirectory
if endsWith(img, {'.png', '.jpg'}, 'IgnoreCase', true)
    outputSubdir = fullfile(outputDir, 'testImages', renderType, img);
    pathName = outputSubdir;
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
triLMSPath = fullfile(outputSubdir, 'triLMSCalFormat.mat');
triRGBPath = fullfile(outputSubdir, 'triRGBCalFormat.mat');
dispPath   = fullfile(outputSubdir, 'Disp.mat');
if endsWith(img, {'.png', '.jpg'}, 'IgnoreCase', true)
    imageBaseName = img;
else
    imageBaseName = [img, '.png'];
end

imageOutputPath = fullfile(outputSubdir, imageBaseName);

% If all outputs exist, load them and return
if exist(triLMSPath, 'file') && exist(triRGBPath, 'file') && exist(dispPath, 'file') && exist(imageOutputPath, 'file')
    fprintf('Found precomputed LMS data for %s', img);
    load(triLMSPath, 'triLMSCalFormat');
    load(triRGBPath, 'trirgbLinCalFormat');
    load(dispPath,   'Disp');
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

    % 1 -> gray with missing cone mod
    % 2 -> background random inside with missing cone mod
    % 3 -> LS background, M inside
    % 4 -> like 2 but constrained between .3 and .7 colors so more room for
    %      modulation

    [insideColors, outsideColors] = chooseIshiharaColors(renderType,imgParams.plateType,Disp);

    % Generate plate now that you have the correct colors
    ishiharaRGB = generateIshiharaPlate('74', insideColors, outsideColors,imgParams.m);
    ishiharaRGB = im2double(ishiharaRGB);
    triRGBImage = ishiharaRGB;
    % Get linear rgb from gamma corrected RGB
    trirgbLinCalFormat = RGB2rgbLin(ishiharaRGB,Disp,imgParams);

    % Plot modified RGB Image 
    % figure();imagesc(ishiharaRGB)
    % axis square;

    % Put modified image into LMS 
    triLMSCalFormat      = Disp.M_rgb2cones * trirgbLinCalFormat;
    
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)

    imgRGB = im2double(imread(img));           
    imgRGB = imresize(imgRGB, [imgParams.m, imgParams.n]);       

    imgrgbLin = RGB2rgbLin(imgRGB,Disp);

    % Cal format RGB 
    trirgbLinCalFormat = ImageToCalFormat(imgrgbLin); 
    trirgbLinCalFormat(trirgbLinCalFormat>1) = 1;
    trirgbLinCalFormat(trirgbLinCalFormat<0) = 0;
    % Convert to LMS 
    triLMSCalFormat = Disp.M_rgb2cones * trirgbLinCalFormat;

    % [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, Disp.grayLMS,  constraintWL, renderType, Disp);
    
    % check
    % rgb1 = inv(Disp.M_rgb2cones) * diLMSCalFormat;
    % imageDi = CalFormatToImage(rgb1, Disp.m, Disp.n);
    % figure();imagesc(imageDi);

else
    % I think t_renderHyperspectralImage is only used in order to create
    % LMS values for the gray image with square isochromatic plates added
    % on. Maybe you can simplify this? Not sure. 
    
    [triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate] = t_renderHyperspectralImage(img,setParams.nSquares,modType,Disp);
    triLMSCalFormat = triLMSCalFormat_plate; % do this when you just want to see the isochromatic plate square version (other is just gray)
    trirgbLinCalFormat = Disp.M_cones2rgb * triLMSCalFormat;
    % diLMSCalFormat  = diLMSCalFormat_plate;
end

% Convert and save image
imwrite(triRGBImage, imageOutputPath);

% Save trichromat values
save(triLMSPath, 'triLMSCalFormat');
save(triRGBPath, 'trirgbLinCalFormat');
save(dispPath,   'Disp');

% Save dichromat values
%
%


fprintf('Generated and saved LMS data for %s\n', img);


% Show image pair of original and dichromat simulation
%
%

end