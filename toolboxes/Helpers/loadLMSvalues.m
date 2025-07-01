function [triLMSCalFormat,trirgbLinCalFormat] = loadLMSvalues(img,renderType,Disp)
% loadLMSvalues  Loads or generates an image and converts to LMS for dichromat simulation
%
% Syntax:
%   [triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues(img,renderType,modType,nSquares,constraintWL,plateType,Disp)
%
% Inputs:
%   img:              Either 'ishihara' or a filename ('.png', '.jpg') or a hyperspectral identifier
%   renderType:       Type of dichromacy to simulate
%                         'Deuteranopia', 'Protanopia', or 'Tritanopia'
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
Disp = loadDisplay('ishihara');
testLMS = loadLMSvalues('ishihara','Deuteranopia',1,Disp);

% Check that behavior has not changed since we declared it good.
if (abs(sum(testLMS(:)) - 525.8024)/525.8024 > 1e-4)
    error('No longer get same LMS image returned by loadLMSValues');
end
%}

projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');


%%%%%%% Check to see if the image already exists %%%%%%%
% Determine output subdirectory
if endsWith(img, {'.png', '.jpg'}, 'IgnoreCase', true)
    outputSubdir = fullfile(outputDir, 'testImages', renderType, img);
else
    outputSubdir = fullfile(outputDir, 'testImages', renderType, img, num2str(Disp.setType));
end
if ~exist(outputSubdir, "dir")
    mkdir(outputSubdir);
end

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

%%%%%%% If the image doesn't exist, do this: %%%%%%%
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
setParams = Disp.setParams;

if strcmp(img,'ishihara')

    % 1 -> gray with missing cone mod
    % 2 -> background random inside with missing cone mod
    % 3 -> LS background, M inside
    % 4 -> like 2 but constrained between .3 and .7 colors so more room for
    %      modulation

    [insideColors, outsideColors] = chooseIshiharaColors(renderType,setParams.plateType,Disp);

    % Generate plate now that you have the correct colors
    ishiharaRGB = generateIshiharaPlate('74', insideColors, outsideColors,Disp.m);
    ishiharaRGB = im2double(ishiharaRGB);

    % Get linear rgb from gamma corrected RGB
    ishiharargbLin = RGB2rgbLin(ishiharaRGB,Disp);

    % Plot modified RGB Image 
    % figure();imagesc(ishiharaRGB)
    % axis square;

    % Put modified image into LMS 
    trirgbLinCalFormat   = ImageToCalFormat(ishiharargbLin);
    triLMSCalFormat      = Disp.M_rgb2cones * trirgbLinCalFormat;
    % 
    % M_plane = triLMSCalFormat(2,:);
    % % Assuming M_plane is 1 x N
    % N = length(M_plane);
    % 
    % % Fill in L and S planes with mean of M (or another estimate)
    % L_plane = mean(M_plane) * zeros(1, N);
    % S_plane = mean(M_plane) * zeros(1, N);
    % LMS_full = [L_plane; M_plane; S_plane];
    % M_cones2rgb = inv(Disp.M_rgb2cones);
    % RGB = M_cones2rgb * LMS_full;
    % 
    % rgbImage = reshape(RGB', Disp.m, Disp.n, 3);
    % figure();
    % imshow(rgbImage);



    % Run the modulated image through the linear dichromat simulation
    % [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, Disp.grayLMS,  constraintWL, renderType, Disp);

    % Check
    % rgb1 = inv(Disp.M_rgb2cones) * diLMSCalFormat;
    % image = CalFormatToImage(rgb1,Disp.m,Disp.n);
    % figure();imagesc(image); axis square;
    
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)

    imgRGB = im2double(imread(img));           
    imgRGB = imresize(imgRGB, [128, 128]);         % Resize to 128x128 pixels

    [rows, cols, ~] = size(imgRGB);
    Disp.m = cols;
    Disp.n = rows;

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
    % Get trichromatic (LMS) image
    % Disp.m = 32;
    % Disp.n = 32;

    % I think t_renderHyperspectralImage is only used in order to create
    % LMS values for the gray image with square isochromatic plates added
    % on. Maybe you can simplify this? Not sure. 
    
    [triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate] = t_renderHyperspectralImage(img,setParams.nSquares,modType,Disp);
    triLMSCalFormat = triLMSCalFormat_plate; % do this when you just want to see the isochromatic plate square version (other is just gray)
    trirgbLinCalFormat = Disp.M_cones2rgb * triLMSCalFormat;
    % diLMSCalFormat  = diLMSCalFormat_plate;
end

% Convert and save image
triRGBImage = CalFormatToImage(trirgbLinCalFormat, Disp.m, Disp.n);
imwrite(triRGBImage, imageOutputPath);

% Save trichromat values
save(triLMSPath, 'triLMSCalFormat');
save(triRGBPath, 'trirgbLinCalFormat');
save(dispPath,   'Disp');

% Save dichromat values
%
%


fprintf('Generated and saved LMS data for %s\n', img);


end