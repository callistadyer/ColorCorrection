%% SCRIPT FOR RUNNING COLOR CORRECTION PROJECT
%
% This script prepares, loads, or generates images for the ColorCorrection
% project. It determines the output directory structure, checks if precomputed
% images exist, loads them if available, or generates new images and saves them.
% The script supports gray w square,  Ishihara plates, and external .png/.jpg images.
%
% HISTORY:
%   06/27/2025  cmd  cleaned up script

%% PARAMETERS TO SET
clc
clear
close all

% Available image options: gray, Ishihara plate, or external images
imageTypes = {'gray','ishihara','flower1.png','flower2.png','flower3.png', ...
    'fruit.png','map.png','painting.png','pool.png','tree.png'};
whichType  = [1:10]; 

% Types of 'setType' usage:
%   gray      - number of squares
%   ishihara  - plate type
%   .png/.jpg - currently unused, but could support e.g. downsampling
setType    = 1;

% Choose dichromat type: 'Deuteranopia', 'Protanopia', 'Tritanopia'
renderType = 'Deuteranopia';

%% Load display calibration
Disp = loadDisplay();

%% Load image parameters
m = 128;
n = 128;

%% Loop through each image and generate all:
for idx = 1:length(whichType)

    typeIdx = whichType(idx);
    imgType = imageTypes{typeIdx};

    imgParams = buildSetParameters(imgType,setType,m,n);

    fprintf('Processing image: %s\n', imgType);


    % Generate LMS/RGB calibration data for the image
    % loadLMSvalues computes and saves trichromat and dichromat LMS and RGB values
    
    % Only clear test images on the first loop
    clearFlag = (idx == 1);

    % Generate LMS/RGB calibration data for the image
    [triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = ...
        loadLMSvalues(imgType, renderType, Disp, imgParams, 'clearTestImages', clearFlag);

end








%% OLD CODE... NOW MOST OF THIS IS DONE INSIDE OF loadLMSvalues
%% DEFINE OUTPUT DIRECTORY
% projectName = 'ColorCorrection';
% myFullPath  = mfilename('fullpath');
% [myPath, myName] = fileparts(myFullPath);
%
% outputDir = getpref(projectName, 'outputDir');
%
% % Determine output subfolder: skip setType folder for png/jpg images
% if endsWith(imageType, {'.png', '.jpg'}, 'IgnoreCase', true)
%     outputSubdir = fullfile(outputDir, 'testImages', dichromatType, imageType);
% else
%     outputSubdir = fullfile(outputDir, 'testImages', dichromatType, imageType, num2str(setType));
% end
%
% % Create output subfolder if it does not exist
% if ~exist(outputSubdir, "dir")
%     mkdir(outputSubdir);
% end

%% CHECK FOR EXISTING IMAGE FILE
% if endsWith(imageType, {'.png', '.jpg'}, 'IgnoreCase', true)
%     imageBaseName = imageType;
% else
%     imageBaseName = [imageType, '.png'];
% end
% imageOutputPath = fullfile(outputSubdir, imageBaseName);
%
% if exist(imageOutputPath, 'file')
%     fprintf('Found precomputed image: %s\n', imageOutputPath);
%     triRGBImage = im2double(imread(imageOutputPath));
%
%     % Optionally load saved calibration data
%     triLMSPath = fullfile(outputSubdir, 'triLMSCalFormat.mat');
%     triRGBPath = fullfile(outputSubdir, 'triRGBCalFormat.mat');
%     dispPath   = fullfile(outputSubdir, 'Disp.mat');
%
%     if exist(triLMSPath, 'file'), load(triLMSPath, 'triLMSCalFormat'); end
%     if exist(triRGBPath, 'file'), load(triRGBPath, 'triRGBCalFormat'); end
%     if exist(dispPath,   'file'), load(dispPath, 'Disp'); end
% else
%% GENERATE AND SAVE NEW IMAGE
% img        = imageType;          % Input image type
% renderType = dichromatType;      % Dichromat simulation type
%
% % Load display calibration
% Disp = loadDisplay(img);
%
% % Generate LMS/RGB calibration data for the image
% [triLMSCalFormat, triRGBCalFormat, Disp] = loadLMSvalues(img, renderType, setType, Disp);

% Convert RGB calibration data to image
% triRGBImage = CalFormatToImage(triRGBCalFormat, Disp.m, Disp.n);
%
% % Save calibration data and image
% save(fullfile(outputSubdir, 'triLMSCalFormat.mat'), 'triLMSCalFormat');
% save(fullfile(outputSubdir, 'triRGBCalFormat.mat'), 'triRGBCalFormat');
% save(fullfile(outputSubdir, 'Disp.mat'), 'Disp');
% imwrite(triRGBImage, imageOutputPath);
% end
