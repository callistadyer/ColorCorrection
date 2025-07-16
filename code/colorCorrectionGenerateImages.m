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

