function [LMSCalFormat, rgbLinCalFormat, LMSCalFormatRendered, rgbLinCalFormatRendered, Disp, imgParams, pathName] = ...
    colorCorrectionGenerateImages(imgType, setType, m, n, dichromatType,clearFlag)
% prepareColorCorrectionImage  Generate LMS/RGB values for an input image used in color correction.
%
% Syntax:
%   [triLMSCalFormat, trirgbLinCalFormat, diLMSCalFormat, dirgbLinCalFormat, pathName] = ...
%       prepareColorCorrectionImage(imgType, setType, m, n, dichromatType)
%
% Description:
%   This function prepares or loads a stimulus image for use in the ColorCorrection
%   project. It loads display calibration, sets image-specific parameters,
%   and generates LMS/RGB values for both trichromat and dichromat observers.
%
% Inputs:
%   imgType         - String specifying image type. Default is 'flower1.png'.
%                     Options:
%                       'gray', 'ishihara', 'flower1.png', 'flower2.png',
%                       'flower3.png', 'fruit.png', 'map.png', 'painting.png',
%                       'pool.png', 'tree.png'
%   setType         - Type-specific index (e.g., 1 for gray squares, plate number for Ishihara).
%                     Default: 1
%   m               - Image height in pixels (used by buildSetParameters). Default: 32
%   n               - Image width in pixels. Default: 32
%   dichromatType   - Type of dichromacy. Options: 'Deuteranopia', 'Protanopia', 'Tritanopia'
%                     Default: 'Deuteranopia'
%   clearFlag       - Boolean. If true, clears any precomputed test images before regeneration.
%                     Default: true
% Outputs:
%   LMSCalFormat               - LMS values for the trichromatic rendering
%   rgbLinCalFormat            - Linear RGB values for the trichromatic rendering
%   LMSCalFormatRendered       - LMS values for the simulated dichromatic rendering
%   rgbLinCalFormatRendered    - Linear RGB values for the simulated dichromatic rendering
%   Disp                       - Display parameters
%   pathName                   - Directory where results are stored
% 
% History:
%   07/30/2025  cmd  Converted script into function with defaults

    % Set defaults if arguments are missing
    if nargin < 1 || isempty(imgType)
        imgType = 'flower1.png';
    end
    if nargin < 2 || isempty(setType)
        setType = 1;
    end
    if nargin < 3 || isempty(m)
        m = 32;
    end
    if nargin < 4 || isempty(n)
        n = 32;
    end
    if nargin < 5 || isempty(dichromatType)
        dichromatType = 'Deuteranopia';
    end

    % Load display calibration
    Disp = loadDisplay();

    % Construct image parameter struct
    imgParams = buildSetParameters(imgType, setType, m, n);

    fprintf('Preparing image: %s\n', imgType);

    % Call loadLMSvalues to get calibrated data
    [LMSCalFormat, rgbLinCalFormat, LMSCalFormatRendered, rgbLinCalFormatRendered, pathName] = ...
        loadLMSvalues(imgType, dichromatType, Disp, imgParams, 'clearTestImages', clearFlag);

end


% %% SCRIPT FOR RUNNING COLOR CORRECTION PROJECT
% %
% % This script prepares, loads, or generates images for the ColorCorrection
% % project. It determines the output directory structure, checks if precomputed
% % images exist, loads them if available, or generates new images and saves them.
% % The script supports gray w square,  Ishihara plates, and external .png/.jpg images.
% %
% % HISTORY:
% %   06/27/2025  cmd  cleaned up script
% 
% %% PARAMETERS TO SET
% clc
% clear
% close all
% 
% % Available image options: gray, Ishihara plate, or external images
% imageTypes = {'gray','ishihara','flower1.png','flower2.png','flower3.png', ...
%     'fruit.png','map.png','painting.png','pool.png','tree.png'};
% whichType  = [3]; 
% 
% % Types of 'setType' usage:
% %   gray      - number of squares
% %   ishihara  - plate type
% %   .png/.jpg - currently unused, but could support e.g. downsampling
% setType    = 1;
% 
% % Choose dichromat type: 'Deuteranopia', 'Protanopia', 'Tritanopia'
% dichromatType = 'Deuteranopia';
% 
% %% Load display calibration
% Disp = loadDisplay();
% 
% %% Load image parameters
% m = 32;
% n = 32;
% 
% %% Loop through each image and generate all:
% for idx = 1:length(whichType)
% 
%     typeIdx = whichType(idx);
%     imgType = imageTypes{typeIdx};
% 
%     imgParams = buildSetParameters(imgType,setType,m,n);
% 
%     fprintf('Processing image: %s\n', imgType);
% 
% 
%     % Generate LMS/RGB calibration data for the image
%     % loadLMSvalues computes and saves trichromat and dichromat LMS and RGB values
% 
%     % Only clear test images on the first loop
%     clearFlag = (idx == 1);
% 
%     % Generate LMS/RGB calibration data for the image
%     [triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = ...
%         loadLMSvalues(imgType, dichromatType, Disp, imgParams, 'clearTestImages', clearFlag);
% 
% end
% 
