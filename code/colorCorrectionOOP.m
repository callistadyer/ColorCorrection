%% How to use daltonize and compute object oriented methods
% This script sets up and runs a daltonization sweep across target information levels.
% It visualizes both the trichromatic renderings and simulated dichromatic views.

% clear;
% close all;

%% Generate input image

% 'Disp' is defined here
colorCorrectionGenerateImages;

%% Define objective functions

% Info: Measures useful information retained after correction
infoFcn = @computeInfo_regress;
infoParams = struct();  

% Distortion: Measures image distortion after correction
distortionFcn = @computeDistortion_squared;
distortionParams = struct();  % Add fields if needed

% Render: Simulates dichromatic vision
renderFcn = @DichromSimulateLinear;
renderParams = struct();  % Add fields if needed

%% Set up daltonizer object
theDaltonizer = daltonize( ...
    infoFcn, infoParams, ...
    distortionFcn, distortionParams, ...
    renderFcn, renderParams, ...
    Disp);

%% Set sweep parameters and run optimization
nSteps = 30;  % Number of interpolation steps in info space
dichromatType = 'Deuteranopia';  % 'Protanopia', 'Deuteranopia', 'Tritanopia'

% Run sweep: for each info target, compute optimized transformation
[triLMSCalFormatOpt, trirgbLinCalFormatOpt, ...
 diLMSCalFormatOpt, dirgbLinCalFormatOpt, ...
 infoVals, infoNormVals, Tmatrices, targetInfoVals] = ...
    theDaltonizer.computeInfoSweep(triLMSCalFormat, imgParams, dichromatType, nSteps);

%% Visualize trichromat renderings
figure('Name', 'Trichromat Renderings');

nCols = ceil(sqrt(nSteps));       
nRows = ceil(nSteps / nCols);

for i = 1:nSteps
    subplot(nRows, nCols, i);
    triRGBCalFormat{i} = rgbLin2RGB(trirgbLinCalFormatOpt{i},Disp);
    rgbTri = CalFormatToImage(triRGBCalFormat{i}, imgParams.m, imgParams.n);
    imshow(rgbTri);
    title(sprintf('Tri Info %.2f', targetInfoVals(i)), 'FontSize', 8);
end

sgtitle(sprintf('Trichromat Renderings — Daltonization Sweep (%s)', dichromatType), ...
    'FontWeight', 'bold');

%% Visualize simulated dichromat views
figure('Name', 'Dichromat Simulations');

for i = 1:nSteps
    subplot(nRows, nCols, i);
    diRGBCalFormat{i} = rgbLin2RGB(dirgbLinCalFormatOpt{i},Disp);
    rgbDi = CalFormatToImage(diRGBCalFormat{i}, imgParams.m, imgParams.n);
    imshow(rgbDi);
    title(sprintf('Dich Sim %.2f', targetInfoVals(i)), 'FontSize', 8);
end

sgtitle(sprintf('Dichromat Simulated View — Daltonization Sweep (%s)', dichromatType), ...
    'FontWeight', 'bold');

%% STEP 7: Single compute call 
% Uncomment and configure if you want to run a single optimization directly.
% The variables below (obj, LMSCalFormat, useLambdaOrTargetInfo, etc.)
% must be defined appropriately if using this block.

% [LMSDaltonizedCalFormat, LMSDaltonizedRenderedCalFormat] = compute(obj, ...
%     LMSCalFormat, imgParams, dichromatType, ...
%     useLambdaOrTargetInfo, lambdaOrTargetInfo, varargin);
