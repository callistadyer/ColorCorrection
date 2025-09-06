%% How to use daltonize and compute object oriented methods
% This script sets up and runs a daltonization sweep across target information levels.
% It visualizes both the trichromatic renderings and simulated dichromatic views.

% clear;
% close all;

%% Generate input image
%
% 'Disp' is defined here
%                       'gray', 'ishihara', 'flower1.png', 'flower2.png',
%                       'flower3.png', 'fruit.png', 'map.png', 'painting.png',
%                       'pool.png', 'tree.png'
imgType = 'ishi45.png';
setType = 1;
dichromatType = 'Deuteranopia';
clearFlag     = 0;
m = 64;
n = 64;

[LMSCalFormat, rgbLinCalFormat, LMSCalFormatRendered, rgbLinCalFormatRendered, Disp, imgParams, pathName] = ...
    colorCorrectionGenerateImages(imgType, setType, m, n, dichromatType,clearFlag);

%% Define objective functions
%
% Info: Measures useful information retained after correction
% infoFcn = @computeInfo_LMdifference;
% infoFcn = @computeInfo_regress;
infoFcn = @computeInfo_Wade;
infoParams = struct();  

% infoParams = struct( ...
%     'predictingWhat',     'M', ...          % options: 'M' | 'L-M' | 'L,M, and S'
%     'predictingFromWhat', 'L and S' ...     % options: 'L and S' | 'deltaL and deltaS' | 'deltaL+M and delta S'
% );

if isequal(func2str(infoFcn), 'computeInfo_regress') && isempty(fieldnames(infoParams))
    error(['If infoFcn is computeInfo_regress, then infoParams cannot be empty.' newline ...
        'Please specify what information you are predicting (predictingWhat) ' ...
        'as well as what information you are using to make the prediction (predictingFromWhat).' newline newline ...
        'Valid options are:' newline ...
        '   predictingWhat:     ''M''   | ''L-M''   | ''L,M, and S''' newline ...
        '   predictingFromWhat: ''L and S''   | ''deltaL and deltaS''   | ''deltaL+M and delta S''' newline newline ...
        'Example call:' newline ...
        '   infoParams = struct(''predictingWhat'', ''M'', ''predictingFromWhat'', ''L and S'');' ]);
end

% Distortion: Measures image distortion after correction
distortionFcn = @computeDistortion_squared;
distortionParams = struct();  % Add fields if needed

% Render: Simulates dichromatic vision
renderFcn = @DichromRenderLinear;
renderParams = struct();  % Add fields if needed

%% Set up daltonizer object
theDaltonizer = daltonize( ...
    infoFcn, infoParams, ...
    distortionFcn, distortionParams, ...
    renderFcn, renderParams, ...
    Disp);

%% Set sweep parameters and run optimization
nSteps = 30;  % Number of interpolation steps in info space

% Run sweep: for e3ach info target, compute optimized transformation
[LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
    LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
    transformRGBmatrixSweep, targetInfoNormalized, infoNormalized, distortionNormalized] = ...
    theDaltonizer.computeInfoSweep(LMSCalFormat, imgParams, dichromatType, nSteps, pathName);


%% Visualize trichromat renderings
figure('Name', 'Trichromat Renderings');

nCols = ceil(sqrt(nSteps));       
nRows = ceil(nSteps / nCols);

for i = 1:nSteps
    subplot(nRows, nCols, i);
    triRGBCalFormat{i} = rgbLin2RGB(rgbLinDaltonizedCalFormatSweep{i},Disp);
    rgbTri = CalFormatToImage(triRGBCalFormat{i}, imgParams.m, imgParams.n);
    imshow(rgbTri);
    title(sprintf('Tri Info %.2f', targetInfoNormalized(i)), 'FontSize', 8);
end

sgtitle(sprintf('Trichromat Renderings — Daltonization Sweep (%s)', dichromatType), ...
    'FontWeight', 'bold');

%% Visualize rendered dichromat views
figure('Name', 'Dichromat Simulations');
for i = 1:nSteps
    subplot(nRows, nCols, i);
    diRGBCalFormat{i} = rgbLin2RGB(rgbLinDaltonizedRenderedCalFormatSweep{i},Disp);
    rgbDi = CalFormatToImage(diRGBCalFormat{i}, imgParams.m, imgParams.n);
    imshow(rgbDi);
    title(sprintf('Dich Sim %.2f', targetInfoNormalized(i)), 'FontSize', 8);
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
