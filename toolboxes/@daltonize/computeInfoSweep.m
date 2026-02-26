function [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
    LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
    transformRGBmatrixSweep, targetInfoNormalized, infoNormalized, distortionNormalized] = computeInfoSweep(obj,...
    LMSCalFormat, imgParams, dichromatType, nSteps,pathName,sweepAxis)
% computeInfoSweep  Sweep through target info values and optimize color correction.
%
% Syntax:
%    [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
%           LMSDaltonizedRenderedCalFiormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
%           transformRGBmatrixSweep, infoNormalized, distortionNormalized, targetInfoNormalized] = computeInfoSweep(obj,...
%           LMSCalFormat, imgParams, dichromatType, nSteps)
% Description:
%   This method sweeps between the info values produced by lambda=0 and lambda=1,
%   then optimizes the transformation for each intermediate targetInfo.
%
% Inputs:
%   LMSCalFormat   - 3 x N LMS values of the input image
%   dichromatType  - type of dichromacy
%                            'Protaniopia'
%                            'Deuteranopia'
%                            'Tritanopia'
%   imgParams      - Struct with image-related params
%   nSteps         - Number of target info steps to interpolate (default: 10).
%   pathName       - where the original photo is stored:
%                   e.g., pathName = 'Deuteranopia/flower2.png/s1_m32_n32'
%
% Outputs:
%   LMSDaltonizedCalFormatSweep            - Cell of optimized LMS images (this is the size of 1xnSteps)
%   rgbLinDaltonizedCalFormatSweep         - Cell of optimized RGB images
%   LMSDaltonizedRenderedCalFormatSweep    - Cell of LMSDaltonizedCalFormatSweep rendered for a dichromat
%   rgbLinDaltonizedRenderedCalFormatSweep - Cell of rgbLinDaltonizedCalFormatSweep rendered for a dichromat
%   transformRGBmatrixSweep                - Cell of 3x3 transformation matrices
%   targetInfoNormalized                   - Vector of target info values used
%   infoNormalized                         - Cell of normalized info values
%   distortionNormalized                   - Cell of normalized distortion values
%

if isempty(nSteps)
    nSteps = 10;
end

if ~exist('sweepAxis','var') || isempty(sweepAxis)
    sweepAxis = 'info';
end

% Check to see if this output already exists, and then load it
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');
saveBase    = fullfile(outputDir, 'testImagesTransformed');

% Build the metric folder 
metricFolder = buildMetricFolderName(obj.infoFcn, obj.infoParams, obj.distortionFcn);

runFolder = sprintf('%dsteps', nSteps);
saveSubdir = fullfile(saveBase, pathName, metricFolder, runFolder);
saveFile   = fullfile(saveSubdir, 'sweepOutputs.mat');

if exist(saveFile, 'file')
    fprintf('[computeInfoSweep] Loading old sweep from: %s\n', saveFile);
    loaded  = load(saveFile);
    outputs = loaded.outputs;

    if isfield(outputs{1}, 'useLambdaOrTargetInfo')
        useLambdaOrTargetInfo = outputs{1}.useLambdaOrTargetInfo;  % 'lambda' or 'targetInfo'
    else
        useLambdaOrTargetInfo = 'targetInfo'; % If it doesn't specify, it's gonna be target info because I only added this after incorporating lambda
    end

    % Extract these important outputs
    nSteps = numel(outputs);
    LMSDaltonizedCalFormatSweep            = cell(1,nSteps);
    rgbLinDaltonizedCalFormatSweep         = cell(1,nSteps);
    LMSDaltonizedRenderedCalFormatSweep    = cell(1,nSteps);
    rgbLinDaltonizedRenderedCalFormatSweep = cell(1,nSteps);
    transformRGBmatrixSweep                = cell(1,nSteps);
    targetInfoNormalized                   = zeros(1,nSteps);
    infoNormalized                         = cell(1,nSteps);
    distortionNormalized                   = cell(1,nSteps);
    useLambdaOrTargetInfo                  = cell(1,nSteps);

    for i = 1:nSteps
        LMSDaltonizedCalFormatSweep{i}            = outputs{i}.LMSDaltonizedCalFormat;
        rgbLinDaltonizedCalFormatSweep{i}         = outputs{i}.rgbLinDaltonizedCalFormat;
        LMSDaltonizedRenderedCalFormatSweep{i}    = outputs{i}.LMSDaltonizedRenderedCalFormat;
        rgbLinDaltonizedRenderedCalFormatSweep{i} = outputs{i}.rgbLinDaltonizedRenderedCalFormat;
        transformRGBmatrixSweep{i}                = outputs{i}.transformRGBmatrix;
        targetInfoNormalized(i)                   = outputs{i}.targetInfoNormalized;
        infoNormalized{i}                         = outputs{i}.infoNormalized;
        distortionNormalized{i}                   = outputs{i}.distortionNormalized;
        useLambdaOrTargetInfo{i}                  = outputs{1}.useLambdaOrTargetInfo;

    end
    return;
end

% If it doesn't already exist, go ahead and compute the transformation as normal:

% Get display
Disp = obj.Disp;
% Get the contrast LMS values
LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

% Get the dichromat renderings of the original, which is used for the normalizer
[calFormatLMS_prot, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Protanopia', Disp);
[calFormatLMS_deut, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Deuteranopia', Disp);
[calFormatLMS_trit, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Tritanopia', Disp);
LMSCalFormat_new = [calFormatLMS_prot(1,:); calFormatLMS_deut(2,:); calFormatLMS_trit(3,:)];
LMSContrastCalFormat_new = (LMSCalFormat_new - Disp.grayLMS) ./ Disp.grayLMS;

% Normalize distortion and info
normalizerValueToGetRawValue = 1;

infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
% distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, normalizerValueToGetRawValue, obj.distortionParams);
distortionNormalizer = obj.distortionFcn(LMSCalFormat, LMSCalFormat_new, imgParams, normalizerValueToGetRawValue, obj.distortionParams);

% Get info values for lambda = 0 and lambda = 1
[~,~,~,info_0,infoNormalized_0,distortion0, distortionNormalized0] = colorCorrectionOptimize("lambda", 0, [], [], LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

[~,~,~,info_1,infoNormalized_1,distortion1, distortionNormalized1] = colorCorrectionOptimize("lambda", 1, [], [], LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

% Get target info values interpolated between lambdas 0 and 1
targetInfoNormalized = linspace(infoNormalized_0, infoNormalized_1, nSteps);

% Preallocate outputs
LMSDaltonizedCalFormatSweep     = cell(1, nSteps);
rgbLinDaltonizedCalFormatSweep  = cell(1, nSteps);
LMSDaltonizedRenderedCalFormatSweep      = cell(1, nSteps);
rgbLinDaltonizedRenderedCalFormatSweep   = cell(1, nSteps);
infoNormalized           = cell(1, nSteps);
distortionNormalized     = cell(1, nSteps);
transformRGBmatrixSweep  = cell(1, nSteps);
targetInfoVsAchievedInfo = cell(1, nSteps);

%% Sweep through target infos
T_I = eye(3);
T_prev = eye(3);

% This part lets us sweep through lambda OR through target info values. I
% want to sweep through info values because it works. I want to sweep
% through lambda values to show people that it doesn't work. 
switch lower(sweepAxis)
    case 'info'
        xVec = targetInfoNormalized;                
        useLambdaOrTargetInfo = "targetInfo";
    case 'lambda'
        xVec = linspace(0,1,nSteps);               
        useLambdaOrTargetInfo = "lambda";
end

for i = 1:nSteps
    % thisTargetInfo = targetInfoNormalized(i);
    thisX = xVec(i);

    if strcmpi(useLambdaOrTargetInfo,'targetInfo')
        % info sweep: lambda=[], targetInfo=thisX
        lamArg = []; tgtInfoArg = thisX; tgtDistArg = [];
    else % 'lambda'
        % lambda sweep: lambda=thisX, targetInfo=[]
        lamArg = thisX; tgtInfoArg = []; tgtDistArg = [];
    end

    % Optimize from identity starting point
    [LMS_TI, rgbLin_TI, T_TI, info_TI, normInfo_TI, distortion_TI, normDistortion_TI] = colorCorrectionOptimize( ...
        useLambdaOrTargetInfo, lamArg, tgtInfoArg, tgtDistArg, ...
        LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        infoNormalizer, distortionNormalizer, Disp, 'T_init', T_I);

    % Optimize from T_prev starting point
    [LMS_Tprev, rgbLin_Tprev, T_Tprev, info_Tprev, normInfo_Tprev, distortion_Tprev, normDistortion_Tprev] = colorCorrectionOptimize( ...
        useLambdaOrTargetInfo, lamArg, tgtInfoArg, tgtDistArg, ...
        LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        infoNormalizer, distortionNormalizer, Disp, 'T_init', T_prev);

    % Compare losses (second arg is the same scalar you optimized on)
    loss_TI    = lossFunction(useLambdaOrTargetInfo, thisX, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
    loss_Tprev = lossFunction(useLambdaOrTargetInfo, thisX, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

    % Select better of the two
    if loss_TI <= loss_Tprev
        LMSDaltonizedCalFormatSweep{i}     = LMS_TI;
        rgbLinDaltonizedCalFormatSweep{i}  = rgbLin_TI;
        infoNormalized{i}          = normInfo_TI;
        distortionNormalized{i}    = normDistortion_TI;
        transformRGBmatrixSweep{i} = T_TI;
        T_prev = T_TI;  % Update T_prev based on previous step

        achievedInfoVal = double(normInfo_TI);
        achievedDistVal = double(normDistortion_TI);

        targetInfoVsAchievedInfo{i} = [achievedInfoVal, thisX];


    else
        LMSDaltonizedCalFormatSweep{i}     = LMS_Tprev;
        rgbLinDaltonizedCalFormatSweep{i}  = rgbLin_Tprev;
        infoNormalized{i}          = normInfo_Tprev;
        distortionNormalized{i}    = normDistortion_Tprev;
        transformRGBmatrixSweep{i} = T_Tprev;
        T_prev = T_Tprev;  % Update T_prev based on previous step

        achievedInfoVal = double(normInfo_Tprev);
        achievedDistVal = double(normDistortion_Tprev);

        targetInfoVsAchievedInfo{i} = [achievedInfoVal, thisX];

    end

    [LMSDaltonizedRenderedCalFormatSweep{i},rgbLinDaltonizedRenderedCalFormatSweep{i},~] = DichromRenderLinear(LMSDaltonizedCalFormatSweep{i},dichromatType,Disp);


    % In lambda mode, collect the achieved info values
    if strcmpi(sweepAxis,'lambda')
        achievedInfoVec(i) = achievedInfoVal;
        achievedDistVec(i) = achievedDistVal;
    end
end

% Just call the achieved info values target so that downstream "info vs. distortion" 
% plots work unchanged
if strcmpi(sweepAxis,'lambda')
    targetInfoNormalized = achievedInfoVec;  
end


outputs = cell(nSteps, 1);  % NSTEPS × 1 CELL ARRAY

for i = 1:nSteps
    outputs{i} = struct( ...
        'LMSDaltonizedCalFormat',            LMSDaltonizedCalFormatSweep{i}, ...
        'rgbLinDaltonizedCalFormat',         rgbLinDaltonizedCalFormatSweep{i}, ...
        'LMSDaltonizedRenderedCalFormat',    LMSDaltonizedRenderedCalFormatSweep{i}, ...
        'rgbLinDaltonizedRenderedCalFormat', rgbLinDaltonizedRenderedCalFormatSweep{i}, ...
        'transformRGBmatrix',                transformRGBmatrixSweep{i}, ...
        'targetInfoNormalized',              targetInfoNormalized(i), ...
        'infoNormalized',                    infoNormalized{i}, ...
        'distortionNormalized',              distortionNormalized{i},...
        'imgParams',                         imgParams,...
        'Disp',                              Disp,...
        'useLambdaOrTargetInfo',             char(useLambdaOrTargetInfo),... 
        'targetInfoVsAchievedInfo',          targetInfoVsAchievedInfo{i});
end

% Save the outputs
saveTransformedOutputs(outputs, pathName, nSteps, obj.infoFcn, obj.infoParams, obj.distortionFcn, Disp);

end
