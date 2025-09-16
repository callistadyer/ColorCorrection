function [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
    LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
    transformRGBmatrixSweep, targetInfoNormalized, infoNormalized, distortionNormalized] = computeInfoSweep(obj,...
    LMSCalFormat, imgParams, dichromatType, nSteps,pathName)
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
%   LMSDaltonizedCalFormatSweep            - Cell array of optimized LMS images (this is
%                                            the size of 1xnSteps)
%   rgbLinDaltonizedCalFormatSweep         - Cell array of optimized RGB images
%   LMSDaltonizedRenderedCalFormatSweep    - Cell array of
%                                            LMSDaltonizedCalFormatSweep rendered for a dichromat
%   rgbLinDaltonizedRenderedCalFormatSweep - Cell array of
%                                            rgbLinDaltonizedCalFormatSweep rendered for a dichromat
%   transformRGBmatrixSweep - Cell array of 3x3 transformation matrices
%   targetInfoNormalized    - Vector of target info values used
%   infoNormalized          - Cell array of normalized info values
%   distortionNormalized    - Cell array of normalized distortion values
%

if isempty(nSteps)
    nSteps = 10;
end

% Check to see if this output already exists, and then load it
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');
saveBase    = fullfile(outputDir, 'testImagesTransformed');

% NEW: build the metric folder EXACTLY like saveTransformedOutputs
metricFolder = buildMetricFolderName(obj.infoFcn, obj.infoParams, obj.distortionFcn);

runFolder = sprintf('%dsteps', nSteps);
saveSubdir = fullfile(saveBase, pathName, metricFolder, runFolder);
saveFile   = fullfile(saveSubdir, 'sweepOutputs.mat');

% --- Optional: backward-compatibility to older folder names (pre-regress change) ---
if ~exist(saveFile, 'file')
    % Old naming variants you used before (adjust if you had others)
    old1_metricFolder = sprintf('%s_%s', func2str(obj.infoFcn), func2str(obj.distortionFcn));   % single underscore
    old2_metricFolder = sprintf('%s__%s', func2str(obj.infoFcn), func2str(obj.distortionFcn));  % double underscore

    old1_saveFile = fullfile(saveBase, pathName, old1_metricFolder, runFolder, 'sweepOutputs.mat');
    old2_saveFile = fullfile(saveBase, pathName, old2_metricFolder, runFolder, 'sweepOutputs.mat');

    if exist(old1_saveFile, 'file')
        saveFile = old1_saveFile;
    elseif exist(old2_saveFile, 'file')
        saveFile = old2_saveFile;
    end
end

if exist(saveFile, 'file')
    fprintf('[computeInfoSweep] Loading cached sweep from: %s\n', saveFile);
    loaded  = load(saveFile);
    outputs = loaded.outputs;

    % Extract individual fields back to output variables
    nSteps = numel(outputs);
    LMSDaltonizedCalFormatSweep            = cell(1,nSteps);
    rgbLinDaltonizedCalFormatSweep         = cell(1,nSteps);
    LMSDaltonizedRenderedCalFormatSweep    = cell(1,nSteps);
    rgbLinDaltonizedRenderedCalFormatSweep = cell(1,nSteps);
    transformRGBmatrixSweep                = cell(1,nSteps);
    targetInfoNormalized                   = zeros(1,nSteps);
    infoNormalized                         = cell(1,nSteps);
    distortionNormalized                   = cell(1,nSteps);

    for i = 1:nSteps
        LMSDaltonizedCalFormatSweep{i}            = outputs{i}.LMSDaltonizedCalFormat;
        rgbLinDaltonizedCalFormatSweep{i}         = outputs{i}.rgbLinDaltonizedCalFormat;
        LMSDaltonizedRenderedCalFormatSweep{i}    = outputs{i}.LMSDaltonizedRenderedCalFormat;
        rgbLinDaltonizedRenderedCalFormatSweep{i} = outputs{i}.rgbLinDaltonizedRenderedCalFormat;
        transformRGBmatrixSweep{i}                = outputs{i}.transformRGBmatrix;
        targetInfoNormalized(i)                   = outputs{i}.targetInfoNormalized;
        infoNormalized{i}                         = outputs{i}.infoNormalized;
        distortionNormalized{i}                   = outputs{i}.distortionNormalized;
    end
    return;
end

% Otherwise, go ahead and compute the transformation as normal
Disp = obj.Disp;
LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

[calFormatLMS_prot, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Protanopia', Disp);
[calFormatLMS_deut, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Deuteranopia', Disp);
[calFormatLMS_trit, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Tritanopia', Disp);
LMSCalFormat_new = [calFormatLMS_prot(1,:); calFormatLMS_deut(2,:); calFormatLMS_trit(3,:)];
LMSContrastCalFormat_new = (LMSCalFormat_new - Disp.grayLMS) ./ Disp.grayLMS;

% Normalize distortion and info
normalizerValueToGetRawValue = 1;

infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, normalizerValueToGetRawValue, obj.distortionParams);

% Get info values for lambda = 0 and lambda = 1
[~,~,~,info_0,infoNormalized_0,distortion0, distortionNormalized0] = colorCorrectionOptimize("lambda", 0, LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

[~,~,~,info_1,infoNormalized_1,distortion1, distortionNormalized1] = colorCorrectionOptimize("lambda", 1, LMSCalFormat, imgParams, dichromatType, ...
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

%% Sweep through target infos
T_I = eye(3);
T_prev = eye(3);
for i = 1:nSteps
    thisTargetInfo = targetInfoNormalized(i);

    % Optimize from identity starting point
    [LMS_TI, rgbLin_TI, T_TI, info_TI,normInfo_TI,distortion_TI, normDistortion_TI] = colorCorrectionOptimize( ...
        "targetInfo", thisTargetInfo, LMSCalFormat, imgParams, ...
        dichromatType, obj.infoFcn, obj.distortionFcn, ...
        infoNormalizer, distortionNormalizer,...
        Disp,'T_init',T_I);

    % Optimize from T_prev starting point
    [LMS_Tprev, rgbLin_Tprev, T_Tprev, info_Tprev, normInfo_Tprev,distortion_Tprev, normDistortion_Tprev] = colorCorrectionOptimize( ...
        "targetInfo", thisTargetInfo, LMSCalFormat, imgParams, ...
        dichromatType, obj.infoFcn, obj.distortionFcn, ...
        infoNormalizer, distortionNormalizer,...
        Disp,'T_init',T_prev);

    % Compare losses
    loss_TI    = lossFunction("targetInfo", thisTargetInfo, T_TI(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
    loss_Tprev = lossFunction("targetInfo", thisTargetInfo, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

    % Select better of the two
    if loss_TI <= loss_Tprev
        LMSDaltonizedCalFormatSweep{i}     = LMS_TI;
        rgbLinDaltonizedCalFormatSweep{i}  = rgbLin_TI;
        infoNormalized{i}          = normInfo_TI;
        distortionNormalized{i}    = normDistortion_TI;
        transformRGBmatrixSweep{i} = T_TI;
        T_prev = T_TI;  % Update T_prev based on previous step
    else
        LMSDaltonizedCalFormatSweep{i}     = LMS_Tprev;
        rgbLinDaltonizedCalFormatSweep{i}  = rgbLin_Tprev;
        infoNormalized{i}          = normInfo_Tprev;
        distortionNormalized{i}    = normDistortion_Tprev;
        transformRGBmatrixSweep{i} = T_Tprev;
        T_prev = T_Tprev;  % Update T_prev based on previous step
    end

    [LMSDaltonizedRenderedCalFormatSweep{i},rgbLinDaltonizedRenderedCalFormatSweep{i},~] = DichromRenderLinear(LMSDaltonizedCalFormatSweep{i},dichromatType,Disp);
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
        'Disp',                              Disp);
end

% Save the outputs
saveTransformedOutputs(outputs, pathName, nSteps, obj.infoFcn, obj.infoParams, obj.distortionFcn, Disp);
% saveTransformedOutputs(outputs, pathName, nSteps, obj.infoFcn, obj.distortionFcn, Disp);

end

function metricFolder = buildMetricFolderName(infoFcn, infoParams, distortionFcn)
    infoFcnName       = func2str(infoFcn);
    distortionFcnName = func2str(distortionFcn);

    if strcmp(infoFcnName, 'computeInfo_regress')
        if ~isfield(infoParams,'predictingWhat') || ~isfield(infoParams,'predictingFromWhat')
            error('infoParams must include predictingWhat and predictingFromWhat for computeInfo_regress.');
        end
        paramsStrRaw = sprintf('%s-from-%s', infoParams.predictingWhat, infoParams.predictingFromWhat);
        paramsSlug   = regexprep(paramsStrRaw, '[^A-Za-z0-9]+', '_');  % keep alnum + underscores
        paramsSlug   = regexprep(paramsSlug, '^_+|_+$', '');           % trim leading/trailing "_"
        metricFolder = sprintf('%s__%s__%s', infoFcnName, paramsSlug, distortionFcnName);
    else
        metricFolder = sprintf('%s__%s', infoFcnName, distortionFcnName);
    end
end
