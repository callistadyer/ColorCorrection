function [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
    LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
    transformRGBmatrixSweep, targetInfoNormalized, targetDistortionNormalized, infoNormalizedAchievedSweep, distortionNormalizedAchievedSweep] = computeSweep(obj,...
    LMSCalFormat, imgParams, dichromatType, nSteps,pathName,sweepAxis)
% computeSweep  Sweep through target info or target distortion or lambda values and optimize color correction.
%
% Syntax:
%   [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
%    LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
%    transformRGBmatrixSweep, targetInfoNormalized, targetDistortionNormalized, infoNormalized, distortionNormalized] = computeSweep(obj,...
%    LMSCalFormat, imgParams, dichromatType, nSteps,pathName,sweepAxis)
%
% Description:
%   This method sweeps between the info values produced by lambda=0 and lambda=1,
%   then optimizes the transformation for each intermediate targetInfo or
%   targetDistortion
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
%

if isempty(nSteps)
    nSteps = 10;
end

if ~exist('sweepAxis','var') || isempty(sweepAxis)
    sweepAxis = 'info';
end

%%%%%%%%%%%%%%%%%%%%%%% Load it %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check to see if this output already exists, and then load it
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');
saveBase    = fullfile(outputDir, 'testImagesTransformed');

% Build the metric folder
metricFolder = buildMetricFolderName(obj.infoFcn, obj.infoParams, obj.distortionFcn);

% New run folder that specifies whether the sweep is info or distortion
runFolder = sprintf('%s_%dsteps', lower(char(sweepAxis)), nSteps);
% runFolder = sprintf('%dsteps', nSteps);

saveSubdir = fullfile(saveBase, pathName, metricFolder, runFolder);
saveFile   = fullfile(saveSubdir, 'sweepOutputs.mat');

if exist(saveFile, 'file')
    fprintf('[computeSweep] Loading old sweep from: %s\n', saveFile);
    loaded  = load(saveFile);
    outputs = loaded.outputs;

    nSteps = numel(outputs);
    LMSDaltonizedCalFormatSweep            = cell(1,nSteps);
    rgbLinDaltonizedCalFormatSweep         = cell(1,nSteps);
    LMSDaltonizedRenderedCalFormatSweep    = cell(1,nSteps);
    rgbLinDaltonizedRenderedCalFormatSweep = cell(1,nSteps);
    transformRGBmatrixSweep                = cell(1,nSteps);

    targetInfoNormalized                   = zeros(1,nSteps);
    targetDistortionNormalized             = zeros(1,nSteps);

    infoNormalizedAchievedSweep            = zeros(1,nSteps);
    distortionNormalizedAchievedSweep      = zeros(1,nSteps);

    for i = 1:nSteps
        LMSDaltonizedCalFormatSweep{i}            = outputs{i}.LMSDaltonizedCalFormat;
        rgbLinDaltonizedCalFormatSweep{i}         = outputs{i}.rgbLinDaltonizedCalFormat;
        LMSDaltonizedRenderedCalFormatSweep{i}    = outputs{i}.LMSDaltonizedRenderedCalFormat;
        rgbLinDaltonizedRenderedCalFormatSweep{i} = outputs{i}.rgbLinDaltonizedRenderedCalFormat;
        transformRGBmatrixSweep{i}                = outputs{i}.transformRGBmatrix;

        targetInfoNormalized(i)                   = outputs{i}.targetInfoNormalized;
        targetDistortionNormalized(i)             = outputs{i}.targetDistortionNormalized;

        infoNormalizedAchievedSweep(i)            = outputs{i}.infoNormalizedAchievedSweep;
        distortionNormalizedAchievedSweep(i)      = outputs{i}.distortionNormalizedAchievedSweep;
    end
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
distortionNormalizer = obj.distortionFcn(LMSCalFormat, LMSCalFormat_new, imgParams, normalizerValueToGetRawValue, obj.distortionParams);

% Get info values for lambda = 0 and lambda = 1
[~,~,~,info_0,infoNormalized_0,distortion0, distortionNormalized0] = colorCorrectionOptimize("lambda", 0, [], [], LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

[~,~,~,info_1,infoNormalized_1,distortion1, distortionNormalized1] = colorCorrectionOptimize("lambda", 1, [], [], LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

% Get target info values interpolated between lambdas 0 and 1
targetInfoNormalized       = linspace(infoNormalized_0, infoNormalized_1, nSteps);
targetDistortionNormalized = linspace(distortionNormalized0, distortionNormalized1, nSteps);
% start with the biggest distortion
% targetDistortionNormalized = fliplr(targetDistortionNormalized);

% Preallocate outputs
LMSDaltonizedCalFormatSweep              = cell(1, nSteps);
rgbLinDaltonizedCalFormatSweep           = cell(1, nSteps);
LMSDaltonizedRenderedCalFormatSweep      = cell(1, nSteps);
rgbLinDaltonizedRenderedCalFormatSweep   = cell(1, nSteps);
transformRGBmatrixSweep                  = cell(1, nSteps);
infoNormalizedAchievedSweep              = zeros(1, nSteps);
distortionNormalizedAchievedSweep        = zeros(1, nSteps);

%% Sweep through target vals

switch lower(sweepAxis)
    case 'info'
        xVec = targetInfoNormalized;                
    case 'lambda'
        xVec = linspace(0,1,nSteps);               
    case 'distortion'
        xVec = targetDistortionNormalized;               
end

T_I = eye(3);
T_prev = eye(3);

for i = 1:nSteps
    thisX = xVec(i);

    
    % Find a good starting point for the distortion sweep (unconstrained on info) 
    % Minimize (distortionNorm - target)^2 to get a feasible T,
    % then use that T as the starting point for the real constrained maximize-info step.
    if strcmpi(sweepAxis,'distortion')
        feas_fun  = @(t_vec) lossFunction('distortion', thisX, t_vec, ...
                            LMSCalFormat, imgParams, dichromatType, ...
                            obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
        feas_opts = optimoptions('fmincon', 'Algorithm','interior-point', ...
                                 'Display','none','MaxIterations',60, ...
                                 'ConstraintTolerance',1e-10,'StepTolerance',1e-10);

        % Constraints
        [A_total, b_total,~] = buildGamutConstraints(LMSCalFormat, dichromatType, Disp);

        T_feas_vec = fmincon(feas_fun, T_prev(:), A_total, b_total, [], [], [], [], [], feas_opts);

        % Check to see if we achieved the target distortion when there is
        % no info constraint
        [~, ~, ~, ~, distNorm_feas(i)] = lossFunction('lambda', 0.0, T_feas_vec, ...
            LMSCalFormat, imgParams, dichromatType, ...
            obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
        fprintf('[pre-solve] step %2d: targetDist=%.6g  achieved=%.6g  =%.3g\n', ...
            i, thisX, distNorm_feas(i));

        T_prev = reshape(T_feas_vec,3,3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if strcmpi(sweepAxis,'info')
        % info sweep: lambda=[], targetInfo=thisX
        lamArg = []; tgtInfoArg = thisX; tgtDistArg = [];
    elseif strcmpi(sweepAxis,'distortion')
        lamArg = []; tgtInfoArg = []; tgtDistArg = thisX;
    else % 'lambda'
        % lambda sweep: lambda=thisX, targetInfo=[]
        lamArg = thisX; tgtInfoArg = []; tgtDistArg = [];
    end

    % Optimize from identity starting point or at feasible distortion start
    [LMS_TI, rgbLin_TI, T_TI, info_TI, normInfo_TI, distortion_TI, normDistortion_TI] = colorCorrectionOptimize( ...
        sweepAxis, lamArg, tgtInfoArg, tgtDistArg, ...
        LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        infoNormalizer, distortionNormalizer, Disp, 'T_init', T_I);

    % Optimize from T_prev starting point
    [LMS_Tprev, rgbLin_Tprev, T_Tprev, info_Tprev, normInfo_Tprev, distortion_Tprev, normDistortion_Tprev] = colorCorrectionOptimize( ...
        sweepAxis, lamArg, tgtInfoArg, tgtDistArg, ...
        LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        infoNormalizer, distortionNormalizer, Disp, 'T_init', T_prev);

    % Compare losses (second arg is the same scalar you optimized on)
    % loss_TI    = lossFunction(sweepAxis, thisX, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
    % loss_Tprev = lossFunction(sweepAxis, thisX, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
    
    switch lower(sweepAxis)
        case 'info'
            % Objective was "minimize distortion" -> lambda = 0
            loss_TI    = lossFunction('lambda', 0.0, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
            loss_Tprev = lossFunction('lambda', 0.0, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

        case 'distortion'
            % Objective was "maximize info" -> lambda = 1
            loss_TI    = lossFunction('lambda', 1.0, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
            loss_Tprev = lossFunction('lambda', 1.0, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

        case 'lambda'
            % Objective is the mixed lambda objective
            loss_TI    = lossFunction('lambda', thisX, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
            loss_Tprev = lossFunction('lambda', thisX, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
    end

    % Select better of the two
    if loss_TI <= loss_Tprev

        LMSDaltonizedCalFormatSweep{i}       = LMS_TI;
        rgbLinDaltonizedCalFormatSweep{i}    = rgbLin_TI;
        infoNormalizedAchievedSweep(i)       = normInfo_TI;
        distortionNormalizedAchievedSweep(i) = normDistortion_TI;
        transformRGBmatrixSweep{i}           = T_TI;
        % Update T_prev based on previous step
        T_prev = T_TI;  

        targetInfoVsAchievedInfo{i} = [normInfo_TI, thisX];
        targetDistVsAchievedDist{i} = [normInfo_TI, thisX];


    else

        LMSDaltonizedCalFormatSweep{i}        = LMS_Tprev;
        rgbLinDaltonizedCalFormatSweep{i}     = rgbLin_Tprev;
        infoNormalizedAchievedSweep(i)        = normInfo_Tprev;
        distortionNormalizedAchievedSweep(i)  = normDistortion_Tprev;
        transformRGBmatrixSweep{i}            = T_Tprev;
        % Update T_prev based on previous step
        T_prev = T_Tprev; 
    end

    [LMSDaltonizedRenderedCalFormatSweep{i},rgbLinDaltonizedRenderedCalFormatSweep{i},~] = DichromRenderLinear(LMSDaltonizedCalFormatSweep{i},dichromatType,Disp);

end
% figure(); plot(distNorm_feas,xVec,'-o')
% disp([distNorm_feas(:), xVec(:)]);


outputs = cell(nSteps, 1);  

for i = 1:nSteps
    outputs{i} = struct( ...
        'LMSDaltonizedCalFormat',            LMSDaltonizedCalFormatSweep{i}, ...
        'rgbLinDaltonizedCalFormat',         rgbLinDaltonizedCalFormatSweep{i}, ...
        'LMSDaltonizedRenderedCalFormat',    LMSDaltonizedRenderedCalFormatSweep{i}, ...
        'rgbLinDaltonizedRenderedCalFormat', rgbLinDaltonizedRenderedCalFormatSweep{i}, ...
        'transformRGBmatrix',                transformRGBmatrixSweep{i}, ...
        'targetInfoNormalized',              targetInfoNormalized(i), ...
        'targetDistortionNormalized',        targetDistortionNormalized(i), ...
        'infoNormalizedAchievedSweep',       infoNormalizedAchievedSweep(i), ...
        'distortionNormalizedAchievedSweep', distortionNormalizedAchievedSweep(i),...
        'imgParams',                         imgParams,...
        'Disp',                              Disp,...
        'sweepAxis',                         char(sweepAxis));
end

% Save the outputs
% saveTransformedOutputs(outputs, pathName, nSteps, obj.infoFcn, obj.infoParams, obj.distortionFcn, Disp);
saveTransformedOutputs(outputs, pathName, nSteps, obj.infoFcn, obj.infoParams, obj.distortionFcn, Disp, 'sweepAxis', sweepAxis);
end
