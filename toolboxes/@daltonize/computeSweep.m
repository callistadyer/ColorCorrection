function [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
    transformRGBmatrixSweep, targetInfoNormalized, targetDistortionNormalized, infoNormalizedAchievedSweep, distortionNormalizedAchievedSweep] = computeSweep(obj,...
    LMSCalFormat, imgParams, dichromatType, nSteps,pathName,sweepAxis,varargin)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load it %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectName   = 'ColorCorrection';
outputDir     = getpref(projectName, 'outputDir');
saveBase      = fullfile(outputDir, 'testImagesTransformed');
metricFolder  = buildMetricFolderName(obj.infoFcn, obj.infoParams, obj.distortionFcn);
runFolder     = sprintf('%s_%dsteps', lower(char(sweepAxis)), nSteps);
saveSubdir    = fullfile(saveBase, pathName, metricFolder, runFolder);
saveFile      = fullfile(saveSubdir, 'sweepOutputs.mat');
% Best folder (final chosen per-step winners)
bestSubdir    = fullfile(saveSubdir, 'best');
if ~exist(bestSubdir, 'dir'); mkdir(bestSubdir); end
bestSaveFile  = fullfile(bestSubdir, 'sweepOutputs_best.mat');   % separate from main
if ~exist(saveSubdir, 'dir'); mkdir(saveSubdir); end
forwardDir    = fullfile(saveSubdir, 'forward');
backwardDir   = fullfile(saveSubdir, 'backward');
bestDir       = fullfile(saveSubdir, 'best');

if ~exist(forwardDir,'dir');  mkdir(forwardDir);  end
if ~exist(backwardDir,'dir'); mkdir(backwardDir); end
if ~exist(bestDir,'dir');     mkdir(bestDir);     end

% Defaults
outputs   = cell(1, nSteps);
startStep = 1;

% Load function
[outputs, startStep, allDone, loaded] = loadSweepOutputs(saveFile, nSteps);

if allDone
    LMSDaltonizedCalFormatSweep            = loaded.LMSDaltonizedCalFormatSweep;
    rgbLinDaltonizedCalFormatSweep         = loaded.rgbLinDaltonizedCalFormatSweep;
    LMSDaltonizedRenderedCalFormatSweep    = loaded.LMSDaltonizedRenderedCalFormatSweep;
    rgbLinDaltonizedRenderedCalFormatSweep = loaded.rgbLinDaltonizedRenderedCalFormatSweep;
    transformRGBmatrixSweep                = loaded.transformRGBmatrixSweep;

    targetInfoNormalized                   = loaded.targetInfoNormalized;
    targetDistortionNormalized             = loaded.targetDistortionNormalized;
    infoNormalizedAchievedSweep            = loaded.infoNormalizedAchievedSweep;
    distortionNormalizedAchievedSweep      = loaded.distortionNormalizedAchievedSweep;
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If it doesn't already exist, go ahead and compute the transformation as normal:
% Get display
Disp = obj.Disp;
% Get the contrast LMS values
LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

% Get the dichromat renderings of the original, which is used for the normalizer
%%%% MAYBE JUST TAKE THE DEUTERANOPE VERSION TO COMPARE WITH THE ORIGINAL
[calFormatLMS_prot, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Protanopia', Disp);
[calFormatLMS_deut, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Deuteranopia', Disp);
[calFormatLMS_trit, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Tritanopia', Disp);
LMSCalFormat_new = [calFormatLMS_prot(1,:); calFormatLMS_deut(2,:); calFormatLMS_trit(3,:)];
% LMSCalFormat_new = [LMSCalFormat(1,:); calFormatLMS_deut(2,:); LMSCalFormat(3,:)];

% LMSCalFormat_new = calFormatLMS_deut;
LMSContrastCalFormat_new = (LMSCalFormat_new - Disp.grayLMS) ./ Disp.grayLMS;

% Normalize distortion and info
normalizerValueToGetRawValue = 1;
% Info takes in the contrast image
infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
infoNormalizer       = abs(infoNormalizer);
% Distortion takes in the regular LMS image (not contrast)
distortionNormalizer = obj.distortionFcn(LMSCalFormat, LMSCalFormat_new, imgParams, normalizerValueToGetRawValue, Disp, obj.distortionParams);

% Uncomment below if you change it so that distortion takes in contrast image
% distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, normalizerValueToGetRawValue, obj.distortionParams);

% Get info values for lambda = 0 and lambda = 1
% [~,~,~,info_0,infoNormalized_0,distortion0, distortionNormalized0] = colorCorrectionOptimize("lambda", 0, [], [], LMSCalFormat, imgParams, dichromatType, ...
%     obj.infoFcn, obj.distortionFcn, obj.infoParams, obj.distortionParams, infoNormalizer, distortionNormalizer, Disp);
% 
% [~,~,~,info_1,infoNormalized_1,distortion1, distortionNormalized1] = colorCorrectionOptimize("lambda", 1, [], [], LMSCalFormat, imgParams, dichromatType, ...
%     obj.infoFcn, obj.distortionFcn, obj.infoParams, obj.distortionParams, infoNormalizer, distortionNormalizer, Disp);

[infoNormalized_0, infoNormalized_1, distortionNormalized0, distortionNormalized1] = computeEndpoints(saveSubdir, LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp);

% Get target info values interpolated between lambdas 0 and 1
targetInfoNormalized       = linspace(infoNormalized_0, infoNormalized_1, nSteps);
targetDistortionNormalized = linspace(distortionNormalized0, distortionNormalized1, nSteps);

% Information range between the two endpoints
infoRange = abs(infoNormalized_1 - infoNormalized_0);
distRange = abs(distortionNormalized1 - distortionNormalized0);
% Allow 1% deviation relative to the information range
infoTolerance = 0.01 * infoRange;
distTolerance = 0.01 * distRange;
fprintf('Info range: %.6f\n', infoRange);
fprintf('Info tolerance: %.6f\n', infoTolerance);

% Preallocate outputs
LMSDaltonizedCalFormatSweep              = cell(1, nSteps);
rgbLinDaltonizedCalFormatSweep           = cell(1, nSteps);
LMSDaltonizedRenderedCalFormatSweep      = cell(1, nSteps);
rgbLinDaltonizedRenderedCalFormatSweep   = cell(1, nSteps);
transformRGBmatrixSweep                  = cell(1, nSteps);
infoNormalizedAchievedSweep              = nan(1, nSteps);
distortionNormalizedAchievedSweep        = nan(1, nSteps);

% Make sure outputs is at least nSteps long (safe even if already true)
if numel(outputs) < nSteps
    outputs{nSteps} = [];
end

% Only populate steps that are complete
isDoneStep = @(s) isstruct(s) ...
    && isfield(s,'stepCompleted') && isequal(s.stepCompleted,true) ...
    && isfield(s,'transformRGBmatrix') && ~isempty(s.transformRGBmatrix);

for ii = 1:nSteps
    % If there is no saved output at all, obviously can't load it
    if isempty(outputs{ii})
        continue;
    end

    % Just recompute this step if it never finished running
    if ~isDoneStep(outputs{ii})
        continue;
    end

    % Restore optimized outputs for this step
    LMSDaltonizedCalFormatSweep{ii}            = outputs{ii}.LMSDaltonizedCalFormat;
    rgbLinDaltonizedCalFormatSweep{ii}         = outputs{ii}.rgbLinDaltonizedCalFormat;
    LMSDaltonizedRenderedCalFormatSweep{ii}    = outputs{ii}.LMSDaltonizedRenderedCalFormat;
    rgbLinDaltonizedRenderedCalFormatSweep{ii} = outputs{ii}.rgbLinDaltonizedRenderedCalFormat;
    transformRGBmatrixSweep{ii}                = outputs{ii}.transformRGBmatrix;

    % Restore achieved metrics if they were saved 
    if isfield(outputs{ii}, 'infoNormalizedAchievedSweep')
        infoNormalizedAchievedSweep(ii) = outputs{ii}.infoNormalizedAchievedSweep;
    end
    if isfield(outputs{ii}, 'distortionNormalizedAchievedSweep')
        distortionNormalizedAchievedSweep(ii) = outputs{ii}.distortionNormalizedAchievedSweep;
    end
end


% Sweep through target vals
% sweepAxis = 'lambda'
switch lower(sweepAxis)
    case 'info'
        xVec = targetInfoNormalized;
    case 'lambda'
        xVec = linspace(0,1,nSteps);
    case 'distortion'
        xVec = targetDistortionNormalized;
end

% Compute candidate start points for the optimizer
modes = {'infoSoln','findDesiredDist','findDesiredInfo'};   % pick what you want
startPts = findStartPoints(modes, sweepAxis, nSteps, startStep, xVec, transformRGBmatrixSweep, saveBase, pathName, metricFolder, LMSCalFormat, imgParams, dichromatType, infoNormalizer, distortionNormalizer, Disp, obj, saveSubdir);

T_I    = startPts.T_I;
T_prev = startPts.T_prev0;
T_infoSeeds       = startPts.T_infoSeeds;
T_findDesiredDist = startPts.T_findDesiredDist;
T_findDesiredInfo = startPts.T_findDesiredInfo;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORWARD PASS (startStep -> nSteps)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%`
forwardBest = repmat(struct('name','','T',[],'LMS',[],'rgbLin',[], 'infoNorm',NaN,'distNorm',NaN,'loss',Inf,'nonLinSatisfied',false,'violation',Inf), 1, nSteps);

T_prev_fwd = T_prev;   % your existing startPts.T_prev0

for i = startStep:nSteps

    thisX = xVec(i);

    % Skip optimizing step 1 (we force identity at commit anyway)
    if i == 1
        % Don't worry, the info and distortion will be computed before save
        forwardBest(i) = struct('name','ForcedIdentity','T',eye(3), 'LMS',[],'rgbLin',[],'infoNorm',NaN,'distNorm',NaN,'loss',NaN, 'nonLinSatisfied',true,'violation',0);
        T_prev_fwd = eye(3);
        continue;
    end

    [bestF, ~] = chooseStartPoints(i, sweepAxis, thisX, T_I, T_prev_fwd, T_infoSeeds, T_findDesiredDist, T_findDesiredInfo, LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp, saveSubdir, 'forward',infoTolerance,distTolerance);

    forwardBest(i) = bestF;

    % Chain T_prev forward
    if ~isempty(bestF.T)
        T_prev_fwd = bestF.T;
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACKWARD PASS (nSteps-1 -> startStep)
% Initialize T_prev using the *forward* end-point
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
backwardBest = repmat(struct('name','','T',[],'LMS',[],'rgbLin',[], 'infoNorm',NaN,'distNorm',NaN,'loss',Inf,'nonLinSatisfied',false,'violation',Inf), 1, nSteps);

% When going backwards, don't start from identity.
% Start at the second-to-last point, using T_prev from the previous solution (the forward final point).
T_prev_bwd = forwardBest(nSteps).T;

for i = (nSteps-1):-1:startStep

    thisX = xVec(i);

    % Skip optimizing step 1 (we force identity at commit anyway)
    if i == 1
        % Don't worry, the info and distortion will be computed before save
        backwardBest(i) = struct('name','ForcedIdentity','T',eye(3), ...
            'LMS',[],'rgbLin',[],'infoNorm',NaN,'distNorm',NaN,'loss',NaN, ...
            'nonLinSatisfied',true,'violation',0);
        T_prev_bwd = eye(3);
        continue;
    end

    [bestB, startPtSolns] = chooseStartPoints(i, sweepAxis, thisX, T_I, T_prev_bwd, T_infoSeeds, T_findDesiredDist, T_findDesiredInfo, ...
        LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp,...
        saveSubdir, 'backward',infoTolerance,distTolerance);

    backwardBest(i) = bestB;

    % Chain T_prev backward
    if ~isempty(bestB.T)
        T_prev_bwd = bestB.T;
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose per-step winner between forward vs backward
% Rule: constraints first, then lower loss.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = startStep:nSteps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Force step 1 (least distortion endpoint) to be identity
    if i == 1
        T_I = eye(3);

        % Apply identity transform exactly (no optimizer)
        triRGBCalFormat         = Disp.M_cones2rgb * LMSCalFormat;
        triRGBContrastCalFormat = (triRGBCalFormat - Disp.grayRGB) ./ Disp.grayRGB;

        triRGBContrastCalFormat_T = T_I * triRGBContrastCalFormat;      % identity
        rgbLin_I = (triRGBContrastCalFormat_T .* Disp.grayRGB) + Disp.grayRGB;

        % Safety clip (identity should already be in gamut)
        % rgbLin_I(rgbLin_I > 1) = 1;
        % rgbLin_I(rgbLin_I < 0) = 0;

        LMS_I = Disp.M_rgb2cones * rgbLin_I;

        % Evaluate metrics for bookkeeping (same path as optimizer uses)
        [~, ~, infoNorm_I, ~, distNorm_I] = lossFunction('lambda', 0.0, T_I(:), ...
            LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
            infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

        best = struct( ...
            'name',            'FORCED_IDENTITY', ...
            'pass',            'forced', ...
            'T',               T_I, ...
            'LMS',             LMS_I, ...
            'rgbLin',          rgbLin_I, ...
            'infoNorm',        infoNorm_I, ...
            'distNorm',        distNorm_I, ...
            'loss',            NaN, ...
            'nonLinSatisfied', true, ...
            'violation',       0 );

        % Commit forced step 1
        LMSDaltonizedCalFormatSweep{i}       = best.LMS;
        rgbLinDaltonizedCalFormatSweep{i}    = best.rgbLin;
        transformRGBmatrixSweep{i}           = best.T;
        infoNormalizedAchievedSweep(i)       = best.infoNorm;
        distortionNormalizedAchievedSweep(i) = best.distNorm;

        fprintf('step %2d | winner=%-12s | forced identity\n', i, best.name);

        [LMSDaltonizedRenderedCalFormatSweep{i}, rgbLinDaltonizedRenderedCalFormatSweep{i}, ~] = ...
            DichromRenderLinear(LMSDaltonizedCalFormatSweep{i}, dichromatType, Disp);

        outputs{i} = struct( ...
            'LMSDaltonizedCalFormat',            LMSDaltonizedCalFormatSweep{i}, ...
            'rgbLinDaltonizedCalFormat',         rgbLinDaltonizedCalFormatSweep{i}, ...
            'LMSDaltonizedRenderedCalFormat',    LMSDaltonizedRenderedCalFormatSweep{i}, ...
            'rgbLinDaltonizedRenderedCalFormat', rgbLinDaltonizedRenderedCalFormatSweep{i}, ...
            'transformRGBmatrix',                transformRGBmatrixSweep{i}, ...
            'targetInfoNormalized',              targetInfoNormalized(i), ...
            'targetDistortionNormalized',        targetDistortionNormalized(i), ...
            'infoNormalizedAchievedSweep',       infoNormalizedAchievedSweep(i), ...
            'distortionNormalizedAchievedSweep', distortionNormalizedAchievedSweep(i), ...
            'imgParams',                         imgParams, ...
            'Disp',                              Disp, ...
            'sweepAxis',                         char(sweepAxis), ...
            'stepCompleted',                     true, ...
            'timestamp',                         datestr(now));

        save(saveFile, 'outputs', '-v7.3');

        bestStepDir = fullfile(bestSubdir, sprintf('step_%03d', i));
        if ~exist(bestStepDir,'dir'); mkdir(bestStepDir); end
        bestStepFile = fullfile(bestStepDir, 'best.mat');
        thisX = xVec(i);  
        save(bestStepFile, 'best', 'thisX', 'sweepAxis', 'i', '-v7.3');
        fprintf('  [best] saved -> %s\n', bestStepFile);

        continue;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If it's not the endpoint, go ahead and get the best...
    thisX = xVec(i);

    bestF = forwardBest(i);
    bestB = backwardBest(i);

    % If backward wasn't computed for i (e.g., i==nSteps), treat it as missing
    hasB = ~isempty(bestB.T);

    if ~hasB
        best = bestF;
        best.pass = 'forward';
    else
        best = pickBetterSolution(bestF, bestB);
    end

    % Commit chosen solution for this step
    LMSDaltonizedCalFormatSweep{i}       = best.LMS;
    rgbLinDaltonizedCalFormatSweep{i}    = best.rgbLin;
    transformRGBmatrixSweep{i}           = best.T;
    infoNormalizedAchievedSweep(i)       = best.infoNorm;
    distortionNormalizedAchievedSweep(i) = best.distNorm;

    fprintf('step %2d | winner=%-12s | from=%-8s | nonLinOK=%d | loss=%.4g | viol=%.4g\n', i, best.name, best.pass, best.nonLinSatisfied, best.loss, best.violation);

    [LMSDaltonizedRenderedCalFormatSweep{i}, rgbLinDaltonizedRenderedCalFormatSweep{i}, ~] = DichromRenderLinear(LMSDaltonizedCalFormatSweep{i}, dichromatType, Disp);

    % DOES THIS WORK?
    % [LMSDaltonizedRenderedCalFormat{i}, rgbLinDaltonizedRenderedCalFormat{i}, ~] = obj.renderFcn(LMSDaltonizedCalFormat{i}, dichromatType, Disp);

    outputs{i} = struct( ...
        'LMSDaltonizedCalFormat',            LMSDaltonizedCalFormatSweep{i}, ...
        'rgbLinDaltonizedCalFormat',         rgbLinDaltonizedCalFormatSweep{i}, ...
        'LMSDaltonizedRenderedCalFormat',    LMSDaltonizedRenderedCalFormatSweep{i}, ...
        'rgbLinDaltonizedRenderedCalFormat', rgbLinDaltonizedRenderedCalFormatSweep{i}, ...
        'transformRGBmatrix',                transformRGBmatrixSweep{i}, ...
        'targetInfoNormalized',              targetInfoNormalized(i), ...
        'targetDistortionNormalized',        targetDistortionNormalized(i), ...
        'infoNormalizedAchievedSweep',       infoNormalizedAchievedSweep(i), ...
        'distortionNormalizedAchievedSweep', distortionNormalizedAchievedSweep(i), ...
        'imgParams',                         imgParams, ...
        'Disp',                              Disp, ...
        'sweepAxis',                         char(sweepAxis), ...
        'stepCompleted',                     true, ...
        'timestamp',                         datestr(now));

    save(saveFile, 'outputs', '-v7.3');
    % Best checkpoint in best/ folder 
    % save(bestSaveFile, 'outputs', '-v7.3');
    % fprintf('  [best] saved step %d/%d -> %s\n', i, nSteps, bestSaveFile);
    bestStepDir = fullfile(bestSubdir, sprintf('step_%03d', i));
    if ~exist(bestStepDir,'dir'); mkdir(bestStepDir); end
    bestStepFile = fullfile(bestStepDir, 'best.mat');
    save(bestStepFile, 'best', 'thisX', 'sweepAxis', 'i', '-v7.3');
    fprintf('  [best] saved -> %s\n', bestStepFile);
end



% figure(); plot(distNorm_feas,xVec,'-o')
% disp([distNorm_feas(:), xVec(:)]);

% for i = 1:nSteps
%     rgbLin = rgbLinDaltonizedRenderedCalFormatSweep{i};
%     rgb = rgbLin2RGB(rgbLin,Disp);
%     rgbImg = CalFormatToImage(rgb,imgParams.m,imgParams.n);
%     figure();
%     imagesc(rgbImg)
%     axis square
% end

% Save the outputs
saveTransformedOutputs(outputs, pathName, nSteps, obj.infoFcn, obj.infoParams, obj.distortionFcn, Disp, 'sweepAxis', sweepAxis);

end

% Choose between two candidate solutions (nested)
function best = pickBetterSolution(bestF, bestB)
    % Add pass labels for logging
    bestF.pass = 'forward';
    bestB.pass = 'backward';
    
    % 1) If only one satisfies nonlinear constraints -> take it
    if bestF.nonLinSatisfied && ~bestB.nonLinSatisfied
        best = bestF; return;
    elseif bestB.nonLinSatisfied && ~bestF.nonLinSatisfied
        best = bestB; return;
    end
    
    % 2) If both satisfy constraints -> choose lower loss
    if bestF.nonLinSatisfied && bestB.nonLinSatisfied
        if bestF.loss <= bestB.loss
            best = bestF;
        else
            best = bestB;
        end
        return;
    end
    
    % 3) If neither satisfies constraints -> smaller violation wins, then lower loss
    if bestF.violation < bestB.violation
        best = bestF;
    elseif bestB.violation < bestF.violation
        best = bestB;
    else
        % tie on violation: choose lower loss
        if bestF.loss <= bestB.loss
            best = bestF;
        else
            best = bestB;
        end
    end
end