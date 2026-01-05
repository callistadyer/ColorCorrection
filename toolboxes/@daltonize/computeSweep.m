function [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
    LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
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

p = inputParser;
p.addParameter('rerunStep', [], @(x) isempty(x) || (isscalar(x) && x>=1));
p.parse(varargin{:});
rerunStep = p.Results.rerunStep;

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

saveSubdir = fullfile(saveBase, pathName, metricFolder, runFolder);
saveFile   = fullfile(saveSubdir, 'sweepOutputs.mat');

% new
if ~exist(saveSubdir, 'dir')
    mkdir(saveSubdir);
end

% Load/resume/rerun logic
outputs = {};
startStep = 1;
stopAfterThisStep = false;

% Helper to detect which steps are "done"
% if the struct exists and contains the 'transformRGBmatrix', and if that field is not empty
isDoneStep = @(s) isstruct(s) && isfield(s,'transformRGBmatrix') && ~isempty(s.transformRGBmatrix);

if exist(saveFile,'file')
    fprintf('computeSweep: Found existing sweep file:\n  %s\n', saveFile);
    S = load(saveFile);

    % If there is a thing called outputs and there's somethin in it...
    if isfield(S,'outputs') && ~isempty(S.outputs)
        outputs = S.outputs;

        % Initialize a array (one per sweep step) marking all steps as not done
        doneMask = false(1, numel(outputs));
        % Check inside each step - is there something in each? 
        for ii = 1:numel(outputs)
            doneMask(ii) = isDoneStep(outputs{ii});
        end

        % Identify the most recent sweep step marked as completed so we know
        % where to resume the sweep
        lastDone = find(doneMask, 1, 'last');

        % True if all nSteps have completed successfully
        allDone  = (numel(outputs) >= nSteps) && all(doneMask(1:nSteps));

        % Type A: rerun of a specific step
        if ~isempty(rerunStep)
            startStep = rerunStep;
            stopAfterThisStep = true;

            % Ensure outputs has at least nSteps entries by padding with empty cells
            if numel(outputs) < nSteps
                outputs{nSteps} = [];
            end

            % Just clear out the desired step so we can re-do
            outputs{startStep} = [];

            fprintf('computeSweep: RERUN mode -> step %d only (will overwrite step %d in sweepOutputs.mat)\n', ...
                startStep, startStep);

            % Type B: no rerunStep, but run is already complete -> load and return 
        elseif allDone
            fprintf('computeSweep: All %d steps already complete. Loading and returning.\n', nSteps);

            LMSDaltonizedCalFormatSweep            = cell(1,nSteps);
            rgbLinDaltonizedCalFormatSweep         = cell(1,nSteps);
            LMSDaltonizedRenderedCalFormatSweep    = cell(1,nSteps);
            rgbLinDaltonizedRenderedCalFormatSweep = cell(1,nSteps);
            transformRGBmatrixSweep                = cell(1,nSteps);

            targetInfoNormalized                   = nan(1,nSteps);
            targetDistortionNormalized             = nan(1,nSteps);
            infoNormalizedAchievedSweep            = nan(1,nSteps);
            distortionNormalizedAchievedSweep      = nan(1,nSteps);

            for i = 1:nSteps
                LMSDaltonizedCalFormatSweep{i}            = outputs{i}.LMSDaltonizedCalFormat;
                rgbLinDaltonizedCalFormatSweep{i}         = outputs{i}.rgbLinDaltonizedCalFormat;
                LMSDaltonizedRenderedCalFormatSweep{i}    = outputs{i}.LMSDaltonizedRenderedCalFormat;
                rgbLinDaltonizedRenderedCalFormatSweep{i} = outputs{i}.rgbLinDaltonizedRenderedCalFormat;
                transformRGBmatrixSweep{i}                = outputs{i}.transformRGBmatrix;

                if isfield(outputs{i},'targetInfoNormalized');           targetInfoNormalized(i) = outputs{i}.targetInfoNormalized; end
                if isfield(outputs{i},'targetDistortionNormalized');     targetDistortionNormalized(i) = outputs{i}.targetDistortionNormalized; end
                if isfield(outputs{i},'infoNormalizedAchievedSweep');    infoNormalizedAchievedSweep(i) = outputs{i}.infoNormalizedAchievedSweep; end
                if isfield(outputs{i},'distortionNormalizedAchievedSweep'); distortionNormalizedAchievedSweep(i) = outputs{i}.distortionNormalizedAchievedSweep; end
            end
            return;

            % Type C: partial run -> resume from next step 
        else
            % Maybe nothing is done...
            if isempty(lastDone)
                startStep = 1;
                fprintf('  No completed steps detected. Starting from step 1.\n');
            else % start right after the last one finished
                startStep = lastDone + 1;
                fprintf('  Completed steps: 1–%d. Resuming at step %d.\n', lastDone, startStep);
            end

            if numel(outputs) < nSteps
                outputs{nSteps} = [];
            end
        end

    else
        % File exists but doesn't contain outputs
        outputs = cell(1,nSteps);
        startStep = 1;
        fprintf('  sweepOutputs.mat exists but no usable outputs found. Starting fresh.\n');
    end

else
    % No prior file means start fresh
    outputs = cell(1,nSteps);
    startStep = 1;
    fprintf('computeSweep: No existing sweep file. Starting fresh.\n');
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
% LMSCalFormat_new = calFormatLMS_deut;
LMSContrastCalFormat_new = (LMSCalFormat_new - Disp.grayLMS) ./ Disp.grayLMS;

% Normalize distortion and info
normalizerValueToGetRawValue = 1;
% Info takes in the contrast image
infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
% Distortion takes in the regular LMS image (not contrast)
distortionNormalizer = obj.distortionFcn(LMSCalFormat, LMSCalFormat_new, imgParams, normalizerValueToGetRawValue, Disp, obj.distortionParams);
% Uncomment below if you change it so that distortion takes in contrast image
% distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, normalizerValueToGetRawValue, obj.distortionParams);

% Get info values for lambda = 0 and lambda = 1
[~,~,~,info_0,infoNormalized_0,distortion0, distortionNormalized0] = colorCorrectionOptimize("lambda", 0, [], [], LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, obj.infoParams, obj.distortionParams, infoNormalizer, distortionNormalizer, Disp);

[~,~,~,info_1,infoNormalized_1,distortion1, distortionNormalized1] = colorCorrectionOptimize("lambda", 1, [], [], LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, obj.infoParams, obj.distortionParams, infoNormalizer, distortionNormalizer, Disp);

% Get target info values interpolated between lambdas 0 and 1
targetInfoNormalized       = linspace(infoNormalized_0, infoNormalized_1, nSteps);
targetDistortionNormalized = linspace(distortionNormalized0, distortionNormalized1, nSteps);

% Try starting with the biggest distortion
% targetDistortionNormalized = fliplr(targetDistortionNormalized);

% Preallocate outputs
LMSDaltonizedCalFormatSweep              = cell(1, nSteps);
rgbLinDaltonizedCalFormatSweep           = cell(1, nSteps);
LMSDaltonizedRenderedCalFormatSweep      = cell(1, nSteps);
rgbLinDaltonizedRenderedCalFormatSweep   = cell(1, nSteps);
transformRGBmatrixSweep                  = cell(1, nSteps);
infoNormalizedAchievedSweep              = nan(1, nSteps);
distortionNormalizedAchievedSweep        = nan(1, nSteps);


% Initialize outputs cell
if numel(outputs) < nSteps
    outputs{nSteps} = [];
end

% If resuming, populate arrays from all of the prev steps outputs{1:startStep-1}
% for ii = 1:startStep-1
%     if isempty(outputs{ii}); continue; end
% 
%     LMSDaltonizedCalFormatSweep{ii}            = outputs{ii}.LMSDaltonizedCalFormat;
%     rgbLinDaltonizedCalFormatSweep{ii}         = outputs{ii}.rgbLinDaltonizedCalFormat;
%     LMSDaltonizedRenderedCalFormatSweep{ii}    = outputs{ii}.LMSDaltonizedRenderedCalFormat;
%     rgbLinDaltonizedRenderedCalFormatSweep{ii} = outputs{ii}.rgbLinDaltonizedRenderedCalFormat;
%     transformRGBmatrixSweep{ii}                = outputs{ii}.transformRGBmatrix;
% 
%     infoNormalizedAchievedSweep(ii)            = outputs{ii}.infoNormalizedAchievedSweep;
%     distortionNormalizedAchievedSweep(ii)      = outputs{ii}.distortionNormalizedAchievedSweep;
% end
% Populate arrays from already saved steps 
for ii = 1:nSteps
    if ii > numel(outputs) || isempty(outputs{ii})
        continue; 
    end
    if ~isDoneStep(outputs{ii})
        continue;
    end

    LMSDaltonizedCalFormatSweep{ii}            = outputs{ii}.LMSDaltonizedCalFormat;
    rgbLinDaltonizedCalFormatSweep{ii}         = outputs{ii}.rgbLinDaltonizedCalFormat;
    LMSDaltonizedRenderedCalFormatSweep{ii}    = outputs{ii}.LMSDaltonizedRenderedCalFormat;
    rgbLinDaltonizedRenderedCalFormatSweep{ii} = outputs{ii}.rgbLinDaltonizedRenderedCalFormat;
    transformRGBmatrixSweep{ii}                = outputs{ii}.transformRGBmatrix;

    if isfield(outputs{ii}, 'infoNormalizedAchievedSweep')
        infoNormalizedAchievedSweep(ii) = outputs{ii}.infoNormalizedAchievedSweep;
    end
    if isfield(outputs{ii}, 'distortionNormalizedAchievedSweep')
        distortionNormalizedAchievedSweep(ii) = outputs{ii}.distortionNormalizedAchievedSweep;
    end

    if isfield(outputs{ii}, 'targetInfoNormalized')
        targetInfoNormalized(ii) = outputs{ii}.targetInfoNormalized;
    end
    if isfield(outputs{ii}, 'targetDistortionNormalized')
        targetDistortionNormalized(ii) = outputs{ii}.targetDistortionNormalized;
    end
end


% Sweep through target vals
switch lower(sweepAxis)
    case 'info'
        xVec = targetInfoNormalized;
    case 'lambda'
        xVec = linspace(0,1,nSteps);
    case 'distortion'
        xVec = targetDistortionNormalized;
end


% If doing a distortion sweep, try to load and use the final T_prev solution
% from the info sweep
T_infoSeeds = [];
if strcmpi(sweepAxis,'distortion')

    infoRunFolder  = sprintf('info_%dsteps', nSteps);
    infoSaveFile   = fullfile(saveBase, pathName, metricFolder, infoRunFolder, 'sweepOutputs.mat');

    if exist(infoSaveFile,'file')
        S = load(infoSaveFile,'outputs');
        outputsInfo = S.outputs;   % nSteps x 1 cell

        % Extract all transform matrices
        T_infoSeeds = cell(nSteps,1);
        for i = 1:nSteps
            T_infoSeeds{i} = outputsInfo{i}.transformRGBmatrix;
        end

        fprintf('[seed] Loaded %d info-sweep T matrices\n', nSteps);
    else
        fprintf('[seed] No info sweep found at:\n  %s\n', infoSaveFile);
    end
end

% Starting points for search over matrices
T_I = eye(3);

if startStep > 1 && ~isempty(transformRGBmatrixSweep{startStep-1})
    T_prev = transformRGBmatrixSweep{startStep-1};
else
    T_prev = eye(3);
end

violationStep = false(1, nSteps);
% for i = 1:nSteps
if stopAfterThisStep
    steps = startStep;          % rerun mode: do ONLY this step
else
    steps = startStep:nSteps;   % normal mode: resume to end
end

for i = steps

    thisX = xVec(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find a good starting point for the distortion sweep (unconstrained on info)
    % The distortion sweep seems to have a harder time finding
    % solutions at each step that all together form a smooth curve. To try
    % and address this, we have been trying "better" starting points for
    % the distortion sweep.
    startPt = 'infoSoln';
    if strcmpi(sweepAxis,'distortion')

        % 1) Try starting the search at the transformation matrix found at that
        % step of the info sweep (that is, grab the solution for the info
        % sweep, and start there)
        if strcmp(startPt,'infoSoln')
            
            T_prev = T_infoSeeds{i};

            % 2) Try starting the search at a transformation matrix that
            % succesfully finds the target distortion (with no attention paid
            % to info)
            % Minimize (distortionNorm - target)^2 to get a feasible T,
            % then use that T as the starting point for the real constrained maximize-info step.
        elseif strcmp(startPt,'findDesiredDist')

            % Define the loss function as essentially (distortionNorm - target)^2
            feas_fun  = @(t_vec) lossFunction('distortion', thisX, t_vec, ...
                LMSCalFormat, imgParams, dichromatType, ...
                obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp,obj.infoParams);

            % Optimization options
            feas_opts = optimoptions('fmincon', 'Algorithm','interior-point','Display','none','MaxIterations',60,'ConstraintTolerance',1e-10,'StepTolerance',1e-10);

            % Gamut constraints
            [A_total, b_total,~] = buildGamutConstraints(LMSCalFormat, dichromatType, Disp);

            % Minimize
            T_feas_vec = fmincon(feas_fun, T_prev(:), A_total, b_total, [], [], [], [], [], feas_opts);

            % Check to see if we achieved the target distortion when there is
            % no info constraint
            [~, ~, ~, ~, distNorm_feas(i)] = lossFunction('lambda', 0.0, T_feas_vec, LMSCalFormat, imgParams, dichromatType, ...
                obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp,obj.infoParams);

            fprintf('[pre-solve] step %2d: targetDist=%.6g  achieved=%.6g  =%.3g\n', i, thisX, distNorm_feas(i));

            % Start the search at that transformation matrix
            T_prev = reshape(T_feas_vec,3,3);
        elseif strcmpi(startPt,'prevStep')
            % T_prev = transformRGBmatrixSweep{i-1};
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Depending on the sweep axis, the "target" value is different
    if strcmpi(sweepAxis,'info')
        lamArg = []; tgtInfoArg = thisX; tgtDistArg = [];
    elseif strcmpi(sweepAxis,'distortion')
        lamArg = []; tgtInfoArg = []; tgtDistArg = thisX;
    else % 'lambda'
        lamArg = thisX; tgtInfoArg = []; tgtDistArg = [];
    end

    % Optimize from identity starting point or at feasible distortion start
    [LMS_TI, rgbLin_TI, T_TI, info_TI, normInfo_TI, distortion_TI, normDistortion_TI, nlconSatisfied_TI] = colorCorrectionOptimize( ...
        sweepAxis, lamArg, tgtInfoArg, tgtDistArg, ...
        LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        obj.infoParams, obj.distortionParams,...
        infoNormalizer, distortionNormalizer, Disp, 'T_init', T_I);

    % Optimize from T_prev starting point
    [LMS_Tprev, rgbLin_Tprev, T_Tprev, info_Tprev, normInfo_Tprev, distortion_Tprev, normDistortion_Tprev,nlconSatisfied_Tprev] = colorCorrectionOptimize( ...
        sweepAxis, lamArg, tgtInfoArg, tgtDistArg, ...
        LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        obj.infoParams, obj.distortionParams,...
        infoNormalizer, distortionNormalizer, Disp, 'T_init', T_prev);

    % Compare losses
    switch lower(sweepAxis)
        case 'info'
            % Objective was "minimize distortion" -> lambda = 0
            loss_TI    = lossFunction('lambda', 0.0, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
            loss_Tprev = lossFunction('lambda', 0.0, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

        case 'distortion'
            % Objective was "maximize info" -> lambda = 1
            loss_TI    = lossFunction('lambda', 1.0, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
            loss_Tprev = lossFunction('lambda', 1.0, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

        case 'lambda'
            % Objective is the mixed lambda objective
            loss_TI    = lossFunction('lambda', thisX, T_TI(:),    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
            loss_Tprev = lossFunction('lambda', thisX, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
    end

    useTI = true;

    switch lower(sweepAxis)
        case 'distortion'
            viol_TI    = abs(normDistortion_TI    - thisX);
            viol_Tprev = abs(normDistortion_Tprev - thisX);
        case 'info'
            viol_TI    = abs(normInfo_TI    - thisX);
            viol_Tprev = abs(normInfo_Tprev - thisX);
        otherwise
            viol_TI    = NaN;
            viol_Tprev = NaN;
    end

    % If nonlinear constraint not satisfied from identity starting point -> use Tprev
    if ~nlconSatisfied_TI && nlconSatisfied_Tprev
        useTI = false;
        % If nonlinear constraint not satisfied from T_prev starting point -> use identity start
    elseif ~nlconSatisfied_Tprev && nlconSatisfied_TI
        useTI = true;
        % If nonlinear constraint satisfied from both starting points -> use whichever has lower loss
    elseif nlconSatisfied_Tprev && nlconSatisfied_TI
        %%%%%%% NOT SURE IF THIS IS A GOOD DECISION?
        useTI = (loss_TI <= loss_Tprev);
        % useTI = false;
        % useTI = true;
        % If both starting points fail to satisfy the constraint...
    elseif ~nlconSatisfied_Tprev && ~nlconSatisfied_TI
        % error('Nonlinear constraints never satisfied')
        fprintf(['\n[WARNING] Nonlinear constraint NOT satisfied\n' ...
            '  sweep step        : %d\n' ...
            '  sweep axis        : %s\n' ...
            '  target distortion : %.6g\n' ...
            '  T_I    -> achieved distortion = %.6g | |Δ| = %.3g\n' ...
            '  T_prev -> achieved distortion = %.6g | |Δ| = %.3g\n'], ...
            i, lower(char(sweepAxis)), thisX, ...
            normDistortion_TI,    viol_TI, ...
            normDistortion_Tprev, viol_Tprev);

        useTI = (viol_TI <= viol_Tprev);          % neither feasible -> use smaller violation???
    end

    % Select better of the two losses
    if useTI

        LMSDaltonizedCalFormatSweep{i}       = LMS_TI;
        rgbLinDaltonizedCalFormatSweep{i}    = rgbLin_TI;
        infoNormalizedAchievedSweep(i)       = normInfo_TI;
        distortionNormalizedAchievedSweep(i) = normDistortion_TI;
        transformRGBmatrixSweep{i}           = T_TI;
        % Update T_prev based on previous step
        T_prev = T_TI;

        targetInfoVsAchievedInfo{i} = [normInfo_TI, thisX];
        targetDistVsAchievedDist{i} = [normDistortion_TI, thisX];


    else

        LMSDaltonizedCalFormatSweep{i}        = LMS_Tprev;
        rgbLinDaltonizedCalFormatSweep{i}     = rgbLin_Tprev;
        infoNormalizedAchievedSweep(i)        = normInfo_Tprev;
        distortionNormalizedAchievedSweep(i)  = normDistortion_Tprev;
        transformRGBmatrixSweep{i}            = T_Tprev;
        % Update T_prev based on previous step
        T_prev = T_Tprev;
    end

    % Violation check: evaluate the constraint band on the
    % chosen transform
    pctTol      = 0.01;   % 1% of target
    absTolFloor = 1e-3;   % minimum absolute tolerance
    eps = max(absTolFloor, pctTol * thisX);
 
    eps = 1e-3;

    T_chosen = transformRGBmatrixSweep{i};

    switch lower(sweepAxis)
        case 'info'
            [~, ~, infoNorm_chk, ~, distNorm_chk] = lossFunction('lambda', 0.0, T_chosen(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

            % Do you get within the range for the info constraint band? 
            violationStep(i) = (round(abs(infoNorm_chk - thisX),3) > (eps + 1e-12));

        case 'distortion'
            [~, ~, infoNorm_chk, ~, distNorm_chk] = lossFunction('lambda', 0.0, T_chosen(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

            % Do you get within the range for the distortion constraint band? 
            violationStep(i) = (round(abs(distNorm_chk - thisX),3) > (eps + 1e-12));

        otherwise
            violationStep(i) = false;
    end

    % If violation, let's take a look
    if violationStep(i)

        % Identify which initialization the SELECTED solution used
        if useTI
            chosenInit = 'identity start (T_I)';
        else
            chosenInit = 'previous start (T_prev)';
        end

        % Get the relevent info
        switch lower(sweepAxis)
            case 'info'
                targetVal   = thisX;
                achievedVal = infoNorm_chk;
                deltaVal    = abs(achievedVal - targetVal);
                metricName  = 'infoNorm';
            case 'distortion'
                targetVal   = thisX;
                achievedVal = distNorm_chk;
                deltaVal    = abs(achievedVal - targetVal);
                metricName  = 'distNorm';
        end

        fprintf('\n================= CONSTRAINT VIOLATION =================\n');
        fprintf('step=%d | sweep=%s\n', i, lower(char(sweepAxis)));
        fprintf('chosen init: %s\n', chosenInit);
        fprintf('target %s = %.6g\n', metricName, targetVal);
        fprintf('achieved %s = %.6g\n', metricName, achievedVal);
        fprintf('delta = %.6g', deltaVal);

        % Also show the two candidate solutions
        if strcmpi(sweepAxis,'info')
            fprintf('candidate (T_I)    achieved info=%.6g | dist=%.6g\n', normInfo_TI,   normDistortion_TI);
            fprintf('candidate (T_prev) achieved info=%.6g | dist=%.6g\n', normInfo_Tprev, normDistortion_Tprev);
        elseif strcmpi(sweepAxis,'distortion')
            fprintf('candidate (T_I)    achieved dist=%.6g | info=%.6g\n', normDistortion_TI,   normInfo_TI);
            fprintf('candidate (T_prev) achieved dist=%.6g | info=%.6g\n', normDistortion_Tprev, normInfo_Tprev);
        end

        % keyboard;  
        % pause;   
    end
 

    fprintf('step %2d | sweep=%s | target=%.4f | achievedInfo=%.4f | achievedDist=%.4f\n', ...
        i, lower(char(sweepAxis)), thisX, infoNormalizedAchievedSweep(i), distortionNormalizedAchievedSweep(i));

    % Plot Achieved info vs distortion up through this step
    debugPlot_InfoVsDistSingleSweep(i, sweepAxis, distortionNormalizedAchievedSweep, infoNormalizedAchievedSweep, targetDistortionNormalized, targetInfoNormalized, violationStep, thisX, pathName);

    % Get the dichromat rendering
    [LMSDaltonizedRenderedCalFormatSweep{i},rgbLinDaltonizedRenderedCalFormatSweep{i},~] = DichromRenderLinear(LMSDaltonizedCalFormatSweep{i},dichromatType,Disp);

    % Save this step to disk immediately
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
    fprintf('  [checkpoint] saved step %d/%d -> %s\n', i, nSteps, saveFile);

    % Save a per-step snapshot so we can rerun any step instantly
    stepSnapFile = fullfile(saveSubdir, sprintf('STEP_SNAP_%s_%03d.mat', lower(char(sweepAxis)), i));
    save(stepSnapFile, 'obj', 'LMSCalFormat', 'imgParams', 'dichromatType','nSteps', 'pathName', 'sweepAxis', 'i', 'thisX', 'infoNormalizer', 'distortionNormalizer', 'targetInfoNormalized', 'targetDistortionNormalized', 'T_I', 'T_prev', '-v7.3');

    % Stop immediately after this step if requested
    if stopAfterThisStep
        fprintf('computeSweep: Finished rerun of step %d. Exiting.\n', i);
        break;
    end
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

% Helper function to draw the info–distortion curve after each optimization step
function debugPlot_InfoVsDistSingleSweep(i, sweepAxis, ...
    achDist, achInfo, ...
    targetDistortionNormalized, targetInfoNormalized, ...
    violationStep, ...
    thisX, runID)

% Declare persistent maps so each (runID, sweepAxis) pair keeps its own figure
persistent hFigMap hAxMap

% Initialize the persistent maps on first call
if isempty(hFigMap)
    hFigMap = containers.Map('KeyType','char','ValueType','any');
    hAxMap  = containers.Map('KeyType','char','ValueType','any');
end

% Build a unique key identifying this run and sweep type
key = sprintf('%s__%s', char(string(runID)), lower(char(sweepAxis)));

% If no figure exists for this key (or it was closed), create a new one
if ~isKey(hFigMap, key) || ~isvalid(hFigMap(key))
    hFig = figure('Name', sprintf('Info vs Distortion (debug) | %s | sweep=%s', char(string(runID)), lower(char(sweepAxis))), 'NumberTitle','off');

    % Create axes inside the new figure
    hAx  = axes('Parent', hFig);

    % Store the handles so future calls reuse them
    hFigMap(key) = hFig;
    hAxMap(key)  = hAx;
else
    % Reuse the existing axes for this run/sweep
    hAx = hAxMap(key);
end

infoColor = [0.55 0.20 0.80];
distColor = [0.10 0.70 0.85];
lineStyle = ':';
lineWidth = 0.8;

switch sweepAxis
    case 'info'
        sweepColor = infoColor;
    case 'distortion'
        sweepColor = distColor;
    otherwise
        sweepColor = [0.35 0.35 0.35];
end

cla(hAx);
hold(hAx,'on');
grid(hAx,'off');
axis(hAx,'square');

xlabel(hAx,'Distortion (normalized)');
ylabel(hAx,'Info (normalized)');

% Draw horizontal dotted lines for all target info values
if exist('yline','file') && ~isempty(targetInfoNormalized)
    for kk = 1:numel(targetInfoNormalized)
        h = yline(hAx, targetInfoNormalized(kk), lineStyle, 'Color', infoColor, 'LineWidth', lineWidth);
        h.HandleVisibility = 'off';
    end
end

% Draw vertical dotted lines for all target distortion values
if exist('xline','file') && ~isempty(targetDistortionNormalized)
    for kk = 1:numel(targetDistortionNormalized)
        h = xline(hAx, targetDistortionNormalized(kk), lineStyle, 'Color', distColor, 'LineWidth', lineWidth);
        h.HandleVisibility = 'off';
    end
end

% Find steps where both achieved distortion and info are defined
valid = ~isnan(achDist) & ~isnan(achInfo);
kAll = find(valid);

% Indices of steps up through and including the current step
kSoFar = kAll(kAll <= i);

% Plot the achieved info–distortion curve up to the current step
if ~isempty(kSoFar)
    plot(hAx, achDist(kSoFar), achInfo(kSoFar), '-', 'Color', sweepColor, 'LineWidth', 1.5);
end

% Create a dummy plot handle for achieved points (for legend)
hAch = plot(hAx, nan, nan, 'o', 'Color', sweepColor, 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'none');

% Plot achieved points as open circles
if ~isempty(kSoFar)
    plot(hAx, achDist(kSoFar), achInfo(kSoFar), 'o', 'Color', sweepColor, 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', 'none', 'HandleVisibility','off');
end

% Create a dummy plot handle for violation points
% hViol = plot(hAx, nan, nan, 'o', 'MarkerSize', 10, 'MarkerFaceColor', sweepColor,'MarkerEdgeColor', sweepColor);

% Find indices of steps marked as violations
% badIdx = kAll(violationStep(kAll));

% Plot violation points as filled circles
% if ~isempty(badIdx)
%     plot(hAx, achDist(badIdx), achInfo(badIdx), 'o', 'MarkerSize', 12, 'MarkerFaceColor', sweepColor, 'MarkerEdgeColor', sweepColor, 'HandleVisibility','off');
% end

% Title
switch sweepAxis
    case 'info'
        titleStr = sprintf([ ...
            'Info vs Distortion (single sweep) | info | up to step %d\n' ...
            'TargetInfo=%.6g | AchievedInfo=%.6g | delta=%.3g'], ...
            i, thisX, achInfo(i), abs(achInfo(i) - thisX));

    case 'distortion'
        titleStr = sprintf([ ...
            'Info vs Distortion (single sweep) | distortion | up to step %d\n' ...
            'TargetDist=%.6g | AchievedDist=%.6g | delta=%.3g'], ...
            i, thisX, achDist(i), abs(achDist(i) - thisX));
end
title(hAx, titleStr);

% Legend
% legend(hAx, [hAch hViol], {sprintf('Achieved (%s sweep)', sweepAxis), sprintf('Violation (%s sweep)', sweepAxis)}, 'Location','best');

drawnow;

end
