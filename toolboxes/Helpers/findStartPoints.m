function startPts = findStartPoints(modes, sweepAxis, nSteps, startStep, xVec, transformRGBmatrixSweep, saveBase, pathName, metricFolder, LMSCalFormat, imgParams, dichromatType, infoNormalizer, distortionNormalizer, Disp, obj, saveSubdir)
% function startPts = findStartPoints(modes, sweepAxis, nSteps, startStep, xVec, transformRGBmatrixSweep, saveBase, pathName, metricFolder, LMSCalFormat, imgParams, dichromatType, infoNormalizer, distortionNormalizer, Disp, obj)
% 
% Description:
%   This helper builds a set of candidate initialization transforms (T_init)
%   that computeSweep can use when optimizing each sweep step. The goal is
%   to (1) resume smoothly from previous results and/or (2) supply better
%   starting points for difficult constrained problems (especially the
%   distortion sweep).
%
%   Supported starting-point strategies are selected via 'modes':
%     'infoSoln'         - load per-step transform matrices from the INFO sweep
%     'findDesiredDist'  - presolve for each step to achieve a target distortion
%     'findDesiredInfo'  - presolve for each step to achieve a target info
%
% Inputs:
%   modes                  - Cell array or string array of requested strategies.
%                            Example: {'infoSoln','findDesiredDist'}.
%   sweepAxis              - Sweep type currently running in computeSweep:
%                            'info', 'distortion', or 'lambda'.
%   nSteps                 - Total number of sweep steps.
%   startStep              - First step that still needs to be computed (resume index).
%   xVec                   - 1 x nSteps vector of target values for this sweep:
%                            if sweepAxis='info'       -> targetInfoNormalized
%                            if sweepAxis='distortion' -> targetDistortionNormalized
%                            if sweepAxis='lambda'     -> lambda values (presolve modes are skipped).
%   transformRGBmatrixSweep - 1 x nSteps cell array of already-loaded transforms
%                             (may contain empties for unfinished steps).
%   saveBase               - Base directory containing sweep output folders.
%                            Example: fullfile(outputDir,'testImagesTransformed').
%   pathName               - Image/run identifier used to build paths.
%                            Example: 'Deuteranopia/flower2.png/s1_m32_n32'.
%   metricFolder           - Metric folder name (built from info/distortion choices).
%                            Example: buildMetricFolderName(...).
%   LMSCalFormat           - 3 x N LMS values of the input image (cal format).
%   imgParams              - Struct with image-related parameters (m, n, etc.).
%   dichromatType          - Dichromat type string passed to render/constraints:
%                            'Protanopia', 'Deuteranopia', or 'Tritanopia'.
%   infoNormalizer         - Scalar normalizer used by lossFunction/info metric
%                            (same value used by computeSweep/colorCorrectionOptimize).
%   distortionNormalizer   - Scalar normalizer used by lossFunction/distortion metric
%                            (same value used by computeSweep/colorCorrectionOptimize).
%   Disp                   - Display/calibration struct (e.g., contains grayLMS, matrices).
%   obj                    - Daltonizer object (provides obj.infoFcn, obj.distortionFcn,
%                            and obj.infoParams needed by lossFunction).
%
% Outputs (returned in struct startPts):
%   startPts.T_I            - 3 x 3 identity matrix (always available baseline).
%   startPts.T_prev0        - 3 x 3 previous-step transform used for resuming:
%                             transformRGBmatrixSweep{startStep-1} if available, else eye(3).
%   startPts.T_infoSeeds    - nSteps x 1 cell array of transforms loaded from an INFO sweep
%                             (only populated if 'infoSoln' requested and file exists).
%   startPts.T_findDesiredDist - nSteps x 1 cell array of presolved transforms that try to
%                             achieve target distortion values (only if 'findDesiredDist'
%                             requested and sweepAxis='distortion').
%   startPts.T_findDesiredInfo - nSteps x 1 cell array of presolved transforms that try to
%                             achieve target info values (only if 'findDesiredInfo'
%                             requested and sweepAxis='info').
%   startPts.meta           - Struct containing bookkeeping fields:
%                             .modes, .sweepAxis (useful for debugging/logging).
%
% Notes:
%   - Presolve modes ('findDesiredDist' and 'findDesiredInfo') run fmincon once
%     per step, so they can be expensive.

%%%% Normalize and interpret requested modes

if ischar(modes) || isstring(modes)
    modes = cellstr(modes);
end

modes = lower(string(modes(:))');

wantInfoSoln        = any(modes == "infosoln");
wantFindDesiredDist = any(modes == "finddesireddist");
wantFindDesiredInfo = any(modes == "finddesiredinfo");

%%%% Initialize output struct and always-available starting points

startPts = struct();

% Record how these start points were constructed (useful for debugging)
startPts.meta.modes     = cellstr(modes);
startPts.meta.sweepAxis = char(sweepAxis);

% Identity matrix is always a valid baseline initialization
startPts.T_I = eye(3);

% Default previous-step initialization is identity
startPts.T_prev0 = eye(3);

% If resuming and a previous transform exists, use it as T_prev0
if startStep > 1 ...
        && numel(transformRGBmatrixSweep) >= (startStep-1) ...
        && ~isempty(transformRGBmatrixSweep{startStep-1})
    startPts.T_prev0 = transformRGBmatrixSweep{startStep-1};
end

% Initialize optional fields so they always exist
startPts.T_infoSeeds       = [];
startPts.T_findDesiredDist = [];
startPts.T_findDesiredInfo = [];

%%%% infoSoln mode: load transforms from a completed info sweep

if wantInfoSoln && strcmpi(sweepAxis,'distortion')

    % Build path to the info sweep outputs
    infoRunFolder = sprintf('info_%dsteps', nSteps);
    infoSaveFile  = fullfile(saveBase, pathName, metricFolder, infoRunFolder, 'sweepOutputs.mat');

    if exist(infoSaveFile, 'file')

        % Load only the outputs variable
        S = load(infoSaveFile, 'outputs');

        if isfield(S,'outputs') && ~isempty(S.outputs)

            outputsInfo = S.outputs;
            T_infoSeeds = cell(nSteps, 1);

            % Extract transform matrices step-by-step
            for i = 1:nSteps
                if i <= numel(outputsInfo) ...
                        && ~isempty(outputsInfo{i}) ...
                        && isfield(outputsInfo{i}, 'transformRGBmatrix')
                    T_infoSeeds{i} = outputsInfo{i}.transformRGBmatrix;
                else
                    T_infoSeeds{i} = [];
                end
            end

            startPts.T_infoSeeds = T_infoSeeds;
            fprintf('Loaded %d info-sweep T matrices\n', nSteps);

        else
            fprintf('Info sweep file exists but contains no usable outputs:\n  %s\n', infoSaveFile);
        end
    else
        fprintf('No info sweep found at:\n  %s\n', infoSaveFile);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This stuff takes longer. Save and reload if you can
presolveRoot = fullfile(saveSubdir, 'presolve');    % lives alongside sweepOutputs.mat
if ~exist(presolveRoot, 'dir'); mkdir(presolveRoot); end

distCacheDir = fullfile(presolveRoot, 'findDesiredDist');
infoCacheDir = fullfile(presolveRoot, 'findDesiredInfo');

if ~exist(distCacheDir, 'dir'); mkdir(distCacheDir); end
if ~exist(infoCacheDir, 'dir'); mkdir(infoCacheDir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build shared optimization components for presolve strategies

% Linear gamut constraints are reused across all presolve steps
[A_total, b_total, ~] = buildGamutConstraints(LMSCalFormat, dichromatType, Disp);

% Optimization options for presolve (kept relatively lightweight)
feas_opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','none', ...
    'MaxIterations',60, ...
    'ConstraintTolerance',1e-10, ...
    'StepTolerance',1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findDesiredDist mode: presolve to hit target distortion values

if wantFindDesiredDist && strcmpi(sweepAxis,'distortion')

    T_findDesiredDist = cell(nSteps, 1);

    % Chain solutions so each presolve starts from the previous one
    T_seed = startPts.T_prev0;

    for i = startStep:nSteps

        targetDist = xVec(i);

        % Cache file for this presolve step
        stepFile = fullfile(distCacheDir, sprintf('step_%03d.mat', i));

        % Try load
        loadedOK = false;
        if exist(stepFile, 'file')
            S = load(stepFile);
            if isfield(S,'completed') && isequal(S.completed,true) && ...
                    isfield(S,'T_feas_vec') && ~isempty(S.T_feas_vec)
                T_feas_vec = S.T_feas_vec;
                loadedOK = true;
                fprintf('[findDesiredDist] LOADED step %2d from cache: %s\n', i, stepFile);
            end
        end

        if ~loadedOK
            % Objective: squared error between achieved and target distortion
            feas_fun = @(t_vec) localDistObjective(t_vec, targetDist,LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp);

            % Solve constrained presolve optimization
            T_feas_vec = fmincon(feas_fun, T_seed(:), A_total, b_total, [], [], [], [], [], feas_opts);

            % Save immediately so we never rerun unless cache is deleted
            completed = true;
            timestamp = datestr(now);
            save(stepFile, 'T_feas_vec', 'targetDist', 'completed', 'timestamp', '-v7.3');
            fprintf('[findDesiredDist] SAVED  step %2d -> %s\n', i, stepFile);
        end

        % Evaluate achieved distortion for diagnostics (works for loaded or computed)
        [~, ~, ~, ~, distNorm] = lossFunction( ...
            'lambda', 0.0, T_feas_vec, ...
            LMSCalFormat, imgParams, dichromatType, ...
            obj.infoFcn, obj.distortionFcn, ...
            infoNormalizer, distortionNormalizer, ...
            Disp, obj.infoParams);

        fprintf('[findDesiredDist] step %2d: target=%.6g achieved=%.6g delta=%.3g\n', ...
            i, targetDist, distNorm, abs(distNorm - targetDist));

        % Store result and chain seed
        T_findDesiredDist{i} = reshape(T_feas_vec, 3, 3);
        T_seed = T_findDesiredDist{i};
    end

    startPts.T_findDesiredDist = T_findDesiredDist;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findDesiredInfo mode: presolve to hit target info values

if wantFindDesiredInfo && strcmpi(sweepAxis,'info')

    T_findDesiredInfo = cell(nSteps, 1);

    % Chain solutions for stability across steps
    T_seed = startPts.T_prev0;

    for i = startStep:nSteps

        targetInfo = xVec(i);

        % Cache file for this presolve step
        stepFile = fullfile(infoCacheDir, sprintf('step_%03d.mat', i));

        % Try load
        loadedOK = false;
        if exist(stepFile, 'file')
            S = load(stepFile);
            if isfield(S,'completed') && isequal(S.completed,true) && ...
                    isfield(S,'T_feas_vec') && ~isempty(S.T_feas_vec)
                T_feas_vec = S.T_feas_vec;
                loadedOK = true;
                fprintf('[findDesiredInfo] LOADED step %2d from cache: %s\n', i, stepFile);
            end
        end

        if ~loadedOK
            feas_fun = @(t_vec) localInfoObjective(t_vec, targetInfo, LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp);

            T_feas_vec = fmincon(feas_fun, T_seed(:), A_total, b_total, [], [], [], [], [], feas_opts);

            completed = true;
            timestamp = datestr(now);
            save(stepFile, 'T_feas_vec', 'targetInfo', 'completed', 'timestamp', '-v7.3');
            fprintf('[findDesiredInfo] SAVED  step %2d -> %s\n', i, stepFile);
        end

        % Diagnostics
        [~, ~, infoNorm, ~, ~] = lossFunction( ...
            'lambda', 0.0, T_feas_vec, ...
            LMSCalFormat, imgParams, dichromatType, ...
            obj.infoFcn, obj.distortionFcn, ...
            infoNormalizer, distortionNormalizer, ...
            Disp, obj.infoParams);


        fprintf('[findDesiredInfo] step %2d: target=%.6g achieved=%.6g delta=%.3g\n', ...
            i, targetInfo, infoNorm, abs(infoNorm - targetInfo));

        T_findDesiredInfo{i} = reshape(T_feas_vec, 3, 3);
        T_seed = T_findDesiredInfo{i};
    end

    startPts.T_findDesiredInfo = T_findDesiredInfo;
end



end

function f = localDistObjective(t_vec, targetDist, LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp)

% Squared error between achieved and target distortion

[~, ~, ~, ~, distNorm] = lossFunction('lambda', 0.0, t_vec, LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

d = distNorm - targetDist;
f = d * d;
end


function f = localInfoObjective(t_vec, targetInfo, LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp)

% Squared error between achieved and target info

[~, ~, infoNorm, ~, ~] = lossFunction('lambda', 0.0, t_vec, LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

d = infoNorm - targetInfo;
f = d * d;
end
