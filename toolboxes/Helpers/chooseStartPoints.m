function [best, startPtSolns, candNames, candTinit] = chooseStartPoints(i, sweepAxis, thisX, T_I, T_prev, T_infoSeeds, T_findDesiredDist, T_findDesiredInfo, ...
    LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp, saveSubdir, passName)
% chooseStartPoints  Build a list of candidate start points and run the optimizer from each, then pick the best.
%
% Syntax:
%   [best, startPtSolns] = chooseStartPoints(i, sweepAxis, thisX, ...
%       T_I, T_prev, T_infoSeeds, T_findDesiredDist, T_findDesiredInfo, ...
%       LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp, saveSubdir, passName)
%
% Description:
%   This takes all possible start-point sources (identity, previous step, info-sweep seed,
%   and any feasibility-seeking starts) and (1) builds the candidate list for this step i, (2) runs
%   colorCorrectionOptimize from each candidate, and (3) returns the best feasible solution according
%   to the sweep's objective (loss) while respecting constraints.
%
% Inputs:
%   i                  - Current sweep step index (scalar int)
%   sweepAxis          - 'info', 'distortion', or 'lambda' (char/string)
%   thisX              - Current target value for this step (scalar)
%   T_I                - Identity start matrix (3x3)
%   T_prev             - Previous-step solution to chain from (3x3)
%   T_infoSeeds        - Cell array of 3x3 matrices from info sweep (nSteps x 1) or [] if none
%   T_findDesiredDist  - Cell array of 3x3 matrices that hit desired distortion (nSteps x 1) or [] if none
%   T_findDesiredInfo  - Cell array of 3x3 matrices that hit desired info (nSteps x 1) or [] if none
%   LMSCalFormat       - 3 x N LMS values of the input image
%   imgParams          - Struct with image-related params
%   dichromatType      - 'Protanopia' | 'Deuteranopia' | 'Tritanopia'
%   obj                - Daltonizer object holding info/distortion functions and params
%   infoNormalizer     - Scalar normalizer returned by obj.infoFcn(...)
%   distortionNormalizer - Scalar normalizer returned by obj.distortionFcn(...)
%   Disp               - Display struct
%
% Outputs:
%   best         - Struct for the chosen winner (fields: name, T, LMS, rgbLin, infoNorm, distNorm, loss, nonLinSatisfied, violation)
%   startPtSolns - Struct array with one entry per candidate start point (same fields as best)
%   candNames    - Cell array of candidate names that were tested (1 x K)
%   candTinit    - Cell array of candidate 3x3 start matrices (1 x K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% CHOOSE START POINT OPTIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build candidate name list
candNames = {};
candTinit = {};

% Always try identity (baseline)
candNames{end+1} = 'T_I';
candTinit{end+1} = T_I;

% Always try chaining from previous step
candNames{end+1} = 'T_prev';
candTinit{end+1} = T_prev;

% Try info-sweep solution at this step 
if exist('T_infoSeeds','var') && ~isempty(T_infoSeeds) && numel(T_infoSeeds) >= i && ~isempty(T_infoSeeds{i})
    candNames{end+1} = 'infoSoln';
    candTinit{end+1} = T_infoSeeds{i};
end

% Try feasibility-start that hits desired distortion 
if exist('T_findDesiredDist','var') && ~isempty(T_findDesiredDist) && numel(T_findDesiredDist) >= i && ~isempty(T_findDesiredDist{i})
    candNames{end+1} = 'findDesiredDist';
    candTinit{end+1} = T_findDesiredDist{i};
end

% Try feasibility-start that hits desired info
if exist('T_findDesiredInfo','var') && ~isempty(T_findDesiredInfo) && numel(T_findDesiredInfo) >= i && ~isempty(T_findDesiredInfo{i})
    candNames{end+1} = 'findDesiredInfo';
    candTinit{end+1} = T_findDesiredInfo{i};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where are you going to save this stuff? You will want to save it as it
% runs so you can pick up where you left off
% e.g. <saveSubdir>/forward/step_003/start_T_prev.mat
passDir = fullfile(saveSubdir, char(passName));
if ~exist(passDir, 'dir'); mkdir(passDir); end

stepDir = fullfile(passDir, sprintf('step_%03d', i));
if ~exist(stepDir, 'dir'); mkdir(stepDir); end

%%%%%%%%%% Tolerances for constraints %%%%%%%%%%
pctTol      = 0.05;
absTolFloor = 5e-3;

% Run each candidate start and collect results
nCands = numel(candNames);

emptySoln = struct( ...
    'name',            '', ...
    'T',               [], ...
    'LMS',             [], ...
    'rgbLin',          [], ...
    'infoNorm',        NaN, ...
    'distNorm',        NaN, ...
    'loss',            Inf, ...
    'nonLinSatisfied', false, ...
    'violation',       Inf);

% Preallocate 
startPtSolns = repmat(emptySoln, 1, nCands);

% Loop through each of your starting points
for k = 1:nCands
    thisName = candNames{k};
    T_init   = candTinit{k};

    % Force it to be identity if it's the first step
    if i == 1
        T_init = eye(3);
    end

    % Convert the current target into the correct argument position
    if strcmpi(sweepAxis,'info')
        lamArg = []; tgtInfoArg = thisX; tgtDistArg = [];
    elseif strcmpi(sweepAxis,'distortion')
        lamArg = []; tgtInfoArg = []; tgtDistArg = thisX;
    else
        lamArg = thisX; tgtInfoArg = []; tgtDistArg = [];
    end

    % Safe filename for this start (avoid weird characters)
    safeName = regexprep(thisName, '[^A-Za-z0-9_]+', '_');
    candFile = fullfile(stepDir, sprintf('start_%s.mat', safeName));

    % Load if exists, otherwise run + save 
    loadedOK = false;

    if exist(candFile, 'file')
        S = load(candFile);

        % Expect a struct named candSoln and completed=true
        if isfield(S,'completed') && isequal(S.completed,true) && ...
                isfield(S,'candSoln') && isstruct(S.candSoln) && ...
                isfield(S.candSoln,'T') && ~isempty(S.candSoln.T)

            % Fill the struct array entry 
            startPtSolns(k) = S.candSoln;
            loadedOK = true;

            fprintf('  [%s][step %03d] LOADED start %-12s from %s\n', char(passName), i, thisName, candFile);
        end
    end

    % If you dont have it saved somewhere, then do the optimization from
    % that start point
    if ~loadedOK

        forceIdentityThisStep = (i == 1);

        % Run the optimizer from this start
        [LMS_k, rgbLin_k, T_k, ~, infoNorm_k, ~, distNorm_k, nonLinSatisfied_k] = colorCorrectionOptimize(sweepAxis, lamArg, tgtInfoArg, tgtDistArg, LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
            obj.infoParams, obj.distortionParams, ...
            infoNormalizer, distortionNormalizer, Disp, ...
            'T_init', T_init,...
            'skipFmincon', forceIdentityThisStep, ...
            'pctTol', pctTol, ...     
            'absTolFloor', absTolFloor);

        % Compute the loss that corresponds to this sweep objective
        switch lower(sweepAxis)
            case 'info'

                loss_k = lossFunction('lambda', 0.0, T_k(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);

                % Check the violation of constraints, but give it some
                % wiggle room (epsBand)
                % Absolute distance from target value
                rawDelta = abs(infoNorm_k - thisX);

                % Allowed tolerance band around the target
                epsBand = max(absTolFloor, pctTol * thisX);

                if rawDelta <= epsBand
                    % Within allowed tolerance = no violation
                    violation_k = 0;
                else
                    % Outside tolerance = violation equals excess distance
                    violation_k = rawDelta - epsBand;
                end

            case 'distortion'
                loss_k = lossFunction('lambda', 1.0, T_k(:), LMSCalFormat, imgParams, dichromatType, ...
                    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
                % violation_k = abs(distNorm_k - thisX);

                % Absolute distance from target value
                rawDelta = abs(distNorm_k - thisX);

                % Allowed tolerance band around the target
                epsBand = max(absTolFloor, pctTol * thisX);

                if rawDelta <= epsBand
                    % Within allowed tolerance = no violation
                    violation_k = 0;
                else
                    % Outside tolerance = violation equals excess distance
                    violation_k = rawDelta - epsBand;
                end

            otherwise % 'lambda'
                loss_k = lossFunction('lambda', thisX, T_k(:), LMSCalFormat, imgParams, dichromatType, ...
                    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
                violation_k = NaN;
        end

        % Store this candidate solution
        startPtSolns(k).name            = thisName;
        startPtSolns(k).T               = T_k;
        startPtSolns(k).LMS             = LMS_k;
        startPtSolns(k).rgbLin          = rgbLin_k;
        startPtSolns(k).infoNorm        = infoNorm_k;
        startPtSolns(k).distNorm        = distNorm_k;
        startPtSolns(k).loss            = loss_k;
        startPtSolns(k).nonLinSatisfied = logical(nonLinSatisfied_k);
        startPtSolns(k).violation       = violation_k;


        % Save so resume works
        candSoln  = startPtSolns(k);  
        completed = true;           
        timestamp = datestr(now);    
        save(candFile, 'candSoln', 'completed', 'timestamp', 'i', 'k', 'thisX', 'sweepAxis', '-v7.3');

        % debugPlot_CandidatesUsingSweepFigure(i, sweepAxis, thisX, startPtSolns(k), passName, saveSubdir);


        fprintf('  [%s] saved %s -> %s\n', char(passName), thisName, candFile);

    end
end

% Choose the best feasible: first filter by whether it satisfies the nonlinear constraint, then minimize loss
nonLinMask = [startPtSolns.nonLinSatisfied];

if any(nonLinMask)
    % Which solutions satisfied the nonlinear constraint?
    solnsFeas = startPtSolns(nonLinMask);
    % Out of those, which has the lowest loss?
    [~, idx]  = min([solnsFeas.loss]);
    % Grab that one 
    best      = solnsFeas(idx);
else
    % If nothing satisfies nonlinear constraints, fall back to smallest violation, then smallest loss
    [~, idxV] = min([startPtSolns.violation]);
    best      = startPtSolns(idxV);

    % Tie-break on loss if multiple equal violations 
    tied = find([startPtSolns.violation] == best.violation);
    if numel(tied) > 1
        [~, idxL] = min([startPtSolns(tied).loss]);
        best = startPtSolns(tied(idxL));
    end
end

end
