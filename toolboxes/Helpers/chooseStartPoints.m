function [best, startPtSolns, candNames, candTinit] = chooseStartPoints( ...
    i, sweepAxis, thisX, ...
    T_I, T_prev, ...
    T_infoSeeds, T_findDesiredDist, T_findDesiredInfo, ...
    LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp)
% chooseStartPoints  Build a list of candidate start points and run the optimizer from each, then pick the best.
%
% Syntax:
%   [best, startPtSolns] = chooseStartPoints(i, sweepAxis, thisX, ...
%       T_I, T_prev, T_infoSeeds, T_findDesiredDist, T_findDesiredInfo, ...
%       LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp)
%
% Description:
%   This helper takes all possible start-point sources (identity, previous step, info-sweep seed,
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

% Build candidate name list
candNames = {};
candTinit = {};

% Always try identity (baseline)
candNames{end+1} = 'T_I';
candTinit{end+1} = T_I;

% Always try chaining from previous step
candNames{end+1} = 'T_prev';
candTinit{end+1} = T_prev;

% Try info-sweep solution at this step (if provided)
if exist('T_infoSeeds','var') && ~isempty(T_infoSeeds) && numel(T_infoSeeds) >= i && ~isempty(T_infoSeeds{i})
    candNames{end+1} = 'infoSoln';
    candTinit{end+1} = T_infoSeeds{i};
end

% % Try feasibility-start that hits desired distortion (if provided)
% if exist('T_findDesiredDist','var') && ~isempty(T_findDesiredDist) && numel(T_findDesiredDist) >= i && ~isempty(T_findDesiredDist{i})
%     candNames{end+1} = 'findDesiredDist';
%     candTinit{end+1} = T_findDesiredDist{i};
% end
% 
% % Try feasibility-start that hits desired info (if provided)
% if exist('T_findDesiredInfo','var') && ~isempty(T_findDesiredInfo) && numel(T_findDesiredInfo) >= i && ~isempty(T_findDesiredInfo{i})
%     candNames{end+1} = 'findDesiredInfo';
%     candTinit{end+1} = T_findDesiredInfo{i};
% end

% Run each candidate start and collect results
nCands = numel(candNames);

% Template struct (scalar!)
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

% Preallocate a 1 x nCands struct array
startPtSolns = repmat(emptySoln, 1, nCands);
for k = 1:nCands
    thisName = candNames{k};
    T_init   = candTinit{k};

    % Convert the current target into the correct argument position
    if strcmpi(sweepAxis,'info')
        lamArg = []; tgtInfoArg = thisX; tgtDistArg = [];
    elseif strcmpi(sweepAxis,'distortion')
        lamArg = []; tgtInfoArg = []; tgtDistArg = thisX;
    else
        lamArg = thisX; tgtInfoArg = []; tgtDistArg = [];
    end

    % Run the optimizer from this start
    [LMS_k, rgbLin_k, T_k, ~, infoNorm_k, ~, distNorm_k, nonLinSatisfied_k] = colorCorrectionOptimize( ...
        sweepAxis, lamArg, tgtInfoArg, tgtDistArg, ...
        LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        obj.infoParams, obj.distortionParams, ...
        infoNormalizer, distortionNormalizer, Disp, ...
        'T_init', T_init);

    % Compute the loss that corresponds to this sweep objective
    switch lower(sweepAxis)
        case 'info'
            loss_k = lossFunction('lambda', 0.0, T_k(:), LMSCalFormat, imgParams, dichromatType, ...
                obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
            violation_k = abs(infoNorm_k - thisX);

        case 'distortion'
            loss_k = lossFunction('lambda', 1.0, T_k(:), LMSCalFormat, imgParams, dichromatType, ...
                obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, obj.infoParams);
            violation_k = abs(distNorm_k - thisX);

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
end

% Choose the best feasible: first filter by nonlinear satisfaction, then minimize loss
nonLinMask = [startPtSolns.nonLinSatisfied];

if any(nonLinMask)
    solnsFeas = startPtSolns(nonLinMask);
    [~, idx]  = min([solnsFeas.loss]);
    best      = solnsFeas(idx);
else
    % If nothing satisfies nonlinear constraints, fall back to smallest violation, then smallest loss
    [~, idxV] = min([startPtSolns.violation]);
    best      = startPtSolns(idxV);

    % Tie-break on loss if multiple equal violations (rare, but safe)
    tied = find([startPtSolns.violation] == best.violation);
    if numel(tied) > 1
        [~, idxL] = min([startPtSolns(tied).loss]);
        best = startPtSolns(tied(idxL));
    end
end

end
