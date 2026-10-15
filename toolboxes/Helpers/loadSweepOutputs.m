function [outputs, startStep, allDone, loaded] = loadSweepOutputs(saveFile, nSteps)
% loadSweepOutputs
% Load a previous sweepOutputs.mat if present, determine resume step,
% and get final sweep arrays 
%
% Inputs
%   saveFile   - full path to sweepOutputs.mat
%   nSteps     - expected number of steps
%
% Outputs
%   outputs   - 1 x nSteps cell array (completed steps copied in; others empty)
%   startStep - first step index to run (lastDone+1), in [1..nSteps+1]
%   allDone   - true if all steps 1..nSteps are done
%   loaded    - if allDone==true 

% Defaults
outputs   = cell(1, nSteps);
startStep = 1;
allDone   = false;
loaded    = [];

% Is the step done?
isDoneStep = @(s) isstruct(s) ...
    && isfield(s,'stepCompleted') && isequal(s.stepCompleted,true) ...
    && isfield(s,'transformRGBmatrix') && ~isempty(s.transformRGBmatrix) ...
    && isfield(s,'infoNormalizedAchievedSweep') && ~isnan(s.infoNormalizedAchievedSweep) ...
    && isfield(s,'distortionNormalizedAchievedSweep') && ~isnan(s.distortionNormalizedAchievedSweep);

% No file
if ~exist(saveFile, 'file')
    fprintf('computeSweep: No existing sweep file. Starting fresh.\n');
    return;
end

fprintf('computeSweep: Found existing sweep file:\n  %s\n', saveFile);

S = load(saveFile, 'outputs');
if ~isfield(S,'outputs') || isempty(S.outputs)
    fprintf('  sweepOutputs.mat exists but no usable outputs found. Starting fresh.\n');
    return;
end

old = S.outputs;

% Copy completed steps
for ii = 1:min(nSteps, numel(old))
    if ~isempty(old{ii}) && isDoneStep(old{ii})
        outputs{ii} = old{ii};
    end
end

% Determine resume point 
doneMask = false(1, nSteps);
for ii = 1:nSteps
    doneMask(ii) = ~isempty(outputs{ii}) && isDoneStep(outputs{ii});
end

lastDone = find(doneMask, 1, 'last');
if ~isempty(lastDone)
    startStep = lastDone + 1;
end

allDone = all(doneMask);

% Early return if complete 
if allDone
    fprintf('computeSweep: All %d steps already complete. Loading and returning.\n', nSteps);

    loaded.LMSDaltonizedCalFormatSweep            = cell(1,nSteps);
    loaded.rgbLinDaltonizedCalFormatSweep         = cell(1,nSteps);
    loaded.LMSDaltonizedRenderedCalFormatSweep    = cell(1,nSteps);
    loaded.rgbLinDaltonizedRenderedCalFormatSweep = cell(1,nSteps);
    loaded.transformRGBmatrixSweep                = cell(1,nSteps);

    loaded.targetInfoNormalized                   = nan(1,nSteps);
    loaded.targetDistortionNormalized             = nan(1,nSteps);
    loaded.infoNormalizedAchievedSweep            = nan(1,nSteps);
    loaded.distortionNormalizedAchievedSweep      = nan(1,nSteps);

    for ii = 1:nSteps
        s = outputs{ii};

        loaded.LMSDaltonizedCalFormatSweep{ii}            = s.LMSDaltonizedCalFormat;
        loaded.rgbLinDaltonizedCalFormatSweep{ii}         = s.rgbLinDaltonizedCalFormat;
        loaded.LMSDaltonizedRenderedCalFormatSweep{ii}    = s.LMSDaltonizedRenderedCalFormat;
        loaded.rgbLinDaltonizedRenderedCalFormatSweep{ii} = s.rgbLinDaltonizedRenderedCalFormat;
        loaded.transformRGBmatrixSweep{ii}                = s.transformRGBmatrix;

        if isfield(s,'targetInfoNormalized')
            loaded.targetInfoNormalized(ii) = s.targetInfoNormalized;
        end
        if isfield(s,'targetDistortionNormalized')
            loaded.targetDistortionNormalized(ii) = s.targetDistortionNormalized;
        end

        loaded.infoNormalizedAchievedSweep(ii)       = s.infoNormalizedAchievedSweep;
        loaded.distortionNormalizedAchievedSweep(ii) = s.distortionNormalizedAchievedSweep;
    end
    return;
end

% Resume message
if isempty(lastDone)
    fprintf('  No completed steps detected. Starting from step 1.\n');
else
    fprintf('  Completed steps: 1–%d. Resuming at step %d.\n', lastDone, startStep);
end

end
