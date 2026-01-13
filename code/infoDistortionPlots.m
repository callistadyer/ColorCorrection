
function [infoNorm, distortionNorm] = infoDistortionPlots(sweepFile, nSteps, Save, bPLOT)
% infoDistortionPlots  Load sweepOutputs.mat and plot info vs distortion
%
% Inputs
%   sweepFile - full path to sweepOutputs.mat
%   nSteps    - number of steps you expect to plot/return
%   Save      - true/false save png next to sweepOutputs.mat
%   bPLOT     - 1/0 whether to plot
%
% Outputs
%   infoNorm
%   distortionNorm (1 x nSteps)

if nargin < 5, bPLOT = 1; end
if nargin < 4, Save  = false; end

sweepFile = char(sweepFile);

if ~exist(sweepFile,'file')
    error('infoDistortionPlots:FileNotFound', 'Could not find sweep file:\n  %s', sweepFile);
end

[folder, ~, ~] = fileparts(sweepFile);

fprintf('[infoDistortionPlots] Loading:\n  %s\n', sweepFile);

% Load outputs
[outputs,~,~,~] = loadSweepOutputs(sweepFile, nSteps, 'raw');

infoNorm       = nan(1, nSteps);
distortionNorm = nan(1, nSteps);

% Extract info / distortion values
for i = 1:min(nSteps, numel(outputs))
    s = outputs{i};
    if isempty(s) || ~isstruct(s)
        continue;
    end

    % distortion
    if isfield(s,'distortionNormalizedAchievedSweep') && ~isempty(s.distortionNormalizedAchievedSweep)
        distortionNorm(i) = s.distortionNormalizedAchievedSweep;
    elseif isfield(s,'distortionNormalized') && ~isempty(s.distortionNormalized)
        distortionNorm(i) = s.distortionNormalized;
    end

    % info
    if isfield(s,'infoNormalizedAchievedSweep') && ~isempty(s.infoNormalizedAchievedSweep)
        infoNorm(i) = s.infoNormalizedAchievedSweep;
    elseif isfield(s,'infoNormalized') && ~isempty(s.infoNormalized)
        infoNorm(i) = s.infoNormalized;
    end
end

% Plot
if bPLOT == 1
    figure('Color','w','Visible','on');

    plot(distortionNorm, infoNorm, '-o','LineWidth',2, 'MarkerSize',8, 'MarkerFaceColor','w', 'Color','k');
    
    % Make it pretty
    grid on;
    axis square;
    xlabel('distortionNormalized');
    ylabel('infoNormalized');
end

% Save figure next to sweepOutputs.mat
if Save
    outPath = fullfile(folder, 'info_vs_distortion.png');
    try
        exportgraphics(gca, outPath, ...
            'Resolution',200,'BackgroundColor','w');
    catch
        set(gcf,'PaperPositionMode','auto');
        print(gcf, outPath, '-dpng','-r200');
    end
    fprintf('[infoDistortionPlots] wrote: %s\n', outPath);
end

end
