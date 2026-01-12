
function [infoNorm, distortionNorm] = infoDistortionPlots(sweepFile, nSteps, Live, Save, bPLOT)
% infoDistortionPlots  Load ONE specific sweepOutputs.mat and plot info vs distortion
% (uses RAW loading mode: no isDoneStep logic, no recursion, no guessing)
%
% Inputs
%   sweepFile : full path to sweepOutputs.mat
%   nSteps    : number of steps you expect to plot/return
%   Live      : true/false show figure
%   Save      : true/false save png next to sweepOutputs.mat
%   bPLOT     : 1/0 whether to plot
%
% Outputs
%   infoNorm, distortionNorm : 1 x nSteps arrays (NaN where missing)

if nargin < 5, bPLOT = 1; end
if nargin < 4, Save  = false; end
if nargin < 3, Live  = true; end

sweepFile = char(sweepFile);

if ~exist(sweepFile,'file')
    error('infoDistortionPlots:FileNotFound', ...
        'Could not find sweep file:\n  %s', sweepFile);
end

[folder, name, ext] = fileparts(sweepFile);
if ~strcmpi([name ext], 'sweepOutputs.mat')
    error('infoDistortionPlots:NotSweepOutputs', ...
        'Expected a file named sweepOutputs.mat, got:\n  %s', sweepFile);
end

fprintf('[infoDistortionPlots] Loading:\n  %s\n', sweepFile);

% ------------------------------------------------------------
% Load outputs using RAW mode (no isDoneStep logic)
% ------------------------------------------------------------
[outputs,~,~,~] = loadSweepOutputs(sweepFile, nSteps, 'raw');

% ------------------------------------------------------------
% Preallocate
% ------------------------------------------------------------
infoNorm       = nan(1, nSteps);
distortionNorm = nan(1, nSteps);

% ------------------------------------------------------------
% Extract info / distortion values if present
% ------------------------------------------------------------
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

% ------------------------------------------------------------
% Plot
% ------------------------------------------------------------
if bPLOT == 1
    figVis = 'on';
    if ~Live
        figVis = 'off';
    end

    figure('Color','w','Visible',figVis);
    plot(distortionNorm, infoNorm, '-o', ...
        'LineWidth',2, 'MarkerSize',8, ...
        'MarkerFaceColor','w', 'Color','k');
    grid on;
    axis square;
    xlabel('distortionNormalized');
    ylabel('infoNormalized');
    title(sprintf('Loaded: %s', folder), 'Interpreter','none');
end

% ------------------------------------------------------------
% Save figure next to sweepOutputs.mat
% ------------------------------------------------------------
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

% function [infoNorm, distortionNorm] = infoDistortionPlots(sweepFile, nSteps, Live, Save, bPLOT)
% % infoDistortionPlots_strict  Load ONE specific sweepOutputs.mat (no recursion, no guessing)
% %
% % Inputs
% %   sweepFile : full path to sweepOutputs.mat
% %   nSteps    : number of steps you expect to plot/return
% %   Live      : true/false show figure
% %   Save      : true/false save png next to sweepOutputs.mat
% %   bPLOT     : 1/0 whether to plot
% %
% % Outputs
% %   infoNorm, distortionNorm : 1 x nSteps arrays (NaN where missing)
% 
% if nargin < 5, bPLOT = 1; end
% if nargin < 4, Save  = false; end
% if nargin < 3, Live  = true; end
% 
% sweepFile = char(sweepFile);
% 
% if ~exist(sweepFile,'file')
%     error('infoDistortionPlots:FileNotFound', ...
%         'Could not find sweep file:\n  %s', sweepFile);
% end
% 
% [folder, name, ext] = fileparts(sweepFile);
% if ~strcmpi([name ext], 'sweepOutputs.mat')
%     error('infoDistortionPlots:NotSweepOutputs', ...
%         'Expected a file named sweepOutputs.mat, got:\n  %s', sweepFile);
% end
% 
% fprintf('[infoDistortionPlots_strict] Loading:\n  %s\n', sweepFile);
% 
% % ---- Load outputs ----
% S = load(sweepFile, 'outputs');
% if ~isfield(S,'outputs') || isempty(S.outputs)
%     error('infoDistortionPlots:BadSweepFile', ...
%         'File does not contain a nonempty variable named "outputs":\n  %s', sweepFile);
% end
% outputs = S.outputs;
% 
% % ---- Preallocate ----
% infoNorm       = nan(1, nSteps);
% distortionNorm = nan(1, nSteps);
% 
% % ---- Extract
% for i = 1:min(nSteps, numel(outputs))
%     s = outputs{i};
%     if isempty(s) || ~isstruct(s), continue; end
% 
%     % distortion
%     if isfield(s,'distortionNormalizedAchievedSweep') && ~isempty(s.distortionNormalizedAchievedSweep)
%         distortionNorm(i) = s.distortionNormalizedAchievedSweep;
%     elseif isfield(s,'distortionNormalized') && ~isempty(s.distortionNormalized)
%         distortionNorm(i) = s.distortionNormalized;
%     end
% 
%     % info
%     if isfield(s,'infoNormalizedAchievedSweep') && ~isempty(s.infoNormalizedAchievedSweep)
%         infoNorm(i) = s.infoNormalizedAchievedSweep;
%     elseif isfield(s,'infoNormalized') && ~isempty(s.infoNormalized)
%         infoNorm(i) = s.infoNormalized;
%     end
% end
% 
% % ---- Plot ----
% if bPLOT == 1
%     figVis = 'on';
%     if ~Live, figVis = 'off'; end
% 
%     figure('Color','w','Visible',figVis);
%     plot(distortionNorm, infoNorm, '-o', 'LineWidth',2, 'MarkerSize',8, ...
%         'MarkerFaceColor','w', 'Color','k');
%     grid on; axis square;
%     xlabel('distortionNormalized');
%     ylabel('infoNormalized');
%     title(sprintf('Loaded: %s', folder), 'Interpreter','none');
% end
% 
% % ---- Save next to the sweepOutputs.mat ----
% if Save
%     outPath = fullfile(folder, 'info_vs_distortion.png');
%     try
%         exportgraphics(gca, outPath, 'Resolution',200,'BackgroundColor','w');
%     catch
%         set(gcf,'PaperPositionMode','auto');
%         print(gcf, outPath, '-dpng','-r200');
%     end
%     fprintf('[infoDistortionPlots_strict] wrote: %s\n', outPath);
% end
% 
% end
