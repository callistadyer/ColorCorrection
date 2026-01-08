
% APPLY SAVED BEST TRANSFORMS FROM A SOURCE RESOLUTION TO A TARGET RESOLUTION
% This script reuses best.T learned on a low-resolution image and applies it to a higher-resolution version of the same image without re-optimizing.

clear;
clc;


% %%%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%

% Name of the image file to process
imgType = 'fruit.png';

% Type of dichromacy to simulate
dichromatType = 'Deuteranopia';

% Set type (usually 1 unless using Ishihara images)
setType = 1;

% Which sweep was run originally: 'info' or 'distortion' (must match the saved run folder)
sweepAxis = 'distortion';

% Number of steps in the saved sweep
nSteps = 20;

% Which sweep steps to visualize (can be 1:5, 1:20, [1 3 7], etc.)
stepsToShow = 1:20;

% Size of the image used during optimization
mSource = 61;
nSource = 61;

% Size of the image to apply the transform to
mTarget = 128;
nTarget = 128;

% Optional extra folder name if your saved runs include e.g. "FreshStart" (leave '' if none)
extraRunFolder = '';

% Whether to clamp RGB values to [0,1] after applying the transform (helps when the 61x61 solution is slightly out of gamut on 128x128)
doClamp = true;
printGamutStats = true;


% %%%%%%%%%%%%%%%%%%%%% Metrics (must match the saved run) %%%%%%%%%%%%%%%%%%%%% 

% Handle to the information metric function (MUST match the saved run folder naming)
infoFcn = @computeInfo_regressCIE;

% Parameters used by the information metric (MUST match the saved run folder naming)
infoParams = struct('predictingWhat','L,M,S','predictingFromWhat','L and S');

% Handle to the distortion metric function (MUST match the saved run folder naming)
distortionFcn = @computeDistortion_DE2000;

% Parameters used by the distortion metric (typically empty, but keep consistent with saved run)
distortionParams = struct();

% Rendering function used to simulate dichromat perception (used only for visualization here)
renderFcn = @DichromRenderLinear;


%%%%%%%%%%%%%%%%%%%%% Which steps to show? %%%%%%%%%%%%%%%%%%%%%

% Ensure step list is a row vector with unique values
stepsToShow = unique(stepsToShow(:).');

% Remove any invalid step indices
stepsToShow = stepsToShow(stepsToShow >= 1 & stepsToShow <= nSteps);

% Error if nothing remains after validation
assert(~isempty(stepsToShow), 'stepsToShow is empty after validation.');

% Number of steps that will actually be displayed
nShow = numel(stepsToShow);


% %%%%%%%%%%%%%%%%%%%%% Build paths to the saved run %%%%%%%%%%%%%%%%%%%%%

% Name of the project (used for MATLAB preferences)
projectName = 'ColorCorrection';

% Base output directory used by computeSweep
outputDir = getpref(projectName, 'outputDir');

% Root directory where transformed images are stored
saveBase = fullfile(outputDir, 'testImagesTransformed');

% Folder name encoding the metric choices (must match computeSweep)
metricFolder = buildMetricFolderName(infoFcn, infoParams, distortionFcn);

% Folder name encoding the sweep axis and number of steps (must match computeSweep)
runFolder = sprintf('%s_%dsteps', lower(char(sweepAxis)), nSteps);

% Construct the source path name (matching computeSweep conventions)
if isempty(extraRunFolder)

    % Path without extra nesting
    pathNameSource = fullfile(dichromatType, imgType, sprintf('s%d_m%d_n%d', setType, mSource, nSource));

else

    % Path with extra nesting (e.g. FreshStart)
    pathNameSource = fullfile(dichromatType, imgType, extraRunFolder, sprintf('s%d_m%d_n%d', setType, mSource, nSource));

end

% Full directory of the saved sweep run
saveSubdirSource = fullfile(saveBase, pathNameSource, metricFolder, runFolder);

% Directory containing per-step best solutions
bestSubdirSource = fullfile(saveSubdirSource, 'best');

% Error if the best directory does not exist
assert(exist(bestSubdirSource, 'dir') == 7, sprintf('Could not find saved best directory:\n%s', bestSubdirSource));


% %%%%%%%%%%%%%%%%%%%%% Generate the target (high res) image %%%%%%%%%%%%%%%%%%%%%

% Do not clear existing cached data inside colorCorrectionGenerateImages
clearFlag = 0;

% Generate the higher-resolution target image (no optimization is run here)
[LMS_target, ~, ~, ~, Disp, imgParams_target, ~] = colorCorrectionGenerateImages(imgType, setType, mTarget, nTarget, dichromatType, clearFlag);

% Convert target LMS to linear RGB (CalFormat, 3xN)
RGBCalFormat_target = Disp.M_cones2rgb * LMS_target;

% Convert linear RGB to RGB contrast (THIS is where T is applied)
RGBContrast_target = (RGBCalFormat_target - Disp.grayRGB) ./ Disp.grayRGB;

% Convert target LMS to LMS contrast (used by info metric)
LMSContrast_target = (LMS_target - Disp.grayLMS) ./ Disp.grayLMS;


% %%%%%%%%%%%%%%%%%%%%% Compute normalizers to get the achieved metrics %%%%%%%%%%%%%%%%%%%%%

[protLMS, ~, ~] = DichromRenderLinear(LMS_target, 'Protanopia', Disp);
[deutLMS, ~, ~] = DichromRenderLinear(LMS_target, 'Deuteranopia', Disp);
[tritLMS, ~, ~] = DichromRenderLinear(LMS_target, 'Tritanopia', Disp);

% Build a reference LMS image by mixing the available cone signals from each dichromat simulation
LMS_ref = [protLMS(1, :); deutLMS(2, :); tritLMS(3, :)];

% Convert the reference LMS image into LMS contrast space
LMSContrast_ref = (LMS_ref - Disp.grayLMS) ./ Disp.grayLMS;

% Compute the information normalizer (passing normalizingValue = 1 returns the raw scale)
infoNormalizer = infoFcn(LMSContrast_target, LMSContrast_ref, imgParams_target, dichromatType, 1, Disp, infoParams);

% Compute the distortion normalizer (passing normalizingValue = 1 returns the raw scale)
distortionNormalizer = distortionFcn(LMS_target, LMS_ref, imgParams_target, 1, Disp, distortionParams);


% %%%%%%%%%%%%%%%%%%%%% Preallocate %%%%%%%%%%%%%%%%%%%%%

% Storage for trichromat images that will be displayed (m x n x 3 x nShow)
rgbTriAll = zeros(mTarget, nTarget, 3, nShow);
% Storage for dichromat-rendered images that will be displayed (m x n x 3 x nShow)
rgbDiAll = zeros(mTarget, nTarget, 3, nShow);
tileLabels = strings(1, nShow);


% %%%%%%%%%%%%%%%%%%%%% Apply the saved transforms to the target image %%%%%%%%%%%%%%%%%%%%%

% Loop over each step
for jj = 1:nShow

    % Convert from "display index" jj to the actual sweep step i
    i = stepsToShow(jj);

    % Build the path to the saved best.mat file for this step
    bestFile = fullfile(bestSubdirSource, sprintf('step_%03d', i), 'best.mat');

    % Verify the best.mat file exists
    assert(exist(bestFile, 'file') == 2, sprintf('Missing best.mat:\n%s', bestFile));

    % Load the saved structure (contains S.best and usually S.thisX)
    S = load(bestFile);

    % Extract the saved 3x3 transform matrix that was optimized on the 61x61 image
    T = S.best.T;

    % Apply the transform to RGB CONTRAST (critical: this matches lossFunction and colorCorrectionOptimize)
    newRGBContrast = T * RGBContrast_target;

    % Convert back to linear RGB by re-adding gray (undo contrast normalization)
    newRGBCalFormat = newRGBContrast .* Disp.grayRGB + Disp.grayRGB;

    % Optionally print out-of-gamut summary diagnostics before clamping
    if printGamutStats
        fprintf('step %02d: min=%g max=%g frac<0=%.3f frac>1=%.3f\n', i, min(newRGBCalFormat, [], 'all'), max(newRGBCalFormat, [], 'all'), mean(newRGBCalFormat(:) < 0), mean(newRGBCalFormat(:) > 1));
    end

    % Optionally clamp RGB into [0,1] to keep the image in gamut for visualization
    if doClamp
        newRGBCalFormat(newRGBCalFormat < 0) = 0;
        newRGBCalFormat(newRGBCalFormat > 1) = 1;
    end

    % Convert transformed RGB back to LMS (so distortion uses the correct domain)
    newLMS = Disp.M_rgb2cones * newRGBCalFormat;

    % Convert transformed LMS to LMS contrast (so info uses the correct domain)
    newLMSContrast = (newLMS - Disp.grayLMS) ./ Disp.grayLMS;

    % Render the transformed LMS image for the dichromat (returns linear RGB)
    [~, rgbLin_di, ~] = renderFcn(newLMS, dichromatType, Disp);

    % Convert the transformed trichromat linear RGB into display RGB for viewing
    rgbTriDisp = rgbLin2RGB(newRGBCalFormat, Disp);

    % Store the trichromat image as an m x n x 3 image for montage/tiled display
    rgbTriAll(:, :, :, jj) = CalFormatToImage(rgbTriDisp, mTarget, nTarget);

    % Convert the dichromat-rendered linear RGB into display RGB for viewing
    rgbDiDisp = rgbLin2RGB(rgbLin_di, Disp);

    % Store the dichromat image as an m x n x 3 image for montage/tiled display
    rgbDiAll(:, :, :, jj) = CalFormatToImage(rgbDiDisp, mTarget, nTarget);

    % If the saved sweep was an INFO sweep, label target info and achieved info for this high-res application
    if strcmpi(sweepAxis, 'info')

        % Pull the target info value that was used at this sweep step (saved in best.mat)
        targetVal = S.thisX;

        % Recompute achieved info at the target resolution using the same normalizer definition
        [~, achievedVal] = infoFcn(LMSContrast_target, newLMSContrast, imgParams_target, dichromatType, infoNormalizer, Disp, infoParams);

        % Construct the tile label (step index + target + achieved)
        tileLabels(jj) = sprintf('step %d\ninfo target=%.3g\ninfo achieved=%.3g', i, targetVal, achievedVal);

    else

        % Pull the target distortion value that was used at this sweep step (saved in best.mat)
        targetVal = S.thisX;

        % Recompute achieved distortion at the target resolution using the same normalizer definition
        [~, achievedVal] = distortionFcn(LMS_target, newLMS, imgParams_target, distortionNormalizer, Disp, distortionParams);

        % Construct the tile label (step index + target + achieved)
        tileLabels(jj) = sprintf('step %d\ndist target=%.3g\ndist achieved=%.3g', i, targetVal, achievedVal);

    end
end


%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%

nCols = ceil(sqrt(nShow));
nRows = ceil(nShow / nCols);

infoName = func2str(infoFcn);
distName = func2str(distortionFcn);

% Initialize the printed info-parameter string
infoParamStr = "";

% Only attempt to print params if infoParams is a struct
if isstruct(infoParams)

    % If predictingWhat is provided, include it in the figure title
    if isfield(infoParams, 'predictingWhat')
        infoParamStr = infoParamStr + " predictingWhat=" + string(infoParams.predictingWhat);
    end

    % If predictingFromWhat is provided, include it in the figure title
    if isfield(infoParams, 'predictingFromWhat')
        infoParamStr = infoParamStr + " predictingFromWhat=" + string(infoParams.predictingFromWhat);
    end
end

% Create the trichromat figure
figure('Color', 'w', 'Name', 'Trichromat labeled grid', 'NumberTitle', 'off');

% Create the tiled layout for the trichromat figure
tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add figure-level title describing metrics and sweep settings (no clamp info printed)
sgtitle(sprintf('%s | sweep=%s | %dx%d → %dx%d\nINFO: %s%s\nDIST: %s', dichromatType, lower(string(sweepAxis)), mSource, nSource, mTarget, nTarget, infoName, infoParamStr, distName), 'Interpreter', 'none');

% Plot each trichromat image with its label
for jj = 1:nShow
    nexttile;
    imagesc(rgbTriAll(:, :, :, jj));
    axis image off;
    title(tileLabels(jj), 'Interpreter', 'none', 'FontSize', 9);
end

% Create the dichromat figure
figure('Color', 'w', 'Name', 'Dichromat labeled grid', 'NumberTitle', 'off');

% Create the tiled layout for the dichromat figure
tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add figure-level title describing metrics and sweep settings (no clamp info printed)
sgtitle(sprintf('%s | sweep=%s | %dx%d → %dx%d\nINFO: %s%s\nDIST: %s', dichromatType, lower(string(sweepAxis)), mSource, nSource, mTarget, nTarget, infoName, infoParamStr, distName), 'Interpreter', 'none');

% Plot each dichromat image with its label
for jj = 1:nShow
    nexttile;
    imagesc(rgbDiAll(:, :, :, jj));
    axis image off;
    title(tileLabels(jj), 'Interpreter', 'none', 'FontSize', 9);
end

