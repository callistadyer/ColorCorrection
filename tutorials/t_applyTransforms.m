
% Compare a single sweep of info or distortion across multiple metrics

clear;
clc;

% Settings
imgType = 'flower1.png';
dichromatType = 'Deuteranopia';
setType = 1;
% Choose which sweep you are comparing across metrics:
%   'distortion' = distortion sweep results across info metrics
%   'info'       = info sweep results across info metrics
sweepAxis = 'distortion';

% How many steps?
nSteps = 20;

% Which sweep steps to show (columns)
stepsToShow = 1:7;

% Resolution used during optimization (must match saved folder)
mSource = 61;
nSource = 61;

% Resolution you want to apply transforms to
mTarget = 128*2;
nTarget = 128*2;

extraRunFolder = '';

% Clamp RGB into [0,1] after applying transform 
% Gotta do this because when we apply the low res transform to the high res
% image, sometimes there are pixels that go out of gamut
doClamp = true;

% Print min/max and out-of-gamut fraction each step
printGamutStats = false;

%%%%%%%%%%%%% Metrics %%%%%%%%%%%%%
% Which info metrics to compare
infoFcns = { ...
    @computeInfo_regressCIE, ...
    @computeInfo_regressSquared, ...
    @computeInfo_LMdifference ...

};

% Info params
infoParams = struct('predictingWhat','L,M,S','predictingFromWhat','L and S');

% Distortion metric
distortionFcn = @computeDistortion_DE2000;
distortionParams = struct();

% Render function (for visualization only)
% If you want, we can use the original brettel here. We will just need to
% do some clipping again
renderFcn = @DichromRenderLinear;

% Check that the steps you want to show actually exist
stepsToShow = unique(stepsToShow(:).');
stepsToShow = stepsToShow(stepsToShow >= 1 & stepsToShow <= nSteps);
assert(~isempty(stepsToShow), 'stepsToShow is empty after validation.');

nShow    = numel(stepsToShow);
nMetrics = numel(infoFcns);
assert(nMetrics >= 1, 'infoFcns is empty.');

%%%%%%%%%%%%%%%%%%%%%% Paths for grabing the files %%%%%%%%%%%%%%%%%%%%%%
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');

% '/Users/callista/Aguirre-Brainard Lab Dropbox/Callista Dyer/DALT_analysis/ColorCorrection/testImagesTransformed'
saveBase    = fullfile(outputDir, 'testImagesTransformed');

% e.g., 'distortion_20steps'
runFolder   = sprintf('%s_%dsteps', lower(char(sweepAxis)), nSteps);

% e.g., 'Deuteranopia/ishi45.png/s1_m61_n61'
if isempty(extraRunFolder)
    pathNameSource = fullfile(dichromatType, imgType, sprintf('s%d_m%d_n%d', setType, mSource, nSource));
else
    pathNameSource = fullfile(dichromatType, imgType, extraRunFolder, sprintf('s%d_m%d_n%d', setType, mSource, nSource));
end

% For each metric, computeSweep stored results under:
%   saveBase / pathNameSource / metricFolder / runFolder / best / step_### / best.mat
metricFolders = cell(1, nMetrics);
bestSubdirs = cell(1, nMetrics);

for k = 1:nMetrics
    % e.g., 'computeInfo_regressCIE__L_M_S_from_L_and_S__computeDistortion_DE2000'
    metricFolders{k} = buildMetricFolderName(infoFcns{k}, infoParams, distortionFcn);
    bestSubdirs{k} = fullfile(saveBase, pathNameSource, metricFolders{k}, runFolder, 'best');
    assert(exist(bestSubdirs{k}, 'dir') == 7, sprintf('Missing best folder for metric %d:\n%s', k, bestSubdirs{k}));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the high res image
% Apply saved transforms to this image (no optimization happens here)
clearFlag = 0;
[LMS_target, ~, ~, ~, Disp, imgParams_target, ~] = colorCorrectionGenerateImages(imgType, setType, mTarget, nTarget, dichromatType, clearFlag);

% Convert target LMS -> linear RGB 
RGBCalFormat_target = Disp.M_cones2rgb * LMS_target;

% Convert linear RGB -> RGB contrast (where T is applied)
RGBContrast_target = (RGBCalFormat_target - Disp.grayRGB) ./ Disp.grayRGB;

% Convert LMS -> LMS contrast (info metrics use that)
LMSContrast_target = (LMS_target - Disp.grayLMS) ./ Disp.grayLMS;

% Build the normalizier
[protLMS, ~, ~] = DichromRenderLinear(LMS_target, 'Protanopia', Disp);
[deutLMS, ~, ~] = DichromRenderLinear(LMS_target, 'Deuteranopia', Disp);
[tritLMS, ~, ~] = DichromRenderLinear(LMS_target, 'Tritanopia', Disp);

LMS_ref = [protLMS(1,:); deutLMS(2,:); tritLMS(3,:)];
LMSContrast_ref = (LMS_ref - Disp.grayLMS) ./ Disp.grayLMS;

% Info normalizers are metric-specific
infoNormalizers = cell(1, nMetrics);
for k = 1:nMetrics
    infoNormalizers{k} = infoFcns{k}(LMSContrast_target, LMSContrast_ref, imgParams_target, dichromatType, 1, Disp, infoParams);
end

% Distortion normalizer is shared
distortionNormalizer = distortionFcn(LMS_target, LMS_ref, imgParams_target, 1, Disp, distortionParams);

rgbTriImg = cell(nShow, nMetrics);
rgbDiImg = cell(nShow, nMetrics);
tileLabels = strings(nShow, nMetrics);

% Load and apply the transformations
for jj = 1:nShow
    % Grab the steps you wanna show
    i = stepsToShow(jj);

    for k = 1:nMetrics

        % Load best solution for this metric and this step
        bestFile = fullfile(bestSubdirs{k}, sprintf('step_%03d', i), 'best.mat');
        assert(exist(bestFile,'file')==2, sprintf('Missing best.mat:\n%s', bestFile));

        S = load(bestFile);

        % Extract 3x3 transform
        T = S.best.T;

        % thisX is the target value used in the sweep (info or distortion)
        targetX = S.thisX;

        % Apply transform to RGB contrast
        newRGBContrast = T * RGBContrast_target;

        % Convert back to linear RGB
        newRGBLinCalFormat = newRGBContrast .* Disp.grayRGB + Disp.grayRGB;

        % If you want, show how much we have to clip
        if printGamutStats
            fprintf('step %02d metric %02d: min=%g max=%g frac<0=%.3f frac>1=%.3f\n', i, k, min(newRGBLinCalFormat, [], 'all'), max(newRGBLinCalFormat, [], 'all'), mean(newRGBLinCalFormat(:) < 0), mean(newRGBLinCalFormat(:) > 1));
        end

        % Cut the stuff to keep in gamut
        if doClamp
            newRGBLinCalFormat = min(max(newRGBLinCalFormat, 0), 1);
        end

        % Convert transformed RGB -> LMS
        newLMSCalFormat = Disp.M_rgb2cones * newRGBLinCalFormat;

        % Convert transformed LMS -> LMS contrast (for info)
        newLMSContrastCalFormat = (newLMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

        % Render transformed LMS for the dichromat
        [~, rgbLin_di, ~] = renderFcn(newLMSCalFormat, dichromatType, Disp);

        % Convert to display RGB
        RGBTriCalFormat = rgbLin2RGB(newRGBLinCalFormat, Disp);
        RGBDiCalFormat  = rgbLin2RGB(rgbLin_di, Disp);

        % Store image-formatted versions for imagesc
        rgbTriImg{jj,k} = CalFormatToImage(RGBTriCalFormat, mTarget, nTarget);
        rgbDiImg{jj,k}  = CalFormatToImage(RGBDiCalFormat,  mTarget, nTarget);

        % Compute achieved distortion
        [~, achievedDist] = distortionFcn(LMS_target, newLMSCalFormat, imgParams_target, distortionNormalizer, Disp, distortionParams);

        % Compute achieved info
        [~, achievedInfo] = infoFcns{k}(LMSContrast_target, newLMSContrastCalFormat, imgParams_target, dichromatType, infoNormalizers{k}, Disp, infoParams);

        if strcmpi(sweepAxis,'distortion')
            tileLabels(jj,k) = sprintf('step %d\nDist target=%.4g\nDist achieved=%.4g\nInfo achieved=%.4g', ...
                i, targetX, achievedDist, achievedInfo);
        else
            tileLabels(jj,k) = sprintf('step %d\nInfo target=%.4g\nInfo achieved=%.4g\nDist achieved=%.4g', ...
                i, targetX, achievedInfo, achievedDist);
        end
    end
end

% Make titles for the graph
metricNames = cellfun(@func2str, infoFcns, 'UniformOutput', false);
metricsTitleStr = strjoin(metricNames, ' | ');

% FIGURE 1: Trichromat grid
figTri = figure('Position',[4 420 1714 607],'Color','w','Name','Trichromat','NumberTitle','off');
tTri = tiledlayout(nMetrics, nShow, 'TileSpacing','compact', 'Padding','compact');
tTri.Padding = 'loose';

sgtitle(sprintf('%s | sweep=%s | %dx%d -> %dx%d | %s', dichromatType, lower(string(sweepAxis)), mSource, nSource, mTarget, nTarget, metricsTitleStr), 'Interpreter','none');

for k = 1:nMetrics
    for c = 1:nShow
        ax = nexttile((k-1)*nShow + c);
        imagesc(ax, rgbTriImg{c,k});
        axis(ax,'image');

        ax.XTick = [];
        ax.YTick = [];
        ax.Box = 'off';

        title(ax, tileLabels(c,k), 'Interpreter','none', 'FontSize', 9);
    end
end

for k = 1:nMetrics
    leftTileIdx = (k-1)*nShow + 1;
    ax = nexttile(leftTileIdx);

    ylabel(ax, metricNames{k}, 'FontWeight','bold', 'FontSize', 12, 'Interpreter','none', 'Color','k');

    ax.Clipping = 'off';
    ax.YLabel.Clipping = 'off';

    ax.YLabel.Units = 'normalized';
    ax.YLabel.Position(1) = -0.12;
end

% FIGURE 2: Dichromat grid
figDi = figure('Position',[4 302 1714 607],'Color','w','Name','Dichromat','NumberTitle','off');
tDi = tiledlayout(nMetrics, nShow, 'TileSpacing','compact', 'Padding','compact');
tDi.Padding = 'loose';

sgtitle(sprintf('%s render | sweep=%s | %dx%d -> %dx%d | %s', ...
    dichromatType, lower(string(sweepAxis)), mSource, nSource, mTarget, nTarget, metricsTitleStr), ...
    'Interpreter','none');

for k = 1:nMetrics
    for c = 1:nShow
        ax = nexttile((k-1)*nShow + c);
        imagesc(ax, rgbDiImg{c,k});
        axis(ax,'image');

        ax.XTick = [];
        ax.YTick = [];
        ax.Box = 'off';

        title(ax, tileLabels(c,k), 'Interpreter','none', 'FontSize', 9);
    end
end

for k = 1:nMetrics
    leftTileIdx = (k-1)*nShow + 1;
    ax = nexttile(leftTileIdx);

    ylabel(ax, metricNames{k}, 'FontWeight','bold', 'FontSize', 12, 'Interpreter','none', 'Color','k');

    ax.Clipping = 'off';
    ax.YLabel.Clipping = 'off';

    ax.YLabel.Units = 'normalized';
    ax.YLabel.Position(1) = -0.12;
end
