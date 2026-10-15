function t_compareInfoMetrics

% t_compareInfoMetrics.m
% Compare image enhancement outputs across multiple inf metrics while
% holding distortion approximately constant
%
% Overview
%   For each input image and each info metric:
%     (1) Build a daltonizer object using:
%           - infoFcn        (varies across metrics)
%           - distortionFcn  (fixed: e.g., computeDistortion_DE2000)
%           - renderFcn      (fixed: e.g., DichromRenderLinear)
%     (2) Run a distortion sweep (computeSweep with sweepAxis='distortion').
%         This produces a set of optimized transforms T along with:
%           - achieved distortion (normalized)
%           - achieved info       (normalized)
%           - optimized images (tri) and their dichromat renderings (di)
%     (3) Select the sweep point whose achieved distortion is closest to
%         targetDistToCompare (some chosen normalized distortion level)
%     (4) Store:
%           - the selected transform T
%           - achieved info/distortion at the selected point
%           - full sweep curves (info vs distortion)
%           - output images (tri + dichromat render)
%     (5) Plot:
%           - A montage of the original (tri + di) image plus each metric's
%             optimized output (tri + di)
%           - A plot of achieved info vs achieved distortion for each metric
%
% Key idea
%   This script tries to isolate the effect of the info metrics by matching
%   solutions at (approximately) the same distortion. 
%
% Inputs / configuration (edit in the "Images parameters" section)
%   imageTypes           - cell array of image filenames
%   dichromatType        - e.g., 'Deuteranopia'
%   m,n                  - output image size
%   nSteps               - number of points in the distortion sweep
%   targetDistToCompare  - normalized distortion level used for matching
%   infoFcnList          - list of info metric function handles
%   infoParamsList       - parameter structs (e.g., regress params)
%
% Outputs
%   results(ii) struct per image, containing:
%     .imgType, .dichromatType, .pathName, .targetDistToCompare
%     .metric(kk) with fields:
%        .name, .infoFcn, .infoParams
%        .idxClosest, .distAtTarget, .infoAtTarget, .absErrAtTarget
%        .T
%        .infoNormAch, .distNormAch, .targetDistNorm
%
%
% History
%   2025-12-22  CD  Wrote initial version to compare info metrics at fixed distortion


clear; clc;

%%%%%%%%%%%%%%%%%%%%%
% Images parameters %
%%%%%%%%%%%%%%%%%%%%%
imageTypes = {'flower1.png', 'fruit.png', 'Gaugin.png', 'ishi45.png'};
setType       = 1;
dichromatType = 'Deuteranopia';

% Size of image
m = 68; n = 68;

% How many steps to sample at 
nSteps = 11;

% Clear the previous image data if it existed already in dropbox?
clearFlag = 0;

% Distortion metric
distortionFcn    = @computeDistortion_DE2000;
distortionParams = struct();

% Render function (for dichromat rendering)
renderFcn    = @DichromRenderLinear;
renderParams = struct();

% Info metrics
infoFcnList = {@computeInfo_LMdifference, ...
    @computeInfo_regress, ...
    @computeInfo_Wade };

infoNameList = { ...
    'LMdifference', ...
    'regress', ...
    'Wade' };

% Some info metrics need params... regress does, others can use empty struct
infoParamsList = cell(size(infoFcnList));
for k = 1:numel(infoFcnList)
    if isequal(infoFcnList{k}, @computeInfo_regress)
        infoParamsList{k} = struct( ...
            'predictingWhat',     'L,M,S', ...
            'predictingFromWhat', 'L and S');
    else
        infoParamsList{k} = struct();
    end
end

% Fixed distortion to compare info types
% This is the distortion where we want to compare info and
% visuals across metrics.
%
% In practice I pick the sweep step whose achieved distortion is closest
targetDistToCompare = 0.80;

log = struct();

%% Main loop over images
for ii = 1:numel(imageTypes)

    imgType = imageTypes{ii};
    fprintf('\n====================================================\n');
    fprintf('IMAGE %d/%d: %s | %s\n', ii, numel(imageTypes), imgType, dichromatType);
    fprintf('====================================================\n');

    % 1) Load images and parameters
    [LMSCalFormat, ~, ~, ~, Disp, imgParams, pathName] = colorCorrectionGenerateImages( ...
        imgType, setType, m, n, dichromatType, clearFlag);

    % Log some relevant stuff for later
    log(ii).imgType = imgType;
    log(ii).dichromatType = dichromatType;
    log(ii).pathName = pathName;
    log(ii).targetDistToCompare = targetDistToCompare;

    % 2) Build reference images
    rgbLinCalFormat = Disp.M_cones2rgb * LMSCalFormat;

    tol = 1e-3;
    mn = min(rgbLinCalFormat(:));
    mx = max(rgbLinCalFormat(:));

    if mn < -tol || mx > 1+tol
        warning('Clamping: out of range beyond tol (min=%.6g, max=%.6g).', mn, mx);
    end

    rgbLinCalFormat = min(max(rgbLinCalFormat, 0), 1);

    rgbTriOrig = CalFormatToImage(rgbLin2RGB(rgbLinCalFormat, Disp), m, n);

    % Render original for the dichromat and convert to image format
    [~, rgbLinTriRenderedOrig, ~] = DichromRenderLinear(LMSCalFormat, dichromatType, Disp);
    rgbDiOrig = CalFormatToImage(rgbLin2RGB(rgbLinTriRenderedOrig, Disp), m, n);

    % 3) Preallocate images to compare later
    nMetrics = numel(infoFcnList);
    triOutImgs = zeros(m, n, 3, nMetrics);
    diOutImgs  = zeros(m, n, 3, nMetrics);

    % Store curves per metric for the curve plot
    allInfoCurves = cell(nMetrics,1);
    allDistCurves = cell(nMetrics,1);

    % Loop over info metrics
    for kk = 1:nMetrics

        infoFcn    = infoFcnList{kk};
        infoParams = infoParamsList{kk};
        metricName = infoNameList{kk};

        fprintf('\n  ---- Metric %d/%d: %s ----\n', kk, nMetrics, metricName);

        % 4) Create daltonizer object for this metric
        theDaltonizer = daltonize(infoFcn, infoParams, distortionFcn, distortionParams,renderFcn, renderParams, Disp);

        % 5) Run the distortion sweep as normal
        % This does:
        %   for each target distortion:
        %     maximize info subject to (distortion ~= target)
        sweepAxis_dist = 'distortion';

        [LMSSweep_dist, rgbLinSweep_dist, LMSRenderedSweep_dist, rgbLinRenderedSweep_dist, TSweep_dist, targetInfoNorm_dist, targetDistNorm_dist, infoNormAch_dist, distNormAch_dist] = computeSweep( ...
            theDaltonizer, LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_dist);

        % 6) Pick the fixed distortion solution
        % We pick the sweep point whose achieved distortion is closest
        % to targetDistToCompare
        [minAbsErr, idxClosest] = min(abs(distNormAch_dist - targetDistToCompare));

        infoAtTarget = infoNormAch_dist(idxClosest);
        distAtTarget = distNormAch_dist(idxClosest);

        fprintf('     targetDist=%.3f | closest achievedDist=%.3f (|err|=%.3f) | achievedInfo=%.3f\n', ...
            targetDistToCompare, distAtTarget, minAbsErr, infoAtTarget);

        % 7) Extract the image for that solution
        % Optimized *trichromat* image (linear RGB in CalFormat)
        rgbTriCal = rgbLin2RGB(rgbLinSweep_dist{idxClosest}, Disp);
        rgbTriImg = CalFormatToImage(rgbTriCal, m, n);

        % Dichromat-rendered version of that optimized image
        rgbDiCal  = rgbLin2RGB(rgbLinRenderedSweep_dist{idxClosest}, Disp);
        rgbDiImg  = CalFormatToImage(rgbDiCal, m, n);

        % Store for montage
        triOutImgs(:,:,:,kk) = rgbTriImg;
        diOutImgs(:,:,:,kk)  = rgbDiImg;

        % 8) Store data for the curve plots
        log(ii).metric(kk).name           = metricName;
        log(ii).metric(kk).infoFcn        = func2str(infoFcn);
        log(ii).metric(kk).infoParams     = infoParams;

        log(ii).metric(kk).idxClosest     = idxClosest;
        log(ii).metric(kk).distAtTarget   = distAtTarget;
        log(ii).metric(kk).infoAtTarget   = infoAtTarget;
        log(ii).metric(kk).absErrAtTarget = minAbsErr;

        log(ii).metric(kk).T = TSweep_dist{idxClosest};

        log(ii).metric(kk).infoNormAch    = infoNormAch_dist(:);
        log(ii).metric(kk).distNormAch    = distNormAch_dist(:);
        log(ii).metric(kk).targetDistNorm = targetDistNorm_dist(:);

        allInfoCurves{kk} = infoNormAch_dist(:);
        allDistCurves{kk} = distNormAch_dist(:);

    end % metrics loop

    %% %%%%%%%%%%%%%%%
    % Visualizations %
    %%%%%%%%%%%%%%%%%%

    % First graph: images of transformed outputs for each info metric
    %   Row 1 = Trichromat
    %   Row 2 = Dichromat render
    %   Columns = metrics

    nMetrics = numel(infoNameList);
    nCols = nMetrics + 1;
    nRows = 2; % trichromat and dichromat

    % Title of figure and subtitles on each figure
    figName = sprintf('%s — %s — fixed distortion ~%.2f', imgType, dichromatType, targetDistToCompare);
    fig = figure('Color','w','Name',figName);
    tl  = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    % Top row: original trichromat
    ax = nexttile(tl, 1);  % row1 col1
    imshow(rgbTriOrig, 'Parent', ax);
    title(ax, 'Original (tri)', 'Interpreter','none');

    % Bottom row: original dichromat render
    ax = nexttile(tl, 1 + nCols); % row2 col1
    imshow(rgbDiOrig, 'Parent', ax);
    title(ax, sprintf('Original render (%s)', dichromatType), 'Interpreter','none');

    % Metric columns
    for kk = 1:nMetrics
        thisCol = kk + 1; % plus 1 bc of the original image

        % Top row: optimized trichromat output
        ax = nexttile(tl, thisCol);  % row1 col=thisCol
        imshow(triOutImgs(:,:,:,kk), 'Parent', ax);

        title(ax, sprintf('%s | tri\n(dist=%.3f, info=%.3f)', ...
            infoNameList{kk}, log(ii).metric(kk).distAtTarget, log(ii).metric(kk).infoAtTarget), 'Interpreter','none');

        % Bottom row: dichromat render of that optimized output
        ax = nexttile(tl, thisCol + nCols);  % row2 col=thisCol
        imshow(diOutImgs(:,:,:,kk), 'Parent', ax);
        title(ax, sprintf('%s | di render', infoNameList{kk}), 'Interpreter','none');
    end


    % Second graph: achieved info vs achieved distortion curves
    figure('Color','w','Name',sprintf('%s — curves by metric', imgType));
    hold on;
    for kk = 1:nMetrics
        plot(allDistCurves{kk}, allInfoCurves{kk}, '-o','LineWidth', 1.5, 'DisplayName', infoNameList{kk});
    end
    xline(targetDistToCompare, '--', 'Target dist');
    grid on; axis square;
    xlabel('Achieved Distortion (normalized)');
    ylabel('Achieved Info (normalized)');
    title(sprintf('%s — %s — distortion sweep curves', imgType, dichromatType));
    legend('Location','bestoutside');
    hold off;

end % images loop

%% %%%%%%%%%%%%%
% SAVE RESULTS %
%%%%%%%%%%%%%%%%
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');

saveFile = fullfile(outputDir, sprintf('fixedDist%.2f_%s_%dsteps.mat', ...
    targetDistToCompare, dichromatType, nSteps));
save(saveFile, 'log');

fprintf('\nSaved results to:\n  %s\n', saveFile);

k = 1