function t_compareFixedMetrics

% t_compareFixedMetrics.m
% Compare image enhancement outputs across multiple metrics while
% holding either distortion OR info approximately constant.
%
% Minimal change: add a switch:
%   fixedMode = 'dist' or 'info'
% and set:
%   targetDistToCompare OR targetInfoToCompare
%
% If fixedMode='dist': run distortion sweep, pick idx closest in distortion.
% If fixedMode='info': run info sweep, pick idx closest in info.

clear; clc;

%%%%%%%%%%%%%%%%%%%%%
% Images parameters %
%%%%%%%%%%%%%%%%%%%%%
imageTypes = {'flower1.png', 'fruit.png', 'Gaugin.png', 'ishi45.png'};
setType       = 1;
dichromatType = 'Deuteranopia';

% Size of image
m = 67; n = 67;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW: choose whether to fix distortion or fix info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixedMode = 'dist';   % 'dist' or 'info'

% If fixedMode='dist', we match solutions at this achieved distortion:
targetDistToCompare = 0.4;

% If fixedMode='info', we match solutions at this achieved info:
targetInfoToCompare = 0.4; %???? idk what to set this to

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
    log(ii).fixedMode = fixedMode;
    log(ii).targetDistToCompare = targetDistToCompare;
    log(ii).targetInfoToCompare = targetInfoToCompare;

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
        theDaltonizer = daltonize(infoFcn, infoParams, distortionFcn, distortionParams, renderFcn, renderParams, Disp);

        % 5) Run the sweep (distortion or info) based on fixedMode
        if strcmpi(fixedMode, 'dist')
            sweepAxis = 'distortion';
        else
            sweepAxis = 'info';
        end

        [LMSSweep, rgbLinSweep, LMSRenderedSweep, rgbLinRenderedSweep, TSweep, ...
            targetInfoNorm, targetDistNorm, infoNormAch, distNormAch] = computeSweep( ...
            theDaltonizer, LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis);

        % 6) Pick the fixed-(dist or info) solution
        if strcmpi(fixedMode, 'dist')
            [minAbsErr, idxClosest] = min(abs(distNormAch - targetDistToCompare));
            pickedTarget = targetDistToCompare;
            pickedAch    = distNormAch(idxClosest);
            pickedOther  = infoNormAch(idxClosest);
            fprintf('     targetDist=%.3f | closest achievedDist=%.3f (|err|=%.3f) | achievedInfo=%.3f\n', ...
                pickedTarget, pickedAch, minAbsErr, pickedOther);
        else
            [minAbsErr, idxClosest] = min(abs(infoNormAch - targetInfoToCompare));
            pickedTarget = targetInfoToCompare;
            pickedAch    = infoNormAch(idxClosest);
            pickedOther  = distNormAch(idxClosest);
            fprintf('     targetInfo=%.3f | closest achievedInfo=%.3f (|err|=%.3f) | achievedDist=%.3f\n', ...
                pickedTarget, pickedAch, minAbsErr, pickedOther);
        end

        infoAtTarget = infoNormAch(idxClosest);
        distAtTarget = distNormAch(idxClosest);

        % 7) Extract the image for that solution
        % Optimized *trichromat* image (linear RGB in CalFormat)
        rgbTriCal = rgbLin2RGB(rgbLinSweep{idxClosest}, Disp);
        rgbTriImg = CalFormatToImage(rgbTriCal, m, n);

        % Dichromat-rendered version of that optimized image
        rgbDiCal  = rgbLin2RGB(rgbLinRenderedSweep{idxClosest}, Disp);
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

        log(ii).metric(kk).T = TSweep{idxClosest};

        log(ii).metric(kk).infoNormAch    = infoNormAch(:);
        log(ii).metric(kk).distNormAch    = distNormAch(:);
        log(ii).metric(kk).targetInfoNorm = targetInfoNorm(:);
        log(ii).metric(kk).targetDistNorm = targetDistNorm(:);

        allInfoCurves{kk} = infoNormAch(:);
        allDistCurves{kk} = distNormAch(:);

    end % metrics loop

    %% %%%%%%%%%%%%%%%
    % Visualizations %
    %%%%%%%%%%%%%%%%%%

    nMetrics = numel(infoNameList);
    nCols = nMetrics + 1;
    nRows = 2; % trichromat and dichromat

    if strcmpi(fixedMode,'dist')
        figName = sprintf('%s — %s — fixed dist ~%.2f', imgType, dichromatType, targetDistToCompare);
    else
        figName = sprintf('%s — %s — fixed info ~%.2f', imgType, dichromatType, targetInfoToCompare);
    end

    fig = figure('Color','w','Name',figName);
    tl  = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    % Top row: original trichromat
    ax = nexttile(tl, 1);
    imshow(rgbTriOrig, 'Parent', ax);
    title(ax, 'Original (tri)', 'Interpreter','none');

    % Bottom row: original dichromat render
    ax = nexttile(tl, 1 + nCols);
    imshow(rgbDiOrig, 'Parent', ax);
    title(ax, sprintf('Original render (%s)', dichromatType), 'Interpreter','none');

    % Metric columns
    for kk = 1:nMetrics
        thisCol = kk + 1;

        ax = nexttile(tl, thisCol);
        imshow(triOutImgs(:,:,:,kk), 'Parent', ax);

        if strcmpi(fixedMode,'dist')
            title(ax, sprintf('%s | tri\n(dist=%.3f, info=%.3f)', ...
                infoNameList{kk}, log(ii).metric(kk).distAtTarget, log(ii).metric(kk).infoAtTarget), 'Interpreter','none');
        else
            title(ax, sprintf('%s | tri\n(info=%.3f, dist=%.3f)', ...
                infoNameList{kk}, log(ii).metric(kk).infoAtTarget, log(ii).metric(kk).distAtTarget), 'Interpreter','none');
        end

        ax = nexttile(tl, thisCol + nCols);
        imshow(diOutImgs(:,:,:,kk), 'Parent', ax);
        title(ax, sprintf('%s | di render', infoNameList{kk}), 'Interpreter','none');
    end

    % Second graph: achieved info vs achieved distortion curves
    figure('Color','w','Name',sprintf('%s — curves by metric', imgType));
    hold on;
    for kk = 1:nMetrics
        plot(allDistCurves{kk}, allInfoCurves{kk}, '-o','LineWidth', 1.5, 'DisplayName', infoNameList{kk});
    end

    if strcmpi(fixedMode,'dist')
        xline(targetDistToCompare, '--', 'Target dist');
    else
        yline(targetInfoToCompare, '--', 'Target info');
    end

    grid on; axis square;
    xlabel('Achieved Distortion (normalized)');
    ylabel('Achieved Info (normalized)');

    if strcmpi(fixedMode,'dist')
        title(sprintf('%s — %s — distortion sweep curves', imgType, dichromatType));
    else
        title(sprintf('%s — %s — info sweep curves', imgType, dichromatType));
    end

    legend('Location','bestoutside');
    hold off;

end % images loop

%% %%%%%%%%%%%%%
% SAVE RESULTS %
%%%%%%%%%%%%%%%%
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');

if strcmpi(fixedMode,'dist')
    saveFile = fullfile(outputDir, sprintf('fixedDist%.2f_%s_%dsteps.mat', ...
        targetDistToCompare, dichromatType, nSteps));
else
    saveFile = fullfile(outputDir, sprintf('fixedInfo%.2f_%s_%dsteps.mat', ...
        targetInfoToCompare, dichromatType, nSteps));
end

save(saveFile, 'log');
fprintf('\nSaved results to:\n  %s\n', saveFile);

k = 1;
