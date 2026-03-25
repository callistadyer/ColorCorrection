function t_compareFixedMetrics

% t_compareFixedMetrics.m
%
% Compare image-enhancement outputs across multiple information metrics
% while holding either distortion OR information approximately constant.
%
% Core behavior:
%   - Choose whether to fix distortion or fix information
%   - Run a sweep for each information metric
%   - Pick the solution whose achieved value is closest to the desired target
%   - Visualize source-resolution trichromat and dichromat-rendered outputs
%   - Plot achieved info vs achieved distortion curves
%   - Save a log struct summarizing the comparison
%
% Added behavior:
%   - Optionally apply the selected transform T to a higher-resolution
%     version of the same image
%   - Optionally save ONLY those upsampled outputs into:
%         ColorCorrection/compareMethods/<imageName>/<runFolder>/
%
% Example target save folder:
%   compareMethods/fruit/distortion_20steps_128x128/
%
% Notes:
%   - The original source-resolution comparison logic is preserved.
%   - The added target-resolution block reuses the selected 3x3 T matrix
%     without re-optimizing at the higher resolution.

clear;
clc;


%% ========================= USER SETTINGS =========================

%%%%%%%%%%%%%%%%%%%%%
% Image parameters  %
%%%%%%%%%%%%%%%%%%%%%

imageTypes     = {'flower1.png', 'fruit.png', 'Gaugin.png', 'ishi45.png'};
setType        = 1;
dichromatType  = 'Deuteranopia';

% Source image size used for the actual optimization / sweep
m = 61;
n = 61;

% Number of sweep steps
nSteps = 20;

% Clear previously generated image data if desired
clearFlag = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distortion / render metric %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distortion metric
distortionFcn    = @computeDistortion_DE2000;
distortionParams = struct();

% Rendering function used for dichromat visualization
renderFcn    = @DichromRenderLinear;
renderParams = struct();


%%%%%%%%%%%%%%%%%%%%
% Info metric list %
%%%%%%%%%%%%%%%%%%%%

% Information metrics to compare
infoFcnList = { ...
    @computeInfo_regressSquared, ...
    @computeInfo_regressCIE};

infoNameList = { ...
    'regressSquared', ...
    'regressCIE'};

% Alternate metric list kept here for convenience
% infoFcnList = { ...
%     @computeInfo_LMdifference, ...
%     @computeInfo_regress, ...
%     @computeInfo_Wade };
%
% infoNameList = { ...
%     'LMdifference', ...
%     'regress', ...
%     'Wade'};

% Parameter structs for each information metric.
% Some metrics need extra settings; others use an empty struct.
infoParamsList = cell(size(infoFcnList));
for k = 1:numel(infoFcnList)
    if isequal(infoFcnList{k}, @computeInfo_regressCIE) || ...
            isequal(infoFcnList{k}, @computeInfo_regressSquared)

        infoParamsList{k} = struct( ...
            'predictingWhat',     'L,M,S', ...
            'predictingFromWhat', 'L and S');

    else
        infoParamsList{k} = struct();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose whether to fix distortion or fix information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fixedMode = 'dist';   % 'dist' or 'info'

% If fixedMode = 'dist', match solutions at this achieved distortion
targetDistToCompare = 0.4;

% If fixedMode = 'info', match solutions at this achieved information
targetInfoToCompare = 0.4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional: apply the selected transform to a target resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn on/off the higher-resolution re-application block
doApplyToTargetResolution = true;

% Target resolution to which the selected T should be applied
mTarget = 128;
nTarget = 128;

% Save target-resolution outputs
saveTargetResolutionImages = true;

% Show an additional montage figure at the target resolution
showTargetResolutionMontage = true;

% Clamp target-resolution RGB values into [0, 1] before visualization/saving
doClampTargetRGB = true;

% Print simple out-of-gamut stats before clamping at target resolution
printTargetGamutStats = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output directory / bookkeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');

% Folder where upsampled comparison outputs will be saved
compareMethodsDir = fullfile(outputDir, 'compareMethods');

% Create it if it does not already exist
if exist(compareMethodsDir, 'dir') ~= 7
    mkdir(compareMethodsDir);
end


%% ========================= INITIALIZE LOG =========================

log = struct();


%% ========================= MAIN LOOP OVER IMAGES =========================

for ii = 1:numel(imageTypes)

    imgType = imageTypes{ii};

    fprintf('\n====================================================\n');
    fprintf('IMAGE %d/%d: %s | %s\n', ii, numel(imageTypes), imgType, dichromatType);
    fprintf('====================================================\n');

    %% -------------------------------------------------------------
    % 1) Load the source-resolution image and associated parameters
    % --------------------------------------------------------------
    [LMSCalFormat, ~, ~, ~, Disp, imgParams, pathName] = ...
        colorCorrectionGenerateImages( ...
            imgType, setType, m, n, dichromatType, clearFlag);

    % Store image-level metadata in the log
    log(ii).imgType                    = imgType;
    log(ii).dichromatType              = dichromatType;
    log(ii).pathName                   = pathName;
    log(ii).fixedMode                  = fixedMode;
    log(ii).targetDistToCompare        = targetDistToCompare;
    log(ii).targetInfoToCompare        = targetInfoToCompare;
    log(ii).sourceSize                 = [m n];
    log(ii).doApplyToTargetResolution  = doApplyToTargetResolution;
    log(ii).targetSize                 = [mTarget nTarget];

    %% -------------------------------------------------------------
    % 2) Build source-resolution original images
    % --------------------------------------------------------------
    rgbLinCalFormat = Disp.M_cones2rgb * LMSCalFormat;

    tol = 1e-3;
    mn = min(rgbLinCalFormat(:));
    mx = max(rgbLinCalFormat(:));
    if mn < -tol || mx > 1 + tol
        warning('Clamping: out of range beyond tol (min=%.6g, max=%.6g).', mn, mx);
    end
    rgbLinCalFormat = min(max(rgbLinCalFormat, 0), 1);

    % Original trichromat image
    rgbTriOrig = CalFormatToImage(rgbLin2RGB(rgbLinCalFormat, Disp), m, n);

    % Original dichromat-rendered version
    [~, rgbLinTriRenderedOrig, ~] = DichromRenderLinear(LMSCalFormat, dichromatType, Disp);
    rgbDiOrig = CalFormatToImage(rgbLin2RGB(rgbLinTriRenderedOrig, Disp), m, n);

    %% -------------------------------------------------------------
    % 3) Optional: prepare the target-resolution version of this image
    % --------------------------------------------------------------
    if doApplyToTargetResolution

        % Generate the same image at the target resolution
        [LMS_target, ~, ~, ~, DispTarget, imgParamsTarget, ~] = ...
            colorCorrectionGenerateImages( ...
                imgType, setType, mTarget, nTarget, dichromatType, clearFlag);

        % Convert target LMS to linear RGB so that T can later be applied
        % in RGB contrast space
        RGBCalFormat_target = DispTarget.M_cones2rgb * LMS_target;
        RGBContrast_target  = (RGBCalFormat_target - DispTarget.grayRGB) ./ DispTarget.grayRGB;

        % Target LMS contrast is needed for recomputing info metrics
        LMSContrast_target = (LMS_target - DispTarget.grayLMS) ./ DispTarget.grayLMS;

        % Build target-resolution original trichromat image
        rgbLinCalFormat_target_orig = RGBCalFormat_target;

        mnTarget = min(rgbLinCalFormat_target_orig(:));
        mxTarget = max(rgbLinCalFormat_target_orig(:));
        if mnTarget < -tol || mxTarget > 1 + tol
            warning('Target clamp: out of range beyond tol (min=%.6g, max=%.6g).', ...
                mnTarget, mxTarget);
        end
        rgbLinCalFormat_target_orig = min(max(rgbLinCalFormat_target_orig, 0), 1);

        rgbTriOrigTarget = CalFormatToImage( ...
            rgbLin2RGB(rgbLinCalFormat_target_orig, DispTarget), mTarget, nTarget);

        % Build target-resolution original dichromat-rendered image
        [~, rgbLinTriRenderedOrigTarget, ~] = ...
            DichromRenderLinear(LMS_target, dichromatType, DispTarget);

        rgbDiOrigTarget = CalFormatToImage( ...
            rgbLin2RGB(rgbLinTriRenderedOrigTarget, DispTarget), mTarget, nTarget);

        % Build target-resolution reference images for normalized metrics
        [protLMS_target, ~, ~] = DichromRenderLinear(LMS_target, 'Protanopia', DispTarget);
        [deutLMS_target, ~, ~] = DichromRenderLinear(LMS_target, 'Deuteranopia', DispTarget);
        [tritLMS_target, ~, ~] = DichromRenderLinear(LMS_target, 'Tritanopia', DispTarget);

        LMS_ref_target = [ ...
            protLMS_target(1, :); ...
            deutLMS_target(2, :); ...
            tritLMS_target(3, :)];

        LMSContrast_ref_target = (LMS_ref_target - DispTarget.grayLMS) ./ DispTarget.grayLMS;

        distortionNormalizerTarget = distortionFcn( ...
            LMS_target, ...
            LMS_ref_target, ...
            imgParamsTarget, ...
            1, ...
            DispTarget, ...
            distortionParams);

        % Preallocate target-resolution images for the optional montage
        nMetrics = numel(infoFcnList);
        triOutImgsTarget = zeros(mTarget, nTarget, 3, nMetrics);
        diOutImgsTarget  = zeros(mTarget, nTarget, 3, nMetrics);

    end


    %% -------------------------------------------------------------
    % 4) Preallocate source-resolution outputs
    % --------------------------------------------------------------
    nMetrics = numel(infoFcnList);

    % Source-resolution selected images for each metric
    triOutImgs = zeros(m, n, 3, nMetrics);
    diOutImgs  = zeros(m, n, 3, nMetrics);

    % Curves for the info-vs-distortion plot
    allInfoCurves = cell(nMetrics, 1);
    allDistCurves = cell(nMetrics, 1);


    %% ========================= LOOP OVER INFO METRICS =========================
    for kk = 1:nMetrics

        infoFcn    = infoFcnList{kk};
        infoParams = infoParamsList{kk};
        metricName = infoNameList{kk};

        fprintf('\n  ---- Metric %d/%d: %s ----\n', kk, nMetrics, metricName);

        %% ---------------------------------------------------------
        % 5) Create a daltonizer object for this metric
        % ----------------------------------------------------------
        theDaltonizer = daltonize( ...
            infoFcn, infoParams, ...
            distortionFcn, distortionParams, ...
            renderFcn, renderParams, ...
            Disp);

        %% ---------------------------------------------------------
        % 6) Choose which sweep to run based on fixedMode
        % ----------------------------------------------------------
        if strcmpi(fixedMode, 'dist')
            sweepAxis = 'distortion';
        else
            sweepAxis = 'info';
        end

        %% ---------------------------------------------------------
        % 7) Run the sweep at the source resolution
        % ----------------------------------------------------------
        [LMSSweep, rgbLinSweep, LMSRenderedSweep, rgbLinRenderedSweep, TSweep, ...
            targetInfoNorm, targetDistNorm, infoNormAch, distNormAch] = computeSweep( ...
            theDaltonizer, LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis);

        %#ok<NASGU> % LMSSweep and LMSRenderedSweep intentionally retained

        %% ---------------------------------------------------------
        % 8) Pick the solution closest to the requested fixed target
        % ----------------------------------------------------------
        if strcmpi(fixedMode, 'dist')

            [minAbsErr, idxClosest] = min(abs(distNormAch - targetDistToCompare));
            pickedTarget = targetDistToCompare;
            pickedAch    = distNormAch(idxClosest);
            pickedOther  = infoNormAch(idxClosest);

            fprintf( ...
                '     targetDist=%.3f | closest achievedDist=%.3f (|err|=%.3f) | achievedInfo=%.3f\n', ...
                pickedTarget, pickedAch, minAbsErr, pickedOther);

        else

            [minAbsErr, idxClosest] = min(abs(infoNormAch - targetInfoToCompare));
            pickedTarget = targetInfoToCompare;
            pickedAch    = infoNormAch(idxClosest);
            pickedOther  = distNormAch(idxClosest);

            fprintf( ...
                '     targetInfo=%.3f | closest achievedInfo=%.3f (|err|=%.3f) | achievedDist=%.3f\n', ...
                pickedTarget, pickedAch, minAbsErr, pickedOther);

        end

        infoAtTarget = infoNormAch(idxClosest);
        distAtTarget = distNormAch(idxClosest);

        %% ---------------------------------------------------------
        % 9) Extract the selected source-resolution images
        % ----------------------------------------------------------
        % Optimized trichromat image
        rgbTriCal = rgbLin2RGB(rgbLinSweep{idxClosest}, Disp);
        rgbTriImg = CalFormatToImage(rgbTriCal, m, n);

        % Dichromat-rendered version of that optimized image
        rgbDiCal = rgbLin2RGB(rgbLinRenderedSweep{idxClosest}, Disp);
        rgbDiImg = CalFormatToImage(rgbDiCal, m, n);

        % Store for source-resolution montage
        triOutImgs(:, :, :, kk) = rgbTriImg;
        diOutImgs(:, :, :, kk)  = rgbDiImg;

        %% ---------------------------------------------------------
        % 10) Store source-resolution curve/log data
        % ----------------------------------------------------------
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

        %% ---------------------------------------------------------
        % 11) Optional: apply the selected T to the target resolution
        % ----------------------------------------------------------
        if doApplyToTargetResolution

            % Selected 3x3 transform from the source-resolution sweep
            T = TSweep{idxClosest};

            % Metric-specific info normalizer at the target resolution
            infoNormalizerTarget = infoFcn( ...
                LMSContrast_target, ...
                LMSContrast_ref_target, ...
                imgParamsTarget, ...
                dichromatType, ...
                1, ...
                DispTarget, ...
                infoParams);

            % Apply T in RGB contrast space at the target resolution
            newRGBContrastTarget  = T * RGBContrast_target;
            newRGBCalFormatTarget = newRGBContrastTarget .* DispTarget.grayRGB + DispTarget.grayRGB;

            % Optional gamut diagnostics before clamping
            if printTargetGamutStats
                fprintf( ...
                    '     target %dx%d | step %02d | min=%g max=%g frac<0=%.3f frac>1=%.3f\n', ...
                    mTarget, nTarget, ...
                    idxClosest, ...
                    min(newRGBCalFormatTarget, [], 'all'), ...
                    max(newRGBCalFormatTarget, [], 'all'), ...
                    mean(newRGBCalFormatTarget(:) < 0), ...
                    mean(newRGBCalFormatTarget(:) > 1));
            end

            % Optional clamp before visualization / saving
            if doClampTargetRGB
                newRGBCalFormatTarget(newRGBCalFormatTarget < 0) = 0;
                newRGBCalFormatTarget(newRGBCalFormatTarget > 1) = 1;
            end

            % Convert back to LMS for target-resolution metrics/rendering
            newLMSTarget         = DispTarget.M_rgb2cones * newRGBCalFormatTarget;
            newLMSContrastTarget = (newLMSTarget - DispTarget.grayLMS) ./ DispTarget.grayLMS;

            % Render the transformed target-resolution image for the dichromat
            [~, rgbLinRenderedTarget, ~] = renderFcn(newLMSTarget, dichromatType, DispTarget);

            % Convert transformed trichromat image to display RGB
            rgbTriTargetCal = rgbLin2RGB(newRGBCalFormatTarget, DispTarget);
            rgbTriTargetImg = CalFormatToImage(rgbTriTargetCal, mTarget, nTarget);

            % Convert transformed dichromat-rendered image to display RGB
            rgbDiTargetCal = rgbLin2RGB(rgbLinRenderedTarget, DispTarget);
            rgbDiTargetImg = CalFormatToImage(rgbDiTargetCal, mTarget, nTarget);

            % Store for optional target-resolution montage
            triOutImgsTarget(:, :, :, kk) = rgbTriTargetImg;
            diOutImgsTarget(:, :, :, kk)  = rgbDiTargetImg;

            % Recompute normalized achieved values at the target resolution
            [~, infoAtTargetResolution] = infoFcn( ...
                LMSContrast_target, ...
                newLMSContrastTarget, ...
                imgParamsTarget, ...
                dichromatType, ...
                infoNormalizerTarget, ...
                DispTarget, ...
                infoParams);

            [~, distAtTargetResolution] = distortionFcn( ...
                LMS_target, ...
                newLMSTarget, ...
                imgParamsTarget, ...
                distortionNormalizerTarget, ...
                DispTarget, ...
                distortionParams);

            % Log target-resolution metadata/results
            log(ii).metric(kk).targetResolution.applied = true;
            log(ii).metric(kk).targetResolution.size    = [mTarget nTarget];
            log(ii).metric(kk).targetResolution.infoAch = infoAtTargetResolution;
            log(ii).metric(kk).targetResolution.distAch = distAtTargetResolution;
            log(ii).metric(kk).targetResolution.doClamp = doClampTargetRGB;

            % Save ONLY the target-resolution images into compareMethods
            if saveTargetResolutionImages

                % Remove file extension from image name for cleaner folder names
                [~, imgBaseName, ~] = fileparts(imgType);

                % Example:
                %   compareMethods/fruit/distortion_20steps_128x128
                compareRunFolder = sprintf('%s_%dsteps_%dx%d', ...
                    lower(char(sweepAxis)), nSteps, mTarget, nTarget);

                saveSubdirTarget = fullfile(compareMethodsDir, imgBaseName, compareRunFolder);

                if exist(saveSubdirTarget, 'dir') ~= 7
                    mkdir(saveSubdirTarget);
                end

                % Include fixed target + selected step in the filename
                if strcmpi(fixedMode, 'dist')
                    targetTag = sprintf('fixedDist_%0.3f', targetDistToCompare);
                else
                    targetTag = sprintf('fixedInfo_%0.3f', targetInfoToCompare);
                end
                targetTag = strrep(targetTag, '.', 'p');

                fileStem = sprintf('%s_%s_step_%03d', metricName, targetTag, idxClosest);

                triToWrite = max(min(rgbTriTargetImg, 1), 0);
                diToWrite  = max(min(rgbDiTargetImg, 1), 0);

                triPath = fullfile(saveSubdirTarget, ...
                    sprintf('%s_trichromat.png', fileStem));

                diPath = fullfile(saveSubdirTarget, ...
                    sprintf('%s_%s.png', fileStem, lower(char(dichromatType))));

                imwrite(triToWrite, triPath);
                imwrite(diToWrite, diPath);

                log(ii).metric(kk).targetResolution.saveFolder = saveSubdirTarget;
                log(ii).metric(kk).targetResolution.triPath    = triPath;
                log(ii).metric(kk).targetResolution.diPath     = diPath;

            else
                log(ii).metric(kk).targetResolution.saveFolder = '';
                log(ii).metric(kk).targetResolution.triPath    = '';
                log(ii).metric(kk).targetResolution.diPath     = '';
            end

        else
            log(ii).metric(kk).targetResolution.applied = false;
        end

    end % metric loop


    %% ========================= SOURCE-RESOLUTION VISUALIZATIONS =========================

    nMetrics = numel(infoNameList);
    nCols    = nMetrics + 1;
    nRows    = 2;   % trichromat row + dichromat-render row

    if strcmpi(fixedMode, 'dist')
        figName = sprintf('%s — %s — fixed dist ~%.2f', ...
            imgType, dichromatType, targetDistToCompare);
    else
        figName = sprintf('%s — %s — fixed info ~%.2f', ...
            imgType, dichromatType, targetInfoToCompare);
    end

    fig = figure('Color', 'w', 'Name', figName);
    tl  = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Top-left: original trichromat
    ax = nexttile(tl, 1);
    imshow(rgbTriOrig, 'Parent', ax);
    title(ax, 'Original (tri)', 'Interpreter', 'none');

    % Bottom-left: original dichromat render
    ax = nexttile(tl, 1 + nCols);
    imshow(rgbDiOrig, 'Parent', ax);
    title(ax, sprintf('Original render (%s)', dichromatType), 'Interpreter', 'none');

    % One column per metric
    for kk = 1:nMetrics
        thisCol = kk + 1;

        % Top row: trichromat optimized image
        ax = nexttile(tl, thisCol);
        imshow(triOutImgs(:, :, :, kk), 'Parent', ax);

        if strcmpi(fixedMode, 'dist')
            title(ax, sprintf('%s | tri\n(dist=%.3f, info=%.3f)', ...
                infoNameList{kk}, ...
                log(ii).metric(kk).distAtTarget, ...
                log(ii).metric(kk).infoAtTarget), ...
                'Interpreter', 'none');
        else
            title(ax, sprintf('%s | tri\n(info=%.3f, dist=%.3f)', ...
                infoNameList{kk}, ...
                log(ii).metric(kk).infoAtTarget, ...
                log(ii).metric(kk).distAtTarget), ...
                'Interpreter', 'none');
        end

        % Bottom row: dichromat-rendered optimized image
        ax = nexttile(tl, thisCol + nCols);
        imshow(diOutImgs(:, :, :, kk), 'Parent', ax);
        title(ax, sprintf('%s | di render', infoNameList{kk}), 'Interpreter', 'none');
    end


    %% ========================= CURVE PLOT =========================

    figure('Color', 'w', 'Name', sprintf('%s — curves by metric', imgType));
    hold on;

    for kk = 1:nMetrics
        plot(allDistCurves{kk}, allInfoCurves{kk}, '-o', ...
            'LineWidth', 1.5, ...
            'DisplayName', infoNameList{kk});
    end

    if strcmpi(fixedMode, 'dist')
        xline(targetDistToCompare, '--', 'Target dist');
    else
        yline(targetInfoToCompare, '--', 'Target info');
    end

    grid on;
    axis square;
    xlabel('Achieved Distortion (normalized)');
    ylabel('Achieved Info (normalized)');

    if strcmpi(fixedMode, 'dist')
        title(sprintf('%s — %s — distortion sweep curves', imgType, dichromatType));
    else
        title(sprintf('%s — %s — info sweep curves', imgType, dichromatType));
    end

    legend('Location', 'bestoutside');
    hold off;


    %% ========================= TARGET-RESOLUTION MONTAGE =========================

    if doApplyToTargetResolution && showTargetResolutionMontage

        if strcmpi(fixedMode, 'dist')
            figNameTarget = sprintf('%s — %s — fixed dist ~%.2f — %dx%d', ...
                imgType, dichromatType, targetDistToCompare, mTarget, nTarget);
        else
            figNameTarget = sprintf('%s — %s — fixed info ~%.2f — %dx%d', ...
                imgType, dichromatType, targetInfoToCompare, mTarget, nTarget);
        end

        figTarget = figure('Color', 'w', 'Name', figNameTarget);
        tlTarget  = tiledlayout(figTarget, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Top-left: target-resolution original trichromat
        ax = nexttile(tlTarget, 1);
        imshow(rgbTriOrigTarget, 'Parent', ax);
        title(ax, sprintf('Original (tri) %dx%d', mTarget, nTarget), 'Interpreter', 'none');

        % Bottom-left: target-resolution original dichromat render
        ax = nexttile(tlTarget, 1 + nCols);
        imshow(rgbDiOrigTarget, 'Parent', ax);
        title(ax, sprintf('Original render (%s)', dichromatType), 'Interpreter', 'none');

        % One column per metric
        for kk = 1:nMetrics
            thisCol = kk + 1;

            % Top row: transformed trichromat image at target resolution
            ax = nexttile(tlTarget, thisCol);
            imshow(triOutImgsTarget(:, :, :, kk), 'Parent', ax);

            if strcmpi(fixedMode, 'dist')
                title(ax, sprintf('%s | tri\n(dist=%.3f, info=%.3f)', ...
                    infoNameList{kk}, ...
                    log(ii).metric(kk).targetResolution.distAch, ...
                    log(ii).metric(kk).targetResolution.infoAch), ...
                    'Interpreter', 'none');
            else
                title(ax, sprintf('%s | tri\n(info=%.3f, dist=%.3f)', ...
                    infoNameList{kk}, ...
                    log(ii).metric(kk).targetResolution.infoAch, ...
                    log(ii).metric(kk).targetResolution.distAch), ...
                    'Interpreter', 'none');
            end

            % Bottom row: dichromat-rendered transformed image at target resolution
            ax = nexttile(tlTarget, thisCol + nCols);
            imshow(diOutImgsTarget(:, :, :, kk), 'Parent', ax);
            title(ax, sprintf('%s | di render', infoNameList{kk}), 'Interpreter', 'none');
        end
    end

end % image loop


%% ========================= SAVE RESULTS =========================

% Save the summary log struct
if strcmpi(fixedMode, 'dist')
    saveFile = fullfile(outputDir, sprintf('fixedDist%.2f_%s_%dsteps.mat', ...
        targetDistToCompare, dichromatType, nSteps));
else
    saveFile = fullfile(outputDir, sprintf('fixedInfo%.2f_%s_%dsteps.mat', ...
        targetInfoToCompare, dichromatType, nSteps));
end

save(saveFile, 'log');
fprintf('\nSaved results to:\n  %s\n', saveFile);

% Original trailing line kept as-is
k = 1;

% function t_compareFixedMetrics
% 
% % t_compareFixedMetrics.m
% %
% % Compare image-enhancement outputs across multiple information metrics
% % while holding either distortion OR information approximately constant.
% %
% % Original core behavior:
% %   - Choose whether to fix distortion or info using fixedMode
% %   - Run the appropriate sweep for each metric
% %   - Pick the solution whose achieved value is closest to the desired target
% %   - Visualize source-resolution trichromat and dichromat-rendered outputs
% %   - Plot achieved info vs achieved distortion curves
% %   - Save a log struct to disk
% %
% % Added capability:
% %   - Optionally apply the selected transform T to a higher-resolution
% %     version of the same image
% %   - Optionally save those higher-resolution outputs into a folder located
% %     near the sweep results, e.g.
% %         distortion_20steps_128x128
% %         info_20steps_256x256
% %
% % Important:
% %   - The original source-resolution comparison logic is unchanged.
% %   - The target-resolution application is optional and controlled by flags
% %     in the USER SETTINGS section below.
% 
% clear;
% clc;
% 
% 
% %% ========================= USER SETTINGS =========================
% 
% %%%%%%%%%%%%%%%%%%%%%
% % Image parameters  %
% %%%%%%%%%%%%%%%%%%%%%
% 
% imageTypes     = {'flower1.png', 'fruit.png', 'Gaugin.png', 'ishi45.png'};
% setType        = 1;
% dichromatType  = 'Deuteranopia';
% 
% % Source image size used for the actual optimization / sweep
% m = 61;
% n = 61;
% 
% % Number of sweep steps
% nSteps = 20;
% 
% % Clear previously generated image data if desired
% clearFlag = 0;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Distortion / render metric %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Distortion metric
% distortionFcn    = @computeDistortion_DE2000;
% distortionParams = struct();
% 
% % Rendering function used for dichromat visualization
% renderFcn    = @DichromRenderLinear;
% renderParams = struct();
% 
% 
% %%%%%%%%%%%%%%%%%%%%
% % Info metric list %
% %%%%%%%%%%%%%%%%%%%%
% 
% % Information metrics to compare
% infoFcnList = { ...
%     @computeInfo_regressSquared, ...
%     @computeInfo_regressCIE};
% 
% infoNameList = { ...
%     'regressSquared', ...
%     'regressCIE'};
% 
% % Alternate metric list kept here for convenience
% % infoFcnList = {@computeInfo_LMdifference, ...
% %     @computeInfo_regress, ...
% %     @computeInfo_Wade };
% %
% % infoNameList = { ...
% %     'LMdifference', ...
% %     'regress', ...
% %     'Wade' };
% 
% % Parameter structs for each information metric.
% % Some metrics need extra settings; others use an empty struct.
% infoParamsList = cell(size(infoFcnList));
% for k = 1:numel(infoFcnList)
%     if isequal(infoFcnList{k}, @computeInfo_regressCIE) || ...
%             isequal(infoFcnList{k}, @computeInfo_regressSquared)
% 
%         infoParamsList{k} = struct( ...
%             'predictingWhat',     'L,M,S', ...
%             'predictingFromWhat', 'L and S');
% 
%     else
%         infoParamsList{k} = struct();
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Choose whether to fix distortion or fix information
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fixedMode = 'dist';   % 'dist' or 'info'
% 
% % If fixedMode = 'dist', match solutions at this achieved distortion
% targetDistToCompare = 0.4;
% 
% % If fixedMode = 'info', match solutions at this achieved information
% targetInfoToCompare = 0.4;   % choose based on the typical achieved range
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Optional: apply the selected transform to a target resolution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Turn on/off the higher-resolution re-application block
% doApplyToTargetResolution = true;
% 
% % Target resolution to which the selected T should be applied
% mTarget = 128*2;
% nTarget = 128*2;
% 
% % Save the target-resolution outputs near the sweep folders
% saveTargetResolutionImages = true;
% 
% % Show an additional montage figure at the target resolution
% showTargetResolutionMontage = true;
% 
% % Clamp target-resolution RGB values into [0, 1] before visualization/saving
% doClampTargetRGB = true;
% 
% % Print simple out-of-gamut stats before clamping at target resolution
% printTargetGamutStats = true;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Output directory / bookkeeping
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% projectName = 'ColorCorrection';
% outputDir   = getpref(projectName, 'outputDir');
% 
% % Root directory used by the sweep code for transformed images
% saveBase = fullfile(outputDir, 'testImagesTransformed');
% 
% 
% %% ========================= INITIALIZE LOG =========================
% 
% log = struct();
% 
% 
% %% ========================= MAIN LOOP OVER IMAGES =========================
% 
% for ii = 1:numel(imageTypes)
% 
%     imgType = imageTypes{ii};
% 
%     fprintf('\n====================================================\n');
%     fprintf('IMAGE %d/%d: %s | %s\n', ii, numel(imageTypes), imgType, dichromatType);
%     fprintf('====================================================\n');
% 
%     %% -------------------------------------------------------------
%     % 1) Load the source-resolution image and parameters
%     % --------------------------------------------------------------
%     [LMSCalFormat, ~, ~, ~, Disp, imgParams, pathName] = colorCorrectionGenerateImages( ...
%         imgType, setType, m, n, dichromatType, clearFlag);
% 
%     % Store image-level metadata in the log
%     log(ii).imgType               = imgType;
%     log(ii).dichromatType         = dichromatType;
%     log(ii).pathName              = pathName;
%     log(ii).fixedMode             = fixedMode;
%     log(ii).targetDistToCompare   = targetDistToCompare;
%     log(ii).targetInfoToCompare   = targetInfoToCompare;
%     log(ii).sourceSize            = [m n];
%     log(ii).doApplyToTargetResolution = doApplyToTargetResolution;
%     log(ii).targetSize            = [mTarget nTarget];
% 
%     %% -------------------------------------------------------------
%     % 2) Build source-resolution reference/original images
%     % --------------------------------------------------------------
%     rgbLinCalFormat = Disp.M_cones2rgb * LMSCalFormat;
% 
%     tol = 1e-3;
%     mn = min(rgbLinCalFormat(:));
%     mx = max(rgbLinCalFormat(:));
%     if mn < -tol || mx > 1 + tol
%         warning('Clamping: out of range beyond tol (min=%.6g, max=%.6g).', mn, mx);
%     end
%     rgbLinCalFormat = min(max(rgbLinCalFormat, 0), 1);
% 
%     % Original trichromat image
%     rgbTriOrig = CalFormatToImage(rgbLin2RGB(rgbLinCalFormat, Disp), m, n);
% 
%     % Original dichromat-rendered version
%     [~, rgbLinTriRenderedOrig, ~] = DichromRenderLinear(LMSCalFormat, dichromatType, Disp);
%     rgbDiOrig = CalFormatToImage(rgbLin2RGB(rgbLinTriRenderedOrig, Disp), m, n);
% 
%     %% -------------------------------------------------------------
%     % 3) Optional: prepare the target-resolution version of this image
%     % --------------------------------------------------------------
%     if doApplyToTargetResolution
% 
%         % Generate the same image at the target resolution
%         [LMS_target, ~, ~, ~, DispTarget, imgParamsTarget, ~] = colorCorrectionGenerateImages( ...
%             imgType, setType, mTarget, nTarget, dichromatType, clearFlag);
% 
%         % Convert target LMS to linear RGB so we can later apply T in RGB contrast
%         RGBCalFormat_target = DispTarget.M_cones2rgb * LMS_target;
%         RGBContrast_target  = (RGBCalFormat_target - DispTarget.grayRGB) ./ DispTarget.grayRGB;
% 
%         % Target LMS contrast is needed for recomputing target-resolution info
%         LMSContrast_target = (LMS_target - DispTarget.grayLMS) ./ DispTarget.grayLMS;
% 
%         % Build target-resolution original trichromat image
%         rgbLinCalFormat_target_orig = RGBCalFormat_target;
% 
%         mnTarget = min(rgbLinCalFormat_target_orig(:));
%         mxTarget = max(rgbLinCalFormat_target_orig(:));
%         if mnTarget < -tol || mxTarget > 1 + tol
%             warning('Target clamp: out of range beyond tol (min=%.6g, max=%.6g).', mnTarget, mxTarget);
%         end
%         rgbLinCalFormat_target_orig = min(max(rgbLinCalFormat_target_orig, 0), 1);
% 
%         rgbTriOrigTarget = CalFormatToImage( ...
%             rgbLin2RGB(rgbLinCalFormat_target_orig, DispTarget), mTarget, nTarget);
% 
%         % Build target-resolution original dichromat-rendered image
%         [~, rgbLinTriRenderedOrigTarget, ~] = DichromRenderLinear(LMS_target, dichromatType, DispTarget);
%         rgbDiOrigTarget = CalFormatToImage( ...
%             rgbLin2RGB(rgbLinTriRenderedOrigTarget, DispTarget), mTarget, nTarget);
% 
%         % Precompute the target-resolution reference image used by the
%         % normalized metrics. The distortion normalizer is shared across
%         % metrics; the info normalizer is metric-specific and will be
%         % computed inside the metric loop.
%         [protLMS_target, ~, ~] = DichromRenderLinear(LMS_target, 'Protanopia', DispTarget);
%         [deutLMS_target, ~, ~] = DichromRenderLinear(LMS_target, 'Deuteranopia', DispTarget);
%         [tritLMS_target, ~, ~] = DichromRenderLinear(LMS_target, 'Tritanopia', DispTarget);
% 
%         LMS_ref_target = [ ...
%             protLMS_target(1, :); ...
%             deutLMS_target(2, :); ...
%             tritLMS_target(3, :)];
% 
%         LMSContrast_ref_target = (LMS_ref_target - DispTarget.grayLMS) ./ DispTarget.grayLMS;
% 
%         distortionNormalizerTarget = distortionFcn( ...
%             LMS_target, ...
%             LMS_ref_target, ...
%             imgParamsTarget, ...
%             1, ...
%             DispTarget, ...
%             distortionParams);
% 
%         % Preallocate target-resolution output arrays for an optional montage
%         nMetrics = numel(infoFcnList);
%         triOutImgsTarget = zeros(mTarget, nTarget, 3, nMetrics);
%         diOutImgsTarget  = zeros(mTarget, nTarget, 3, nMetrics);
% 
%     end
% 
% 
%     %% -------------------------------------------------------------
%     % 4) Preallocate source-resolution outputs
%     % --------------------------------------------------------------
%     nMetrics = numel(infoFcnList);
% 
%     % Source-resolution images selected for each metric
%     triOutImgs = zeros(m, n, 3, nMetrics);
%     diOutImgs  = zeros(m, n, 3, nMetrics);
% 
%     % Curves for the info-vs-distortion plot
%     allInfoCurves = cell(nMetrics, 1);
%     allDistCurves = cell(nMetrics, 1);
% 
% 
%     %% ========================= LOOP OVER INFO METRICS =========================
%     for kk = 1:nMetrics
% 
%         infoFcn    = infoFcnList{kk};
%         infoParams = infoParamsList{kk};
%         metricName = infoNameList{kk};
% 
%         fprintf('\n  ---- Metric %d/%d: %s ----\n', kk, nMetrics, metricName);
% 
%         %% ---------------------------------------------------------
%         % 5) Create a daltonizer object for this metric
%         % ----------------------------------------------------------
%         theDaltonizer = daltonize( ...
%             infoFcn, infoParams, ...
%             distortionFcn, distortionParams, ...
%             renderFcn, renderParams, ...
%             Disp);
% 
%         %% ---------------------------------------------------------
%         % 6) Choose which sweep to run based on fixedMode
%         % ----------------------------------------------------------
%         if strcmpi(fixedMode, 'dist')
%             sweepAxis = 'distortion';
%         else
%             sweepAxis = 'info';
%         end
% 
%         %% ---------------------------------------------------------
%         % 7) Run the sweep at the source resolution
%         % ----------------------------------------------------------
%         [LMSSweep, rgbLinSweep, LMSRenderedSweep, rgbLinRenderedSweep, TSweep, ...
%             targetInfoNorm, targetDistNorm, infoNormAch, distNormAch] = computeSweep( ...
%             theDaltonizer, LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis);
% 
%         %#ok<NASGU> % LMSSweep and LMSRenderedSweep are intentionally retained
% 
%         %% ---------------------------------------------------------
%         % 8) Pick the solution closest to the requested fixed target
%         % ----------------------------------------------------------
%         if strcmpi(fixedMode, 'dist')
% 
%             [minAbsErr, idxClosest] = min(abs(distNormAch - targetDistToCompare));
%             pickedTarget = targetDistToCompare;
%             pickedAch    = distNormAch(idxClosest);
%             pickedOther  = infoNormAch(idxClosest);
% 
%             fprintf( ...
%                 '     targetDist=%.3f | closest achievedDist=%.3f (|err|=%.3f) | achievedInfo=%.3f\n', ...
%                 pickedTarget, pickedAch, minAbsErr, pickedOther);
% 
%         else
% 
%             [minAbsErr, idxClosest] = min(abs(infoNormAch - targetInfoToCompare));
%             pickedTarget = targetInfoToCompare;
%             pickedAch    = infoNormAch(idxClosest);
%             pickedOther  = distNormAch(idxClosest);
% 
%             fprintf( ...
%                 '     targetInfo=%.3f | closest achievedInfo=%.3f (|err|=%.3f) | achievedDist=%.3f\n', ...
%                 pickedTarget, pickedAch, minAbsErr, pickedOther);
% 
%         end
% 
%         infoAtTarget = infoNormAch(idxClosest);
%         distAtTarget = distNormAch(idxClosest);
% 
%         %% ---------------------------------------------------------
%         % 9) Extract the selected source-resolution images
%         % ----------------------------------------------------------
%         % Optimized trichromat image (linear RGB in CalFormat -> display RGB)
%         rgbTriCal = rgbLin2RGB(rgbLinSweep{idxClosest}, Disp);
%         rgbTriImg = CalFormatToImage(rgbTriCal, m, n);
% 
%         % Dichromat-rendered version of that optimized image
%         rgbDiCal = rgbLin2RGB(rgbLinRenderedSweep{idxClosest}, Disp);
%         rgbDiImg = CalFormatToImage(rgbDiCal, m, n);
% 
%         % Store for the source-resolution montage
%         triOutImgs(:, :, :, kk) = rgbTriImg;
%         diOutImgs(:, :, :, kk)  = rgbDiImg;
% 
%         %% ---------------------------------------------------------
%         % 10) Store source-resolution curve/log data
%         % ----------------------------------------------------------
%         log(ii).metric(kk).name           = metricName;
%         log(ii).metric(kk).infoFcn        = func2str(infoFcn);
%         log(ii).metric(kk).infoParams     = infoParams;
% 
%         log(ii).metric(kk).idxClosest     = idxClosest;
%         log(ii).metric(kk).distAtTarget   = distAtTarget;
%         log(ii).metric(kk).infoAtTarget   = infoAtTarget;
%         log(ii).metric(kk).absErrAtTarget = minAbsErr;
% 
%         log(ii).metric(kk).T = TSweep{idxClosest};
% 
%         log(ii).metric(kk).infoNormAch    = infoNormAch(:);
%         log(ii).metric(kk).distNormAch    = distNormAch(:);
%         log(ii).metric(kk).targetInfoNorm = targetInfoNorm(:);
%         log(ii).metric(kk).targetDistNorm = targetDistNorm(:);
% 
%         allInfoCurves{kk} = infoNormAch(:);
%         allDistCurves{kk} = distNormAch(:);
% 
%         %% ---------------------------------------------------------
%         % 11) Optional: apply the selected T to the target resolution
%         % ----------------------------------------------------------
%         if doApplyToTargetResolution
% 
%             % The selected 3x3 transform from the source-resolution sweep
%             T = TSweep{idxClosest};
% 
%             % Metric-specific info normalizer at the target resolution
%             infoNormalizerTarget = infoFcn( ...
%                 LMSContrast_target, ...
%                 LMSContrast_ref_target, ...
%                 imgParamsTarget, ...
%                 dichromatType, ...
%                 1, ...
%                 DispTarget, ...
%                 infoParams);
% 
%             % Apply T in RGB contrast space at the target resolution
%             newRGBContrastTarget  = T * RGBContrast_target;
%             newRGBCalFormatTarget = newRGBContrastTarget .* DispTarget.grayRGB + DispTarget.grayRGB;
% 
%             % Optional gamut diagnostics before clamping
%             if printTargetGamutStats
%                 fprintf( ...
%                     '     target %dx%d | step %02d | min=%g max=%g frac<0=%.3f frac>1=%.3f\n', ...
%                     mTarget, nTarget, ...
%                     idxClosest, ...
%                     min(newRGBCalFormatTarget, [], 'all'), ...
%                     max(newRGBCalFormatTarget, [], 'all'), ...
%                     mean(newRGBCalFormatTarget(:) < 0), ...
%                     mean(newRGBCalFormatTarget(:) > 1));
%             end
% 
%             % Optional clamp before visualization / saving
%             if doClampTargetRGB
%                 newRGBCalFormatTarget(newRGBCalFormatTarget < 0) = 0;
%                 newRGBCalFormatTarget(newRGBCalFormatTarget > 1) = 1;
%             end
% 
%             % Convert back to LMS for target-resolution metrics/rendering
%             newLMSTarget         = DispTarget.M_rgb2cones * newRGBCalFormatTarget;
%             newLMSContrastTarget = (newLMSTarget - DispTarget.grayLMS) ./ DispTarget.grayLMS;
% 
%             % Render the target-resolution transformed image for the dichromat
%             [~, rgbLinRenderedTarget, ~] = renderFcn(newLMSTarget, dichromatType, DispTarget);
% 
%             % Convert transformed trichromat image to display RGB
%             rgbTriTargetCal = rgbLin2RGB(newRGBCalFormatTarget, DispTarget);
%             rgbTriTargetImg = CalFormatToImage(rgbTriTargetCal, mTarget, nTarget);
% 
%             % Convert transformed dichromat-rendered image to display RGB
%             rgbDiTargetCal = rgbLin2RGB(rgbLinRenderedTarget, DispTarget);
%             rgbDiTargetImg = CalFormatToImage(rgbDiTargetCal, mTarget, nTarget);
% 
%             % Store for an optional target-resolution montage
%             triOutImgsTarget(:, :, :, kk) = rgbTriTargetImg;
%             diOutImgsTarget(:, :, :, kk)  = rgbDiTargetImg;
% 
%             % Recompute normalized achieved values at the target resolution
%             [~, infoAtTargetResolution] = infoFcn( ...
%                 LMSContrast_target, ...
%                 newLMSContrastTarget, ...
%                 imgParamsTarget, ...
%                 dichromatType, ...
%                 infoNormalizerTarget, ...
%                 DispTarget, ...
%                 infoParams);
% 
%             [~, distAtTargetResolution] = distortionFcn( ...
%                 LMS_target, ...
%                 newLMSTarget, ...
%                 imgParamsTarget, ...
%                 distortionNormalizerTarget, ...
%                 DispTarget, ...
%                 distortionParams);
% 
%             % Log target-resolution metadata/results
%             log(ii).metric(kk).targetResolution.applied   = true;
%             log(ii).metric(kk).targetResolution.size      = [mTarget nTarget];
%             log(ii).metric(kk).targetResolution.infoAch   = infoAtTargetResolution;
%             log(ii).metric(kk).targetResolution.distAch   = distAtTargetResolution;
%             log(ii).metric(kk).targetResolution.doClamp   = doClampTargetRGB;
% 
%             % Save target-resolution images into a folder near the sweep folder
%             if saveTargetResolutionImages
% 
%                 metricFolder = buildMetricFolderName(infoFcn, infoParams, distortionFcn);
% 
%                 % Example:
%                 %   distortion_20steps_128x128
%                 %   info_20steps_256x256
%                 saveRunFolderTarget = sprintf('%s_%dsteps_%dx%d', ...
%                     lower(char(sweepAxis)), nSteps, mTarget, nTarget);
% 
%                 saveSubdirTarget = fullfile(saveBase, pathName, metricFolder, saveRunFolderTarget);
% 
%                 if exist(saveSubdirTarget, 'dir') ~= 7
%                     mkdir(saveSubdirTarget);
%                 end
% 
%                 % Use fixed target + picked step in the filename to avoid ambiguity
%                 if strcmpi(fixedMode, 'dist')
%                     targetTag = sprintf('fixedDist_%0.3f', targetDistToCompare);
%                 else
%                     targetTag = sprintf('fixedInfo_%0.3f', targetInfoToCompare);
%                 end
%                 targetTag = strrep(targetTag, '.', 'p');
% 
%                 fileStem = sprintf('%s_%s_step_%03d', metricName, targetTag, idxClosest);
% 
%                 triToWrite = max(min(rgbTriTargetImg, 1), 0);
%                 diToWrite  = max(min(rgbDiTargetImg, 1), 0);
% 
%                 triPath = fullfile(saveSubdirTarget, sprintf('%s_trichromat.png', fileStem));
%                 diPath  = fullfile(saveSubdirTarget, sprintf('%s_%s.png', fileStem, lower(char(dichromatType))));
% 
%                 imwrite(triToWrite, triPath);
%                 imwrite(diToWrite, diPath);
% 
%                 log(ii).metric(kk).targetResolution.saveFolder = saveSubdirTarget;
%                 log(ii).metric(kk).targetResolution.triPath    = triPath;
%                 log(ii).metric(kk).targetResolution.diPath     = diPath;
%             else
%                 log(ii).metric(kk).targetResolution.saveFolder = '';
%                 log(ii).metric(kk).targetResolution.triPath    = '';
%                 log(ii).metric(kk).targetResolution.diPath     = '';
%             end
% 
%         else
%             log(ii).metric(kk).targetResolution.applied = false;
%         end
% 
%     end % metric loop
% 
% 
%     %% ========================= SOURCE-RESOLUTION VISUALIZATIONS =========================
% 
%     nMetrics = numel(infoNameList);
%     nCols    = nMetrics + 1;
%     nRows    = 2;   % trichromat row + dichromat-render row
% 
%     if strcmpi(fixedMode, 'dist')
%         figName = sprintf('%s — %s — fixed dist ~%.2f', ...
%             imgType, dichromatType, targetDistToCompare);
%     else
%         figName = sprintf('%s — %s — fixed info ~%.2f', ...
%             imgType, dichromatType, targetInfoToCompare);
%     end
% 
%     fig = figure('Color', 'w', 'Name', figName);
%     tl  = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
%     % Top-left: original trichromat
%     ax = nexttile(tl, 1);
%     imshow(rgbTriOrig, 'Parent', ax);
%     title(ax, 'Original (tri)', 'Interpreter', 'none');
% 
%     % Bottom-left: original dichromat render
%     ax = nexttile(tl, 1 + nCols);
%     imshow(rgbDiOrig, 'Parent', ax);
%     title(ax, sprintf('Original render (%s)', dichromatType), 'Interpreter', 'none');
% 
%     % One column per metric
%     for kk = 1:nMetrics
%         thisCol = kk + 1;
% 
%         % Top row: trichromat optimized image
%         ax = nexttile(tl, thisCol);
%         imshow(triOutImgs(:, :, :, kk), 'Parent', ax);
% 
%         if strcmpi(fixedMode, 'dist')
%             title(ax, sprintf('%s | tri\n(dist=%.3f, info=%.3f)', ...
%                 infoNameList{kk}, ...
%                 log(ii).metric(kk).distAtTarget, ...
%                 log(ii).metric(kk).infoAtTarget), ...
%                 'Interpreter', 'none');
%         else
%             title(ax, sprintf('%s | tri\n(info=%.3f, dist=%.3f)', ...
%                 infoNameList{kk}, ...
%                 log(ii).metric(kk).infoAtTarget, ...
%                 log(ii).metric(kk).distAtTarget), ...
%                 'Interpreter', 'none');
%         end
% 
%         % Bottom row: dichromat-rendered optimized image
%         ax = nexttile(tl, thisCol + nCols);
%         imshow(diOutImgs(:, :, :, kk), 'Parent', ax);
%         title(ax, sprintf('%s | di render', infoNameList{kk}), 'Interpreter', 'none');
%     end
% 
% 
%     %% ========================= CURVE PLOT =========================
% 
%     figure('Color', 'w', 'Name', sprintf('%s — curves by metric', imgType));
%     hold on;
% 
%     for kk = 1:nMetrics
%         plot(allDistCurves{kk}, allInfoCurves{kk}, '-o', ...
%             'LineWidth', 1.5, ...
%             'DisplayName', infoNameList{kk});
%     end
% 
%     if strcmpi(fixedMode, 'dist')
%         xline(targetDistToCompare, '--', 'Target dist');
%     else
%         yline(targetInfoToCompare, '--', 'Target info');
%     end
% 
%     grid on;
%     axis square;
%     xlabel('Achieved Distortion (normalized)');
%     ylabel('Achieved Info (normalized)');
% 
%     if strcmpi(fixedMode, 'dist')
%         title(sprintf('%s — %s — distortion sweep curves', imgType, dichromatType));
%     else
%         title(sprintf('%s — %s — info sweep curves', imgType, dichromatType));
%     end
% 
%     legend('Location', 'bestoutside');
%     hold off;
% 
% 
%     %% ========================= TARGET-RESOLUTION MONTAGE =========================
% 
%     if doApplyToTargetResolution && showTargetResolutionMontage
% 
%         if strcmpi(fixedMode, 'dist')
%             figNameTarget = sprintf('%s — %s — fixed dist ~%.2f — %dx%d', ...
%                 imgType, dichromatType, targetDistToCompare, mTarget, nTarget);
%         else
%             figNameTarget = sprintf('%s — %s — fixed info ~%.2f — %dx%d', ...
%                 imgType, dichromatType, targetInfoToCompare, mTarget, nTarget);
%         end
% 
%         figTarget = figure('Color', 'w', 'Name', figNameTarget);
%         tlTarget  = tiledlayout(figTarget, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
%         % Top-left: target-resolution original trichromat
%         ax = nexttile(tlTarget, 1);
%         imshow(rgbTriOrigTarget, 'Parent', ax);
%         title(ax, sprintf('Original (tri) %dx%d', mTarget, nTarget), 'Interpreter', 'none');
% 
%         % Bottom-left: target-resolution original dichromat render
%         ax = nexttile(tlTarget, 1 + nCols);
%         imshow(rgbDiOrigTarget, 'Parent', ax);
%         title(ax, sprintf('Original render (%s)', dichromatType), 'Interpreter', 'none');
% 
%         % One column per metric
%         for kk = 1:nMetrics
%             thisCol = kk + 1;
% 
%             % Top row: transformed trichromat image at target resolution
%             ax = nexttile(tlTarget, thisCol);
%             imshow(triOutImgsTarget(:, :, :, kk), 'Parent', ax);
% 
%             if strcmpi(fixedMode, 'dist')
%                 title(ax, sprintf('%s | tri\n(dist=%.3f, info=%.3f)', ...
%                     infoNameList{kk}, ...
%                     log(ii).metric(kk).targetResolution.distAch, ...
%                     log(ii).metric(kk).targetResolution.infoAch), ...
%                     'Interpreter', 'none');
%             else
%                 title(ax, sprintf('%s | tri\n(info=%.3f, dist=%.3f)', ...
%                     infoNameList{kk}, ...
%                     log(ii).metric(kk).targetResolution.infoAch, ...
%                     log(ii).metric(kk).targetResolution.distAch), ...
%                     'Interpreter', 'none');
%             end
% 
%             % Bottom row: dichromat-rendered transformed image at target resolution
%             ax = nexttile(tlTarget, thisCol + nCols);
%             imshow(diOutImgsTarget(:, :, :, kk), 'Parent', ax);
%             title(ax, sprintf('%s | di render', infoNameList{kk}), 'Interpreter', 'none');
%         end
%     end
% 
% end % image loop
% 
% 
% %% ========================= SAVE RESULTS =========================
% 
% % Keep the original save behavior: save the full log struct to the
% % ColorCorrection output directory.
% projectName = 'ColorCorrection';
% outputDir   = getpref(projectName, 'outputDir');
% 
% if strcmpi(fixedMode, 'dist')
%     saveFile = fullfile(outputDir, sprintf('fixedDist%.2f_%s_%dsteps.mat', ...
%         targetDistToCompare, dichromatType, nSteps));
% else
%     saveFile = fullfile(outputDir, sprintf('fixedInfo%.2f_%s_%dsteps.mat', ...
%         targetInfoToCompare, dichromatType, nSteps));
% end
% 
% save(saveFile, 'log');
% fprintf('\nSaved results to:\n  %s\n', saveFile);
% 
% % Original trailing line kept as-is
% k = 1;
% 
% % function t_compareFixedMetrics
% % 
% % % t_compareFixedMetrics.m
% % % Compare image enhancement outputs across multiple metrics while
% % % holding either distortion OR info approximately constant.
% % %
% % % Minimal change: add a switch:
% % %   fixedMode = 'dist' or 'info'
% % % and set:
% % %   targetDistToCompare OR targetInfoToCompare
% % %
% % % If fixedMode='dist': run distortion sweep, pick idx closest in distortion.
% % % If fixedMode='info': run info sweep, pick idx closest in info.
% % 
% % clear; clc;
% % 
% % %%%%%%%%%%%%%%%%%%%%%
% % % Images parameters %
% % %%%%%%%%%%%%%%%%%%%%%
% % imageTypes = {'flower1.png', 'fruit.png', 'Gaugin.png', 'ishi45.png'};
% % setType       = 1;
% % dichromatType = 'Deuteranopia';
% % 
% % % Size of image
% % m = 61; n = 61;
% % 
% % % How many steps to sample at
% % nSteps = 20;
% % 
% % % Clear the previous image data if it existed already in dropbox?
% % clearFlag = 0;
% % 
% % % Distortion metric
% % distortionFcn    = @computeDistortion_DE2000;
% % distortionParams = struct();
% % 
% % % Render function (for dichromat rendering)
% % renderFcn    = @DichromRenderLinear;
% % renderParams = struct();
% % 
% % % Info metrics
% % infoFcnList = {@computeInfo_regressSquared, ...
% %     @computeInfo_regressCIE};
% % 
% % infoNameList = { ...
% %     'regressSquared', ...
% %     'regressCIE'};
% % 
% % % Info metrics
% % % infoFcnList = {@computeInfo_LMdifference, ...
% % %     @computeInfo_regress, ...
% % %     @computeInfo_Wade };
% % % 
% % % infoNameList = { ...
% % %     'LMdifference', ...
% % %     'regress', ...
% % %     'Wade' };
% % 
% % % Some info metrics need params... regress does, others can use empty struct
% % infoParamsList = cell(size(infoFcnList));
% % for k = 1:numel(infoFcnList)
% %     if isequal(infoFcnList{k}, @computeInfo_regressCIE) || isequal(infoFcnList{k}, @computeInfo_regressSquared)
% %         infoParamsList{k} = struct( ...
% %             'predictingWhat',     'L,M,S', ...
% %             'predictingFromWhat', 'L and S');
% %     else
% %         infoParamsList{k} = struct();
% %     end
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % NEW: choose whether to fix distortion or fix info
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fixedMode = 'dist';   % 'dist' or 'info'
% % 
% % % If fixedMode='dist', we match solutions at this achieved distortion:
% % targetDistToCompare = 0.4;
% % 
% % % If fixedMode='info', we match solutions at this achieved info:
% % targetInfoToCompare = 0.4; %???? idk what to set this to
% % 
% % log = struct();
% % 
% % %% Main loop over images
% % for ii = 1:numel(imageTypes)
% % 
% %     imgType = imageTypes{ii};
% %     fprintf('\n====================================================\n');
% %     fprintf('IMAGE %d/%d: %s | %s\n', ii, numel(imageTypes), imgType, dichromatType);
% %     fprintf('====================================================\n');
% % 
% %     % 1) Load images and parameters
% %     [LMSCalFormat, ~, ~, ~, Disp, imgParams, pathName] = colorCorrectionGenerateImages( ...
% %         imgType, setType, m, n, dichromatType, clearFlag);
% % 
% %     % Log some relevant stuff for later
% %     log(ii).imgType = imgType;
% %     log(ii).dichromatType = dichromatType;
% %     log(ii).pathName = pathName;
% %     log(ii).fixedMode = fixedMode;
% %     log(ii).targetDistToCompare = targetDistToCompare;
% %     log(ii).targetInfoToCompare = targetInfoToCompare;
% % 
% %     % 2) Build reference images
% %     rgbLinCalFormat = Disp.M_cones2rgb * LMSCalFormat;
% % 
% %     tol = 1e-3;
% %     mn = min(rgbLinCalFormat(:));
% %     mx = max(rgbLinCalFormat(:));
% %     if mn < -tol || mx > 1+tol
% %         warning('Clamping: out of range beyond tol (min=%.6g, max=%.6g).', mn, mx);
% %     end
% %     rgbLinCalFormat = min(max(rgbLinCalFormat, 0), 1);
% % 
% %     rgbTriOrig = CalFormatToImage(rgbLin2RGB(rgbLinCalFormat, Disp), m, n);
% % 
% %     % Render original for the dichromat and convert to image format
% %     [~, rgbLinTriRenderedOrig, ~] = DichromRenderLinear(LMSCalFormat, dichromatType, Disp);
% %     rgbDiOrig = CalFormatToImage(rgbLin2RGB(rgbLinTriRenderedOrig, Disp), m, n);
% % 
% %     % 3) Preallocate images to compare later
% %     nMetrics = numel(infoFcnList);
% %     triOutImgs = zeros(m, n, 3, nMetrics);
% %     diOutImgs  = zeros(m, n, 3, nMetrics);
% % 
% %     % Store curves per metric for the curve plot
% %     allInfoCurves = cell(nMetrics,1);
% %     allDistCurves = cell(nMetrics,1);
% % 
% %     % Loop over info metrics
% %     for kk = 1:nMetrics
% % 
% %         infoFcn    = infoFcnList{kk};
% %         infoParams = infoParamsList{kk};
% %         metricName = infoNameList{kk};
% % 
% %         fprintf('\n  ---- Metric %d/%d: %s ----\n', kk, nMetrics, metricName);
% % 
% %         % 4) Create daltonizer object for this metric
% %         theDaltonizer = daltonize(infoFcn, infoParams, distortionFcn, distortionParams, renderFcn, renderParams, Disp);
% % 
% %         % 5) Run the sweep (distortion or info) based on fixedMode
% %         if strcmpi(fixedMode, 'dist')
% %             sweepAxis = 'distortion';
% %         else
% %             sweepAxis = 'info';
% %         end
% % 
% %         [LMSSweep, rgbLinSweep, LMSRenderedSweep, rgbLinRenderedSweep, TSweep, ...
% %             targetInfoNorm, targetDistNorm, infoNormAch, distNormAch] = computeSweep( ...
% %             theDaltonizer, LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis);
% % 
% %         % 6) Pick the fixed-(dist or info) solution
% %         if strcmpi(fixedMode, 'dist')
% %             [minAbsErr, idxClosest] = min(abs(distNormAch - targetDistToCompare));
% %             pickedTarget = targetDistToCompare;
% %             pickedAch    = distNormAch(idxClosest);
% %             pickedOther  = infoNormAch(idxClosest);
% %             fprintf('     targetDist=%.3f | closest achievedDist=%.3f (|err|=%.3f) | achievedInfo=%.3f\n', ...
% %                 pickedTarget, pickedAch, minAbsErr, pickedOther);
% %         else
% %             [minAbsErr, idxClosest] = min(abs(infoNormAch - targetInfoToCompare));
% %             pickedTarget = targetInfoToCompare;
% %             pickedAch    = infoNormAch(idxClosest);
% %             pickedOther  = distNormAch(idxClosest);
% %             fprintf('     targetInfo=%.3f | closest achievedInfo=%.3f (|err|=%.3f) | achievedDist=%.3f\n', ...
% %                 pickedTarget, pickedAch, minAbsErr, pickedOther);
% %         end
% % 
% %         infoAtTarget = infoNormAch(idxClosest);
% %         distAtTarget = distNormAch(idxClosest);
% % 
% %         % 7) Extract the image for that solution
% %         % Optimized *trichromat* image (linear RGB in CalFormat)
% %         rgbTriCal = rgbLin2RGB(rgbLinSweep{idxClosest}, Disp);
% %         rgbTriImg = CalFormatToImage(rgbTriCal, m, n);
% % 
% %         % Dichromat-rendered version of that optimized image
% %         rgbDiCal  = rgbLin2RGB(rgbLinRenderedSweep{idxClosest}, Disp);
% %         rgbDiImg  = CalFormatToImage(rgbDiCal, m, n);
% % 
% %         % Store for montage
% %         triOutImgs(:,:,:,kk) = rgbTriImg;
% %         diOutImgs(:,:,:,kk)  = rgbDiImg;
% % 
% %         % 8) Store data for the curve plots
% %         log(ii).metric(kk).name           = metricName;
% %         log(ii).metric(kk).infoFcn        = func2str(infoFcn);
% %         log(ii).metric(kk).infoParams     = infoParams;
% % 
% %         log(ii).metric(kk).idxClosest     = idxClosest;
% %         log(ii).metric(kk).distAtTarget   = distAtTarget;
% %         log(ii).metric(kk).infoAtTarget   = infoAtTarget;
% %         log(ii).metric(kk).absErrAtTarget = minAbsErr;
% % 
% %         log(ii).metric(kk).T = TSweep{idxClosest};
% % 
% %         log(ii).metric(kk).infoNormAch    = infoNormAch(:);
% %         log(ii).metric(kk).distNormAch    = distNormAch(:);
% %         log(ii).metric(kk).targetInfoNorm = targetInfoNorm(:);
% %         log(ii).metric(kk).targetDistNorm = targetDistNorm(:);
% % 
% %         allInfoCurves{kk} = infoNormAch(:);
% %         allDistCurves{kk} = distNormAch(:);
% % 
% %     end % metrics loop
% % 
% %     %% %%%%%%%%%%%%%%%
% %     % Visualizations %
% %     %%%%%%%%%%%%%%%%%%
% % 
% %     nMetrics = numel(infoNameList);
% %     nCols = nMetrics + 1;
% %     nRows = 2; % trichromat and dichromat
% % 
% %     if strcmpi(fixedMode,'dist')
% %         figName = sprintf('%s — %s — fixed dist ~%.2f', imgType, dichromatType, targetDistToCompare);
% %     else
% %         figName = sprintf('%s — %s — fixed info ~%.2f', imgType, dichromatType, targetInfoToCompare);
% %     end
% % 
% %     fig = figure('Color','w','Name',figName);
% %     tl  = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
% % 
% %     % Top row: original trichromat
% %     ax = nexttile(tl, 1);
% %     imshow(rgbTriOrig, 'Parent', ax);
% %     title(ax, 'Original (tri)', 'Interpreter','none');
% % 
% %     % Bottom row: original dichromat render
% %     ax = nexttile(tl, 1 + nCols);
% %     imshow(rgbDiOrig, 'Parent', ax);
% %     title(ax, sprintf('Original render (%s)', dichromatType), 'Interpreter','none');
% % 
% %     % Metric columns
% %     for kk = 1:nMetrics
% %         thisCol = kk + 1;
% % 
% %         ax = nexttile(tl, thisCol);
% %         imshow(triOutImgs(:,:,:,kk), 'Parent', ax);
% % 
% %         if strcmpi(fixedMode,'dist')
% %             title(ax, sprintf('%s | tri\n(dist=%.3f, info=%.3f)', ...
% %                 infoNameList{kk}, log(ii).metric(kk).distAtTarget, log(ii).metric(kk).infoAtTarget), 'Interpreter','none');
% %         else
% %             title(ax, sprintf('%s | tri\n(info=%.3f, dist=%.3f)', ...
% %                 infoNameList{kk}, log(ii).metric(kk).infoAtTarget, log(ii).metric(kk).distAtTarget), 'Interpreter','none');
% %         end
% % 
% %         ax = nexttile(tl, thisCol + nCols);
% %         imshow(diOutImgs(:,:,:,kk), 'Parent', ax);
% %         title(ax, sprintf('%s | di render', infoNameList{kk}), 'Interpreter','none');
% %     end
% % 
% %     % Second graph: achieved info vs achieved distortion curves
% %     figure('Color','w','Name',sprintf('%s — curves by metric', imgType));
% %     hold on;
% %     for kk = 1:nMetrics
% %         plot(allDistCurves{kk}, allInfoCurves{kk}, '-o','LineWidth', 1.5, 'DisplayName', infoNameList{kk});
% %     end
% % 
% %     if strcmpi(fixedMode,'dist')
% %         xline(targetDistToCompare, '--', 'Target dist');
% %     else
% %         yline(targetInfoToCompare, '--', 'Target info');
% %     end
% % 
% %     grid on; axis square;
% %     xlabel('Achieved Distortion (normalized)');
% %     ylabel('Achieved Info (normalized)');
% % 
% %     if strcmpi(fixedMode,'dist')
% %         title(sprintf('%s — %s — distortion sweep curves', imgType, dichromatType));
% %     else
% %         title(sprintf('%s — %s — info sweep curves', imgType, dichromatType));
% %     end
% % 
% %     legend('Location','bestoutside');
% %     hold off;
% % 
% % end % images loop
% % 
% % %% %%%%%%%%%%%%%
% % % SAVE RESULTS %
% % %%%%%%%%%%%%%%%%
% % projectName = 'ColorCorrection';
% % outputDir   = getpref(projectName, 'outputDir');
% % 
% % if strcmpi(fixedMode,'dist')
% %     saveFile = fullfile(outputDir, sprintf('fixedDist%.2f_%s_%dsteps.mat', ...
% %         targetDistToCompare, dichromatType, nSteps));
% % else
% %     saveFile = fullfile(outputDir, sprintf('fixedInfo%.2f_%s_%dsteps.mat', ...
% %         targetInfoToCompare, dichromatType, nSteps));
% % end
% % 
% % save(saveFile, 'log');
% % fprintf('\nSaved results to:\n  %s\n', saveFile);
% % 
% % k = 1;
