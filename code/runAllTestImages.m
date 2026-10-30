%% COLOR CORRECTION: INFORMATION–DISTORTION OPTIMIZATION
%
% PURPOSE
%   This script evaluates image transformations designed to improve the
%   color information available to dichromatic observers while minimizing
%   changes to the original image.
%
%
% OPTIMIZATION
%   The script supports two types of constrained optimization:
%
%   1. Information sweep ('info'):
%      Minimize distortion subject to achieving a target information score.
%
%   2. Distortion sweep ('distortion'):
%      Maximize information subject to a target distortion level.
%
%   Target values are sampled between the solutions obtained by
%   minimizing distortion (lambda = 0) and maximizing information
%   (lambda = 1).

% Image types
% imageTypes = {'Gaugin.png'};
% imageTypes = {'map.png','tree.png','car.png','cell.png','pool.png','flower1.png', 'fruit.png', 'Gaugin.png', 'ishi45.png'};   % all 4 test images (updates 1-3 at 64x64, 10/29)

% Gray image with one square in the middle.
% The square changes only the missing cone (M for Deuteranopia).
% setType (below) = number of squares, so setType = 1 gives one square.
% imageTypes = {'gray'};

% Gray image with a red-green (L-M) square plus noise in L and M
% (settings in loadLMSvalues): same noise in L and M vs different noise.
imageTypes = {'noiseSquareCorrelated', 'noiseSquareUncorrelated'};
% imageTypes = {'tree.png'};
% imageTypes = {'flower1.png', ...
%     'Gaugin.png'};
% imageTypes = {'Gaugin.png'};          % convex-hull test: compare with info_10steps_centralFiniteDifferences

% Set type (choose 1 unless ishihara)
setType       = 1;

% Dichromacy type
dichromatType = 'Deuteranopia';

% Image size (keep <60 for fast testing)

% m = 100;
% n = 100;
m = 128;
n = 128;
m = 20;
n = 20;


% How many steps do you want in a transformation sweep?
nSteps = 10;

sweepAxis = 'info';
%   'info'        = run info sweep only. this chooses set of info scores and minimizes distortion.
%   'distortion'  = run distortion sweep only
%   'both'        = run both sweeps

% Which update to look at (see updateLog.txt for what each number means)
%   []  = run a new sweep with the current code; results go in info_10steps
%   N   = LOOK UP ONLY: load the finished sweep in info_10steps_update<N> and
%         make the figures. If it doesn't exist, stop; nothing is computed.
updateNumber = [];

% TEST: replace the L cone plane with the M cone plane, so the image is
% [M M S], and look at how well M is predicted. true = run only this test
% (no sweep); false = normal run.
runMMSTest = false;


%% Define metric functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% INFO FUNCTION %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoFcn = @computeInfo_regressCIE;
% infoFcn = @computeInfo_regressSquared;
% infoFcn = @computeInfo_LMdifference;
% infoFcn = @computeInfo_Wade;

% infoParams only necessary when infoFcn is regress
infoParams = struct( ...
    'predictingWhat',     'L-M', ...          % options: 'M', 'L,M,S', 'L-M' (score R^2)
    'predictingFromWhat', 'L and S', ...        % options: 'L and S', 'opponent' (dichromat's perceived L*a*b*; use weightPenaltySD in DE, e.g. 1), 'deltaL and deltaS', 'deltaL+M and delta S'
    'whichWeightPenalty', 'ridge', ...      % options: 'none', 'noise', 'ridge' (penalizes huge regression weights)
    'weightPenaltySD',    .02);            % noise SD (or matching ridge strength), in cone contrast
% infoParams = struct( ...
%     'predictingWhat',     'L,M,S', ...          % options: 'M' or 'L-M' or 'L,M,S'
%     'predictingFromWhat', 'L and S');           % options: 'L and S' or 'deltaL and deltaS' or 'deltaL+M and delta S'
% infoParams = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% DISTORTION FUNCTION %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distortionFcn = @computeDistortion_squared;
distortionFcn = @computeDistortion_DE2000;
distortionParams = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% RENDER FUNCTION %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
renderFcn = @DichromRenderLinear;
renderParams = struct();


%% Loop over all image types
for ii = 1:numel(imageTypes)

    imgType = imageTypes{ii};
    fprintf('Processing image type: %s\n', imgType);

    % Generate input image
    [LMSCalFormat, rgbLinCalFormat, LMSCalFormatRendered, rgbLinCalFormatRendered, Disp, imgParams, pathName] = colorCorrectionGenerateImages(imgType, setType, m, n, dichromatType, 0);

    % ---- CHECK: is there contrast only in the M cone? ----
    LMSContrastCheck = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;
    maxContrastL = max(abs(LMSContrastCheck(1,:)));
    maxContrastM = max(abs(LMSContrastCheck(2,:)));
    maxContrastS = max(abs(LMSContrastCheck(3,:)));
    fprintf('\nLargest contrast in each cone:  L = %.3g   M = %.3g   S = %.3g\n', maxContrastL, maxContrastM, maxContrastS);
    if maxContrastL > 1e-10 || maxContrastS > 1e-10
        warning('There is contrast in L or S, not only in M.');
    end

    % ---- TEST: image = [M M S] (L cone plane replaced by M) ----
    if runMMSTest

        % 1. Original cone contrast
        originalContrast = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

        % 2. Copy the M plane into the L plane: [M M S]
        mmsContrast      = originalContrast;
        mmsContrast(1,:) = 2 * (originalContrast(1,:) - originalContrast(2,:));
        % mmsContrast(3,:) = originalContrast(2,:);        
        mmsLMS           = mmsContrast .* Disp.grayLMS + Disp.grayLMS;

        % 3. Info of the original and of [M M S] (raw info, normalizer = 1)
        infoOriginal = computeInfo_regressCIE(originalContrast, originalContrast, imgParams, dichromatType, 1, Disp, infoParams);
        [infoMMS, ~, mmsPredictedLMS, betaMMS] = computeInfo_regressCIE(originalContrast, mmsContrast, imgParams, dichromatType, 1, Disp, infoParams);

        % 4. Dichromat view of [M M S]
        [~, mmsDiRGBlin] = DichromRenderLinear(mmsLMS, dichromatType, Disp);

        % 5. Linear RGB of each picture, clipped to the display range [0, 1]
        originalRGBlin     = min(max(rgbLinCalFormat, 0), 1);
        originalDiRGBlin   = min(max(rgbLinCalFormatRendered, 0), 1);
        mmsRGBlin          = min(max(Disp.M_cones2rgb * mmsLMS, 0), 1);
        mmsDiRGBlin        = min(max(mmsDiRGBlin, 0), 1);
        mmsPredictedRGBlin = min(max(Disp.M_cones2rgb * mmsPredictedLMS, 0), 1);

        % 6. Convert to displayable images
        originalImage     = CalFormatToImage(rgbLin2RGB(originalRGBlin, Disp), m, n);
        originalDiImage   = CalFormatToImage(rgbLin2RGB(originalDiRGBlin, Disp), m, n);
        mmsImage          = CalFormatToImage(rgbLin2RGB(mmsRGBlin, Disp), m, n);
        mmsDiImage        = CalFormatToImage(rgbLin2RGB(mmsDiRGBlin, Disp), m, n);
        mmsPredictedImage = CalFormatToImage(rgbLin2RGB(mmsPredictedRGBlin, Disp), m, n);

        % 7. Show them side by side
        figure;
        tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        imshow(originalImage);
        title('Original');

        nexttile;
        imshow(originalDiImage);
        title('Original (dichromat)');

        nexttile;
        imshow(mmsImage);
        title('[M M S]');

        nexttile;
        imshow(mmsDiImage);
        title('[M M S] (dichromat)');

        nexttile;
        imshow(mmsPredictedImage);
        title('Predicted from [M M S]');

        % Skip the sweep for this image
        % continue
    end

    % Set up daltonizer object
    theDaltonizer = daltonize(infoFcn, infoParams,distortionFcn, distortionParams,renderFcn, renderParams, Disp);

    %% Info and distortion sweeps
    %   (1) Run an info sweep (minimize distortion at each target info)
    %   (2) Run a distortion sweep (maximize info at each target distortion)

    nCols  = ceil(sqrt(nSteps));                % grid columns (near-square)
    nRows  = ceil(nSteps / nCols);              % grid rows

    ranInfo = false;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1) INFO SWEEP — minimize distortion subject to target info
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(sweepAxis,'info') || strcmp(sweepAxis,'both')
        ranInfo = true;
        sweepAxis_info = 'info';

        % Compute the sweep over info targets!!
        [LMSSweep_info, rgbLinSweep_info, LMSRenderedSweep_info, rgbLinRenderedSweep_info,TSweep_info, targetInfoNorm_info, targetDistNorm_info, infoNormAch_info, distNormAch_info] = ...
            computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType,nSteps, pathName, sweepAxis_info, 'updateNumber', updateNumber);

        % Target vs achieved info
        fprintf('\nInfo sweep (n=%d): target vs achieved info\n', numel(targetInfoNorm_info));
        disp([targetInfoNorm_info(:), infoNormAch_info(:)]);

        % Turn each step's result into pictures (used by the figures below)
        RGBTriImgFormat_optAllInfo = zeros(m, n, 3, nSteps);
        RGBDiImgFormat_optAllInfo  = zeros(m, n, 3, nSteps);

        for i = 1:nSteps
            % Keep RGB values in the display range [0, 1]
            rgbLinSweep_info{i}(rgbLinSweep_info{i}<0) = 0;
            rgbLinSweep_info{i}(rgbLinSweep_info{i}>1) = 1;

            % Trichromat (transformed) image for step i
            RGBTriCalFormat_opt = rgbLin2RGB(rgbLinSweep_info{i}, Disp);
            RGBTriImgFormat_optAllInfo(:,:,:,i) = CalFormatToImage(RGBTriCalFormat_opt,m,n);

            % Dichromat view of the transformed image for step i
            RGBDiCalFormat_opt = rgbLin2RGB(rgbLinRenderedSweep_info{i}, Disp);
            RGBDiImgFormat_optAllInfo(:,:,:,i) = CalFormatToImage(RGBDiCalFormat_opt,m,n);
        end


        % ---- CHECKS FOR EACH STEP: predicted image and regression weights ----
        % The info score rebuilds the full color image from the dichromat's
        % channels using a regression. For each step we:
        %   1. save and show that rebuilt ("predicted") trichromat image
        %   2. save the regression weights (beta), and print the largest
        %      weight and how well M was rebuilt (R^2, 1 = perfect).
        % A huge weight with a near-perfect rebuild = a tiny signal is being
        % scaled up by the regression.

        % Folder to save the checks in (inside this run's results folder)
        runFolder = sprintf('info_%dsteps', nSteps);
        if ~isempty(updateNumber)
            runFolder = sprintf('%s_update%d', runFolder, updateNumber);
        end
        metricFolder = buildMetricFolderName(infoFcn, infoParams, distortionFcn);
        checkDir = fullfile(getpref('ColorCorrection','outputDir'), 'testImagesTransformed', ...
            pathName, metricFolder, runFolder, 'regressionChecks');
        if ~exist(checkDir, 'dir')
            mkdir(checkDir);
        end

        % Original M cone contrast (what the regression tries to rebuild)
        originalContrast = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;
        originalM        = originalContrast(2,:);

        % Places to keep the results of every step
        predictedImages = zeros(m, n, 3, nSteps);
        allBeta         = cell(1, nSteps);
        largestWeight   = zeros(1, nSteps);
        rebuildR2       = zeros(1, nSteps);

        for i = 1:nSteps

            % 1. This step's transform
            T = TSweep_info{i};

            % 2. Run the info score on it. Output 6 = rebuilt LMS image,
            %    output 7 = regression weights. (Normalizers = 1: these
            %    two outputs don't depend on them.)
            [~, ~, ~, ~, ~, predictedLMS, beta] = lossFunction('lambda', 0, T(:), ...
                LMSCalFormat, imgParams, dichromatType, infoFcn, distortionFcn, ...
                1, 1, Disp, infoParams);

            % 3. How well was M rebuilt? (R^2: 1 = perfectly)
            predictedContrast = (predictedLMS - Disp.grayLMS) ./ Disp.grayLMS;
            predictedM        = predictedContrast(2,:);
            errorSum          = sum((originalM - predictedM).^2);
            totalSum          = sum((originalM - mean(originalM)).^2);
            rebuildR2(i)      = 1 - errorSum / max(totalSum, eps);

            % 4. Keep the weights and the largest one
            allBeta{i}       = beta;
            largestWeight(i) = max(abs(beta(:)));

            % 5. Turn the rebuilt LMS image into a picture
            predictedRGBlinear = Disp.M_cones2rgb * predictedLMS;
            predictedRGBlinear = min(max(predictedRGBlinear, 0), 1);   % keep in display range
            predictedRGB       = rgbLin2RGB(predictedRGBlinear, Disp);
            predictedImages(:,:,:,i) = CalFormatToImage(predictedRGB, m, n);

            % 6. Save this step's picture
            imwrite(predictedImages(:,:,:,i), fullfile(checkDir, sprintf('predicted_step_%02d.png', i)));
        end

        % ---- Print the largest regression weight at each step (table) ----
        Step          = (1:nSteps)';
        LargestWeight = largestWeight(:);
        weightTable   = table(Step, LargestWeight);

        fprintf('\nLargest regression weight at each step:\n');
        disp(weightTable)

    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) DISTORTION SWEEP — maximize info subject to target distortion
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ranDist = false;
    if strcmp(sweepAxis,'distortion') || strcmp(sweepAxis,'both')
        ranDist = true;
        sweepAxis_dist = 'distortion';

        % Compute the sweep over info targets!!
        [LMSSweep_dist, rgbLinSweep_dist, LMSRenderedSweep_dist, rgbLinRenderedSweep_dist, TSweep_dist, targetInfoNorm_dist, targetDistNorm_dist, infoNormAch_dist, distNormAch_dist] = ...
            computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_dist, 'updateNumber', updateNumber);

        % Target vs achieved distortion (normalized)
        fprintf('\nDistortion sweep (n=%d): target vs achieved distortion\n', numel(targetDistNorm_dist));
        disp([targetDistNorm_dist(:), distNormAch_dist(:)]);

        % Turn each step's result into pictures (used by the figures below)
        RGBTriImgFormat_optAllDist = zeros(m, n, 3, nSteps);
        RGBDiImgFormat_optAllDist  = zeros(m, n, 3, nSteps);

        for i = 1:nSteps
            % Keep RGB values in the display range [0, 1]
            rgbLinSweep_dist{i}(rgbLinSweep_dist{i}<0) = 0;
            rgbLinSweep_dist{i}(rgbLinSweep_dist{i}>1) = 1;

            % Trichromat (transformed) image for step i
            RGBTriCalFormat_opt = rgbLin2RGB(rgbLinSweep_dist{i}, Disp);
            RGBTriImgFormat_optAllDist(:,:,:,i) = CalFormatToImage(RGBTriCalFormat_opt,m,n);

            % Dichromat view of the transformed image for step i
            RGBDiCalFormat_opt = rgbLin2RGB(rgbLinRenderedSweep_dist{i}, Disp);
            RGBDiImgFormat_optAllDist(:,:,:,i) = CalFormatToImage(RGBDiCalFormat_opt,m,n);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Info-Distortion curve %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure();

    % If you did the info sweep, then go ahead and plot the info-distortion
    % curve for that sweep
    nPlot = 20;
    if ranInfo
        idx = 1:min(nPlot, numel(distNormAch_info));
        plot(distNormAch_info(idx), infoNormAch_info(idx), 'ko-', 'LineWidth', 2,'MarkerFaceColor','w', 'MarkerSize',10 , 'DisplayName', 'Info sweep (minimize distortion)');
        hold on;
        for k = idx
            yline(targetInfoNorm_info(k), ':', 'HandleVisibility','off');
        end
    end

    % If you did the distortion sweep, then go ahead and plot the info-distortion
    % curve for that sweep
    if ranDist
        idx = 1:min(nPlot, numel(distNormAch_dist));
        plot(distNormAch_dist(idx), infoNormAch_dist(idx), 'ks-', 'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerSize',10 ,'DisplayName', 'Distortion sweep (maximize info)');
        hold on;
        for k = idx
            xline(targetDistNorm_dist(k), ':', 'HandleVisibility','off');
        end
    end

    % Make it pretty
    grid off; axis square;
    ax = gca;
    ax.LineWidth = 1.5;
    ax.FontSize = 15;
    xlabel('Achieved Distortion (normalized)','FontSize',15);
    ylabel('Achieved Info (normalized)','FontSize',15);
    title(sprintf('%s — %s — %d steps', imgType, dichromatType, nSteps),'FontSize',20);

    if ranInfo && ranDist
        subtitle('Achieved Info vs Distortion for both sweeps');
    elseif ranInfo
        subtitle('Achieved Info vs Distortion (info sweep)');
    else
        subtitle('Achieved Info vs Distortion (distortion sweep)');
    end

    legend('Location','southeast');


    %% Plot each sweep as 2 rows x nSteps columns
    % Row 1 = Trichromat optimized images
    % Row 2 = Dichromat renderings
    % Above each column = Achieved information score

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFO SWEEP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ranInfo

        % 3 rows: transformed (trichromat), predicted, dichromat
        figure('Color','w', ...
            'Position',[100 100 150*nSteps 520], ...
            'Name',sprintf('%s - %s - Info Sweep',imgType,dichromatType));

        % 'compact' padding leaves room for the row labels on the left
        % 'compact' tile spacing leaves a small gap between rows
        tl = tiledlayout(3,nSteps, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        for i = 1:nSteps

            % Top row: Trichromat (transformed image)
            nexttile(i);
            imshow(RGBTriImgFormat_optAllInfo(:,:,:,i));
            axis image off;

            title(sprintf('%.3f',infoNormAch_info(i)), ...
                'FontSize',15);

            % Row label on the left (axis off hides labels, so turn it on)
            if i == 1
                ylabel('Trichromat', 'FontSize', 15, 'Visible', 'on');
            end

            % Middle row: Predicted image (rebuilt by the info regression,
            % computed in the regression checks above)
            nexttile(nSteps+i);
            imshow(predictedImages(:,:,:,i));
            axis image off;

            if i == 1
                ylabel('Predicted', 'FontSize', 15, 'Visible', 'on');
            end

            % Bottom row: Dichromat
            nexttile(2*nSteps+i);
            imshow(RGBDiImgFormat_optAllInfo(:,:,:,i));
            axis image off;

            if i == 1
                ylabel('Dichromat', 'FontSize', 15, 'Visible', 'on');
            end

        end

        % sgtitle(sprintf('%s - %s - Info Sweep',imgType,dichromatType), ...
        %     'FontSize',14);

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFO SWEEP: info-distortion curve with the images
    %   The images sit in a row along the bottom of the plot, one pair per
    %   chosen step: trichromat (top), dichromat (bottom). A grey line
    %   connects each pair to its point on the curve.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ranInfo

        % Settings
        stepsToShow = unique(round(linspace(1, nSteps, 5)));   % which steps get images (5, incl. first and last)
        thumbSize   = 0.10;             % largest image size, as a fraction of the figure

        % 1. Draw the curve
        figure('Color','w', 'Position',[100 100 1000 900], ...
            'Name',sprintf('%s - %s - Info vs Distortion with images',imgType,dichromatType));
        curveAx = axes('Position',[0.10 0.08 0.82 0.82]);   % [left bottom width height], fraction of the figure
        plot(curveAx, distNormAch_info, infoNormAch_info, 'ko-', ...
            'LineWidth',2, 'MarkerFaceColor','w', 'MarkerSize',10);
        curveAx.LineWidth = 1.5;
        curveAx.FontSize  = 15;
        xlabel(curveAx, 'Achieved Distortion (normalized)');
        ylabel(curveAx, 'Achieved Info (normalized)');
        title(curveAx, sprintf('%s — %s — %d steps', imgType, dichromatType, nSteps), 'FontSize',20);
        subtitle(curveAx, 'Images: trichromat (top), dichromat (bottom)');

        % 2. Image size: shrink if needed so the pairs fit side by side
        axPos     = curveAx.Position;
        nShow     = numel(stepsToShow);
        thumbSize = min(thumbSize, 0.9 * axPos(3) / nShow);

        % 3. Axis limits: leave an empty band at the bottom for the images.
        %    The lowest point on the curve sits just above the band.
        bandHeight = 2*thumbSize + 0.05;                % 2 images + room for labels, fraction of the figure
        bandFrac   = bandHeight / axPos(4);             % same, as a fraction of the axes
        xRange = max(distNormAch_info) - min(distNormAch_info);
        yRange = max(infoNormAch_info) - min(infoNormAch_info);
        yTop    = max(infoNormAch_info) + 0.05*yRange;
        yBottom = (min(infoNormAch_info) - bandFrac*yTop) / (1 - bandFrac);
        xlim(curveAx, [min(distNormAch_info) - 0.05*xRange, max(distNormAch_info) + 0.05*xRange]);
        ylim(curveAx, [yBottom, yTop]);
        xLimits = xlim(curveAx);
        yLimits = ylim(curveAx);

        % 4. One image pair per chosen step, evenly spaced along the bottom
        for k = 1:nShow
            i = stepsToShow(k);

            % Where this pair goes (fraction of the figure)
            pairCenterX = axPos(1) + (k - 0.5) * axPos(3) / nShow;
            imageLeft   = pairCenterX - thumbSize/2;
            dichromatBottom  = axPos(2) + 0.01;
            trichromatBottom = dichromatBottom + thumbSize;

            % Where this step's point is on the curve (fraction of the figure)
            pointX = axPos(1) + axPos(3) * (distNormAch_info(i) - xLimits(1)) / diff(xLimits);
            pointY = axPos(2) + axPos(4) * (infoNormAch_info(i) - yLimits(1)) / diff(yLimits);

            % Grey line from the point to the top of its pair
            annotation('line', [pointX, pairCenterX], [pointY, trichromatBottom + thumbSize + 0.025], ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 1);

            % Trichromat image (top), labeled with its info value
            axes('Position', [imageLeft, trichromatBottom, thumbSize, thumbSize]);
            imshow(RGBTriImgFormat_optAllInfo(:,:,:,i));
            axis image off;
            title(sprintf('%.3f', infoNormAch_info(i)), 'FontSize', 12);

            % Dichromat image (bottom)
            axes('Position', [imageLeft, dichromatBottom, thumbSize, thumbSize]);
            imshow(RGBDiImgFormat_optAllInfo(:,:,:,i));
            axis image off;
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISTORTION SWEEP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ranDist

        figure('Color','w', ...
            'Position',[100 100 150*nSteps 350], ...
            'Name',sprintf('%s - %s - Distortion Sweep',imgType,dichromatType));

        % 'compact' padding leaves room for the row labels on the left
        % 'compact' tile spacing leaves a small gap between rows
        tl = tiledlayout(2,nSteps, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        for i = 1:nSteps

            % Top row: Trichromat
            nexttile(i);
            imshow(RGBTriImgFormat_optAllDist(:,:,:,i));
            axis image off;

            title(sprintf('%.3f',infoNormAch_dist(i)), ...
                'FontSize',15);

            % Row label on the left (axis off hides labels, so turn it on)
            if i == 1
                ylabel('Trichromat', 'FontSize', 15, 'Visible', 'on');
            end

            % Bottom row: Dichromat
            nexttile(nSteps+i);
            imshow(RGBDiImgFormat_optAllDist(:,:,:,i));
            axis image off;

            if i == 1
                ylabel('Dichromat', 'FontSize', 15, 'Visible', 'on');
            end

        end

        % sgtitle(sprintf('%s - %s - Distortion Sweep',imgType,dichromatType), ...
        %     'FontSize',14);

    end


end