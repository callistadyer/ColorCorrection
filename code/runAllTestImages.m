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
imageTypes = {'flower1.png'};
% imageTypes = {'flower1.png', ...
%     'fruit.png', ...
%     'Gaugin.png', ...
%     'ishi45.png'};

% Set type (choose 1 unless ishihara)
setType       = 1;

% Dichromacy type
dichromatType = 'Deuteranopia';

% Image size (keep <60 for fast testing)
m = 30;
n = 30;

% How many steps do you want in a transformation sweep?
nSteps = 20;

sweepAxis = 'info';
%   'info'        = run info sweep only. this chooses set of info scores and minimizes distortion. 
%   'distortion'  = run distortion sweep only
%   'both'        = run both sweeps 


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
    'predictingWhat',     'L,M,S', ...          % options: 'M' or 'L-M' or 'L,M,S'
    'predictingFromWhat', 'L and S');           % options: 'L and S' or 'deltaL and deltaS' or 'deltaL+M and delta S'
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
            computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType,nSteps, pathName, sweepAxis_info);
   
        % Target vs achieved info
        fprintf('\nInfo sweep (n=%d): target vs achieved info\n', numel(targetInfoNorm_info));
        disp([targetInfoNorm_info(:), infoNormAch_info(:)]);

        % Image montage
        RGBTriImgFormat_optAllInfo = zeros(m, n, 3, nSteps);
        RGBDiImgFormat_optAllInfo  = zeros(m, n, 3, nSteps);

        for i = 1:nSteps
            rgbLinSweep_info{i}(rgbLinSweep_info{i}<0) = 0;
            rgbLinSweep_info{i}(rgbLinSweep_info{i}>1) = 1;

            % Trichromat (optimized) image for step i
            RGBTriCalFormat_opt = rgbLin2RGB(rgbLinSweep_info{i}, Disp);
            RGBTriImgFormat_opt = CalFormatToImage(RGBTriCalFormat_opt,m,n);
            RGBTriImgFormat_optAllInfo(:,:,:,i) = RGBTriImgFormat_opt;

            % Dichromat render of the optimized image for step i
            RGBDiCalFormat_opt = rgbLin2RGB(rgbLinRenderedSweep_info{i}, Disp);
            RGBDiImgFormat_opt = CalFormatToImage(RGBDiCalFormat_opt,m,n);
            RGBDiImgFormat_optAllInfo(:,:,:,i) = RGBDiImgFormat_opt;
        end

        % Figure stuff
        figName_info = sprintf('%s — %s — %d-step info sweep', imgType, dichromatType, nSteps);
        figInfo      = figure('Color','w','Name', figName_info);
        tlInfo       = tiledlayout(figInfo, 1, 2, 'TileSpacing','compact', 'Padding','compact');

        % Left panel: Trichromat grid
        triMontageInfo = nexttile(tlInfo, 1);
        montage(RGBTriImgFormat_optAllInfo, 'Size', [nRows nCols], 'Parent', triMontageInfo,'BorderSize', 1, 'BackgroundColor', 'w');
        title(triMontageInfo, 'Trichromat (optimized)');

        % Right panel: Dichromat grid
        diMontageInfo = nexttile(tlInfo, 2);
        montage(RGBDiImgFormat_optAllInfo, 'Size', [nRows nCols], 'Parent', diMontageInfo,'BorderSize', 1, 'BackgroundColor', 'w');
        title(diMontageInfo, sprintf('Dichromat render (%s)', dichromatType));
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
            computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_dist);

        % Target vs achieved distortion (normalized)
        fprintf('\nDistortion sweep (n=%d): target vs achieved distortion\n', numel(targetDistNorm_dist));
        disp([targetDistNorm_dist(:), distNormAch_dist(:)]);

        % Image montage
        RGBTriImgFormat_optAllDist = zeros(m, n, 3, nSteps);
        RGBDiImgFormat_optAllDist  = zeros(m, n, 3, nSteps);

        for i = 1:nSteps
            % rgbLinSweep_dist{i}(rgbLinSweep_dist{i} < 0) = 0;
            % rgbLinSweep_dist{i}(rgbLinSweep_dist{i} > 1) = 1;

            % Trichromat (optimized) image for step i
            RGBTriCalFormat_opt = rgbLin2RGB(rgbLinSweep_dist{i}, Disp);
            RGBTriImgFormat_opt = CalFormatToImage(RGBTriCalFormat_opt,m,n);
            RGBTriImgFormat_optAllDist(:,:,:,i) = RGBTriImgFormat_opt;

            % Dichromat render of the optimized image for step i
            RGBDiCalFormat_opt = rgbLin2RGB(rgbLinRenderedSweep_dist{i}, Disp);
            RGBDiImgFormat_opt = CalFormatToImage(RGBDiCalFormat_opt,m,n);
            RGBDiImgFormat_optAllDist(:,:,:,i) = RGBDiImgFormat_opt;
        end

        % Figure stuff
        figNameDist = sprintf('%s — %s — %d-step distortion sweep', imgType, dichromatType, nSteps);
        figDist = figure('Color','w','Name', figNameDist);
        tlDist  = tiledlayout(figDist, 1, 2, 'TileSpacing','compact', 'Padding','compact');

        % Left panel: Trichromat grid
        triMontageDist = nexttile(tlDist, 1);
        montage(RGBTriImgFormat_optAllDist, 'Size', [nRows nCols], 'Parent', triMontageDist,'BorderSize', 1, 'BackgroundColor', 'w');
        title(triMontageDist, 'Trichromat (optimized)');

        % Right panel: Dichromat grid
        diMontageDist = nexttile(tlDist, 2);
        montage(RGBDiImgFormat_optAllDist, 'Size', [nRows nCols], 'Parent', diMontageDist,'BorderSize', 1, 'BackgroundColor', 'w');
        title(diMontageDist, sprintf('Dichromat render (%s)', dichromatType));

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
        plot(distNormAch_dist(idx), infoNormAch_dist(idx), 'ko-', 'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerSize',10 ,'DisplayName', 'Distortion sweep (maximize info)');
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


end
-