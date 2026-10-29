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
% imageTypes = {'flower1.png', 'fruit.png', 'Gaugin.png', 'ishi45.png'};   % all 4 test images (updates 1-3 at 64x64, 10/29)
imageTypes = {'square.png'};          % update 2 at 100x100 (10/29)

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
m = 60;
n = 60;

% THIS WAS FOR THE OLD NORMALIZING... WHERE IT WAS NO ABS() 
% m = 69;
% n = 69;

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
    'predictingWhat',     'M', ...          % options: 'M' or 'L-M' or 'L,M,S'
    'predictingFromWhat', 'L and S');           % options: 'L and S' or 'deltaL and deltaS' or 'deltaL+M and delta S'
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
    % Set up daltonizer object
    theDaltonizer = daltonize(infoFcn, infoParams,distortionFcn, distortionParams,renderFcn, renderParams, Disp);


    %% Let's do some R squared tests on this image - how much do the L and S cones explain without transforming at all?
    % Cone contrast of the original image (same as lossFunction)
    % C = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;
    % L = C(1,:)'; M = C(2,:)'; S = C(3,:)';
    % 
    % % Same regression as computeInfo_regressCIE: predict M from L and S, no intercept
    % X     = [L S];
    % w     = X \ M;
    % Mhat  = X * w;
    % resid = M - Mhat;
    % R2_noIntercept = 1 - sum(resid.^2) / sum(M.^2);
    % 
    % % Standard centered R^2 (with intercept), so the mean color can't inflate it
    % Xc     = [ones(size(L)) L S];
    % residC = M - Xc * (Xc \ M);
    % R2_centered = 1 - sum(residC.^2) / sum((M - mean(M)).^2);
    % 
    % fprintf('corr(L,M) = %.3f\n', corr(L, M));
    % fprintf('M = %.2f*L + %.2f*S\n', w(1), w(2));
    % fprintf('R2 (no intercept) = %.4f, R2 (centered) = %.4f\n', R2_noIntercept, R2_centered);
    % fprintf('rms of M = %.3f, rms of leftover = %.3f\n', sqrt(mean(M.^2)), sqrt(mean(resid.^2)));
    % 
    % 
    % figure; scatter(L, M, 8, 'filled',MarkerEdgeColor='black',MarkerFaceColor='white'); hold on; axis equal; grid on
    % plot(xlim, xlim, 'k--')
    % axis square
    % xlabel('L cone contrast'); 
    % ylabel('M cone contrast')
    % title(sprintf('corr(L,M) = %.3f', corr(L, M)))
    % 
    % lim = max(abs(M));
    % figure; tiledlayout(1,2)
    % nexttile; imagesc(reshape(M, m, n), [-lim lim]); axis image off; colorbar; title('M contrast')
    % nexttile; imagesc(reshape(resid, m, n), [-lim lim]); axis image off; colorbar; title('M not predictable from L,S')
    % colormap(gray)
    % 
    % Cpred = C; Cpred(2,:) = Mhat';
    % LMSpred = Cpred .* Disp.grayLMS + Disp.grayLMS;
    % 
    % Lab1 = XYZToLab(Disp.M_lms2xyz * LMSCalFormat, Disp.labWhiteXYZ);
    % Lab2 = XYZToLab(Disp.M_lms2xyz * LMSpred,      Disp.labWhiteXYZ);
    % dE   = ComputeDE2000_Lab(Lab1', Lab2');
    % 
    % toImg = @(LMS) CalFormatToImage(rgbLin2RGB(min(max(Disp.M_cones2rgb * LMS, 0), 1), Disp), m, n);
    % figure; tiledlayout(1,3)
    % nexttile; imshow(toImg(LMSCalFormat)); title('original')
    % nexttile; imshow(toImg(LMSpred));      title('M replaced by prediction from L,S')
    % nexttile; imagesc(reshape(dE, m, n)); axis image off; colorbar; title(sprintf('\\DeltaE2000, mean %.1f', mean(dE)))


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
            computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_dist, 'updateNumber', updateNumber);

        % Target vs achieved distortion (normalized)
        fprintf('\nDistortion sweep (n=%d): target vs achieved distortion\n', numel(targetDistNorm_dist));
        disp([targetDistNorm_dist(:), distNormAch_dist(:)]);

        % Image montage
        RGBTriImgFormat_optAllDist = zeros(m, n, 3, nSteps);
        RGBDiImgFormat_optAllDist  = zeros(m, n, 3, nSteps);

        for i = 1:nSteps
            rgbLinSweep_dist{i}(rgbLinSweep_dist{i} < 0) = 0;
            rgbLinSweep_dist{i}(rgbLinSweep_dist{i} > 1) = 1;

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




    %% Save individual sweep images to folders
    %
    % Saves each optimized image and its dichromat rendering as separate PNGs.
    % Creates a folder for each input image, dichromacy type, and sweep type.

    % % Remove extension from image name
    % [~, imgName, ~] = fileparts(imgType);
    % 
    % % Main output folder (inside current MATLAB working directory)
    % outputFolder = fullfile(pwd, 'ColorCorrection_Results', ...
    %     imgName, dichromatType);
    % 
    % if ~exist(outputFolder, 'dir')
    %     mkdir(outputFolder);
    % end

    % %% INFO SWEEP
    % if ranInfo
    % 
    %     infoFolder = fullfile(outputFolder, 'InfoSweep');
    % 
    %     if ~exist(infoFolder, 'dir')
    %         mkdir(infoFolder);
    %     end
    % 
    %     for i = 1:nSteps
    % 
    %         % Extract individual images
    %         imgTri = RGBTriImgFormat_optAllInfo(:,:,:,i);
    %         imgDi  = RGBDiImgFormat_optAllInfo(:,:,:,i);
    % 
    %         % Save trichromat optimized image
    %         imwrite(imgTri, fullfile(infoFolder, ...
    %             sprintf('Step_%02d_Trichromat.png', i)));
    % 
    %         % Save dichromat rendering
    %         imwrite(imgDi, fullfile(infoFolder, ...
    %             sprintf('Step_%02d_Dichromat.png', i)));
    % 
    %     end
    % 
    %     fprintf('Saved info sweep images to:\n%s\n', infoFolder);
    % 
    % end
    % 
    % 
    % %% DISTORTION SWEEP
    % if ranDist
    % 
    %     distFolder = fullfile(outputFolder, 'DistortionSweep');
    % 
    %     if ~exist(distFolder, 'dir')
    %         mkdir(distFolder);
    %     end
    % 
    %     for i = 1:nSteps
    % 
    %         % Extract individual images
    %         imgTri = RGBTriImgFormat_optAllDist(:,:,:,i);
    %         imgDi  = RGBDiImgFormat_optAllDist(:,:,:,i);
    % 
    %         % Save trichromat optimized image
    %         imwrite(imgTri, fullfile(distFolder, ...
    %             sprintf('Step_%02d_Trichromat.png', i)));
    % 
    %         % Save dichromat rendering
    %         imwrite(imgDi, fullfile(distFolder, ...
    %             sprintf('Step_%02d_Dichromat.png', i)));
    % 
    %     end
    % 
    %     fprintf('Saved distortion sweep images to:\n%s\n', distFolder);
    % 
    % end


    %% Plot each sweep as 2 rows x nSteps columns
    % Row 1 = Trichromat optimized images
    % Row 2 = Dichromat renderings
    % Above each column = Achieved information score

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFO SWEEP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ranInfo

        figure('Color','w', ...
            'Position',[100 100 150*nSteps 350], ...
            'Name',sprintf('%s - %s - Info Sweep',imgType,dichromatType));

        tl = tiledlayout(2,nSteps, ...
            'TileSpacing','none', ...
            'Padding','none');

        for i = 1:nSteps

            % Top row: Trichromat
            nexttile(i);
            imshow(RGBTriImgFormat_optAllInfo(:,:,:,i));
            axis image off;

            title(sprintf('%.3f',infoNormAch_info(i)), ...
                'FontSize',15);

            % Bottom row: Dichromat
            nexttile(nSteps+i);
            imshow(RGBDiImgFormat_optAllInfo(:,:,:,i));
            axis image off;

        end

        % sgtitle(sprintf('%s - %s - Info Sweep',imgType,dichromatType), ...
        %     'FontSize',14);

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISTORTION SWEEP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ranDist

        figure('Color','w', ...
            'Position',[100 100 150*nSteps 350], ...
            'Name',sprintf('%s - %s - Distortion Sweep',imgType,dichromatType));

        tl = tiledlayout(2,nSteps, ...
            'TileSpacing','none', ...
            'Padding','none');

        for i = 1:nSteps

            % Top row: Trichromat
            nexttile(i);
            imshow(RGBTriImgFormat_optAllDist(:,:,:,i));
            axis image off;

            title(sprintf('%.3f',infoNormAch_dist(i)), ...
                'FontSize',15);

            % Bottom row: Dichromat
            nexttile(nSteps+i);
            imshow(RGBDiImgFormat_optAllDist(:,:,:,i));
            axis image off;

        end

        % sgtitle(sprintf('%s - %s - Distortion Sweep',imgType,dichromatType), ...
        %     'FontSize',14);

    end

end

