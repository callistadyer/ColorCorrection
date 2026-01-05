% Image types
% imageTypes = {'fruit.png'};
% imageTypes = {'flower1.png', ...
%     'fruit.png', ...
%     'Gaugin.png', ...
%     'ishi45.png'};

imageTypes = {'flower1.png', ...
    'fruit.png', ...
    'Gaugin.png'};

% Set type (choose 1 unless ishihara)
setType       = 1;

% Dichromacy type
dichromatType = 'Deuteranopia';

% Image size (keep <60 for fast testing)
m = 55;
n = 55;

% How many steps do you want in a transformation sweep?
nSteps = 11;

sweepAxis = 'both';
%   'info'  = run info sweep only
%   'distortion'  = run distortion sweep only
%   'both'  = run both sweeps + overlay plot

% Do you want to clear the image transformation results if it already
% exists (0)? Or do you want to load the image result if it exists (1)?
clearFlag     = 0;

% Rerun a certain point
% Set enable=true to rerun one step
rerun.enable = true;        
rerun.step   = 3;         
rerun.which  = 'distortion';  % 'info' or 'distortion'
rerun.img    = 'flower1.png';  % rerun only for this image
if rerun.enable
    sweepAxis = rerun.which; 
end



%% Define objective functions
%%%%%%%%% INFO FUNCTION %%%%%%%%%
infoFcn = @computeInfo_regress;
% infoFcn = @computeInfo_LMdifference;
% infoFcn = @computeInfo_Wade;

% infoParams only necessary when infoFcn is regress
infoParams = struct( ...
    'predictingWhat',     'L,M,S', ...          % options: 'M' or 'L-M' or 'L,M,S'
    'predictingFromWhat', 'L and S');           % options: 'L and S' or 'deltaL and deltaS' or 'deltaL+M and delta S'
% infoParams = struct();

%%%%%%% DISTORTION FUNCTION %%%%%%
% distortionFcn = @computeDistortion_squared;
distortionFcn = @computeDistortion_DE2000;
distortionParams = struct();

%%%%%%%%% RENDER FUNCTION %%%%%%%%
renderFcn = @DichromRenderLinear;
renderParams = struct();

%% Loop over all image types
for ii = 1:numel(imageTypes)

    imgType = imageTypes{ii};
    fprintf('Processing image type: %s\n', imgType);

    % If rerunning, only process the requested image
    if rerun.enable && ~strcmpi(imgType, rerun.img)
        fprintf('  [rerun] skipping (want img=%s)\n', rerun.img);
        continue;
    end

    % Generate input image
    [LMSCalFormat, rgbLinCalFormat, LMSCalFormatRendered, rgbLinCalFormatRendered, ...
        Disp, imgParams, pathName] = colorCorrectionGenerateImages( ...
        imgType, setType, m, n, dichromatType, clearFlag);

    % Set up daltonizer object
    theDaltonizer = daltonize( ...
        infoFcn, infoParams, ...
        distortionFcn, distortionParams, ...
        renderFcn, renderParams, ...
        Disp);

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

        % [LMSSweep_info, rgbLinSweep_info, ...
        %     LMSRenderedSweep_info, rgbLinRenderedSweep_info, ...
        %     TSweep_info, targetInfoNorm_info, targetDistNorm_info, ...
        %     infoNormAch_info, distNormAch_info] = computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType,nSteps, pathName, sweepAxis_info);

        % If rerunning, tell it which step of the sweep to re-do
        extraArgs = {};
        if rerun.enable && strcmpi(rerun.which,'info')
            extraArgs = {'rerunStep', rerun.step};
        end

        [LMSSweep_info, rgbLinSweep_info, ...
            LMSRenderedSweep_info, rgbLinRenderedSweep_info, ...
            TSweep_info, targetInfoNorm_info, targetDistNorm_info, ...
            infoNormAch_info, distNormAch_info] = computeSweep( ...
            theDaltonizer, LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_info, ...
            extraArgs{:});

        % Target vs achieved info
        fprintf('\nInfo sweep (n=%d): target vs achieved info\n', numel(targetInfoNorm_info));
        disp([targetInfoNorm_info(:), infoNormAch_info(:)]);

        % Image montage
        rgbTriImgFormat_optAllInfo = zeros(m, n, 3, nSteps);
        rgbDiImgFormat_optAllInfo  = zeros(m, n, 3, nSteps);

        for i = 1:nSteps
            % Trichromat (optimized) image for step i
            rgbTriCalFormat_opt = rgbLin2RGB(rgbLinSweep_info{i}, Disp);
            rgbTriImgFormat_opt = CalFormatToImage(rgbTriCalFormat_opt,m,n);
            rgbTriImgFormat_optAllInfo(:,:,:,i) = rgbTriImgFormat_opt;

            % Dichromat render of the optimized image for step i
            rgbDiCalFormat_opt = rgbLin2RGB(rgbLinRenderedSweep_info{i}, Disp);
            rgbDiImgFormat_opt = CalFormatToImage(rgbDiCalFormat_opt,m,n);
            rgbDiImgFormat_optAllInfo(:,:,:,i) = rgbDiImgFormat_opt;
        end

        figName_info = sprintf('%s — %s — %d-step info sweep', imgType, dichromatType, nSteps);
        figInfo = figure('Color','w','Name', figName_info);
        tlInfo  = tiledlayout(figInfo, 1, 2, 'TileSpacing','compact', 'Padding','compact');

        % Left panel: Trichromat grid
        triMontageInfo = nexttile(tlInfo, 1);
        montage(rgbTriImgFormat_optAllInfo, 'Size', [nRows nCols], 'Parent', triMontageInfo,'BorderSize', 1, 'BackgroundColor', 'w');
        title(triMontageInfo, 'Trichromat (optimized)');

        % Right panel: Dichromat grid
        diMontageInfo = nexttile(tlInfo, 2);
        montage(rgbDiImgFormat_optAllInfo, 'Size', [nRows nCols], 'Parent', diMontageInfo,'BorderSize', 1, 'BackgroundColor', 'w');
        title(diMontageInfo, sprintf('Dichromat render (%s)', dichromatType));
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) DISTORTION SWEEP — maximize info subject to target distortion
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ranDist = false;
    if strcmp(sweepAxis,'distortion') || strcmp(sweepAxis,'both')
        ranDist = true;
        sweepAxis_dist = 'distortion';

        % [LMSSweep_dist, rgbLinSweep_dist, ...
        %     LMSRenderedSweep_dist, rgbLinRenderedSweep_dist, ...
        %     TSweep_dist, targetInfoNorm_dist, targetDistNorm_dist, ...
        %     infoNormAch_dist, distNormAch_dist] = computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_dist);

        % If rerunning, tell it which step of the sweep to re-do
        extraArgs = {};
        if rerun.enable && strcmpi(rerun.which,'distortion')
            extraArgs = {'rerunStep', rerun.step};
        end

        [LMSSweep_dist, rgbLinSweep_dist, ...
            LMSRenderedSweep_dist, rgbLinRenderedSweep_dist, ...
            TSweep_dist, targetInfoNorm_dist, targetDistNorm_dist, ...
            infoNormAch_dist, distNormAch_dist] = computeSweep( ...
            theDaltonizer, LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_dist, ...
            extraArgs{:});

        % Target vs achieved distortion (normalized)
        fprintf('\nDistortion sweep (n=%d): target vs achieved distortion\n', numel(targetDistNorm_dist));
        disp([targetDistNorm_dist(:), distNormAch_dist(:)]);

        % Image montage
        rgbTriImgFormat_optAllDist = zeros(m, n, 3, nSteps);
        rgbDiImgFormat_optAllDist  = zeros(m, n, 3, nSteps);

        for i = 1:nSteps
            % Trichromat (optimized) image for step i
            rgbTriCalFormat_opt = rgbLin2RGB(rgbLinSweep_dist{i}, Disp);
            rgbTriImgFormat_opt = CalFormatToImage(rgbTriCalFormat_opt,m,n);
            rgbTriImgFormat_optAllDist(:,:,:,i) = rgbTriImgFormat_opt;

            % Dichromat render of the optimized image for step i
            rgbDiCalFormat_opt = rgbLin2RGB(rgbLinRenderedSweep_dist{i}, Disp);
            rgbDiImgFormat_opt = CalFormatToImage(rgbDiCalFormat_opt,m,n);
            rgbDiImgFormat_optAllDist(:,:,:,i) = rgbDiImgFormat_opt;
        end

        figNameDist = sprintf('%s — %s — %d-step distortion sweep', imgType, dichromatType, nSteps);
        figDist = figure('Color','w','Name', figNameDist);
        tlDist  = tiledlayout(figDist, 1, 2, 'TileSpacing','compact', 'Padding','compact');

        % Left panel: Trichromat grid
        triMontageDist = nexttile(tlDist, 1);
        montage(rgbTriImgFormat_optAllDist, 'Size', [nRows nCols], 'Parent', triMontageDist,'BorderSize', 1, 'BackgroundColor', 'w');
        title(triMontageDist, 'Trichromat (optimized)');

        % Right panel: Dichromat grid
        diMontageDist = nexttile(tlDist, 2);
        montage(rgbDiImgFormat_optAllDist, 'Size', [nRows nCols], 'Parent', diMontageDist,'BorderSize', 1, 'BackgroundColor', 'w');
        title(diMontageDist, sprintf('Dichromat render (%s)', dichromatType));

    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Overlay achieved info vs achieved distortion (both sweeps)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure();
    % plot(distNormAch_info, infoNormAch_info, 'o-','LineWidth',1.5,'DisplayName','Info sweep (minimize distortion)'); hold on;
    % plot(distNormAch_dist, infoNormAch_dist, 's-','LineWidth',1.5,'DisplayName','Distortion sweep (maximize info)');
    % grid on;
    % axis square;
    % xlabel('Achieved Distortion (normalized)');
    % ylabel('Achieved Info (normalized)');
    % title(sprintf('%s — %s — %d steps', imgType, dichromatType, nSteps));
    % subtitle('Achieved Info vs Distortion for both sweeps')
    % legend('Location','southeast');

    figure();

    if ranInfo
         
        plot(distNormAch_info, infoNormAch_info, 'o-', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Info sweep (minimize distortion)');
        hold on;

        for k = 1:numel(targetInfoNorm_info)
            yline(targetInfoNorm_info(k), ':', ...
                'HandleVisibility','off');   % keeps legend clean
        end
        for k = 1:numel(targetDistNorm_info)
            xline(targetDistNorm_info(k), ':', ...
                'HandleVisibility','off');
        end

    end

    if ranDist

        plot(distNormAch_dist, infoNormAch_dist, 's-', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Distortion sweep (maximize info)');
        hold on;

        for k = 1:numel(targetInfoNorm_dist)
            yline(targetInfoNorm_info(k), ':', ...
                'HandleVisibility','off');   % keeps legend clean
        end
        for k = 1:numel(targetDistNorm_dist)
            xline(targetDistNorm_info(k), ':', ...
                'HandleVisibility','off');
        end

    end

    grid off; axis square;
    xlabel('Achieved Distortion (normalized)');
    ylabel('Achieved Info (normalized)');
    title(sprintf('%s — %s — %d steps', imgType, dichromatType, nSteps));

    if ranInfo && ranDist
        subtitle('Achieved Info vs Distortion for both sweeps');
    elseif ranInfo
        subtitle('Achieved Info vs Distortion (info sweep only)');
    else
        subtitle('Achieved Info vs Distortion (distortion sweep only)');
    end

    legend('Location','southeast');



end