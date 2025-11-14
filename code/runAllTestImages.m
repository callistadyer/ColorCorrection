function runAllTestImages()
%% List of image types
imageTypes = {'flower1.png'};

% imageTypes = {'flower1.png', ...
%     'fruit.png', ...
%     'Gaugin.png', ...
%     'ishi45.png'};

setType       = 1;
dichromatType = 'Deuteranopia';
clearFlag     = 0;
m = 31;
n = 31;

%% Loop over all image types
for ii = 1:numel(imageTypes)

    imgType = imageTypes{ii};
    fprintf('Processing image type: %s\n', imgType);

    % Generate input image
    [LMSCalFormat, rgbLinCalFormat, LMSCalFormatRendered, rgbLinCalFormatRendered, ...
        Disp, imgParams, pathName] = colorCorrectionGenerateImages( ...
        imgType, setType, m, n, dichromatType, clearFlag);

    %% Define objective functions
    % infoFcn = @computeInfo_LMdifference;
    infoFcn = @computeInfo_regress;
    % infoFcn = @computeInfo_Wade;
    infoParams = struct( ...
        'predictingWhat',     'L,M,S', ...          % options: 'M' or 'L-M' or 'L,M,S'
        'predictingFromWhat', 'L and S');           % options: 'L and S' or 'deltaL and deltaS' or 'deltaL+M and delta S'
    % infoParams = struct();
    distortionFcn = @computeDistortion_squared;
    distortionParams = struct();

    renderFcn = @DichromRenderLinear;
    renderParams = struct();

    %% Set up daltonizer object
    theDaltonizer = daltonize( ...
        infoFcn, infoParams, ...
        distortionFcn, distortionParams, ...
        renderFcn, renderParams, ...
        Disp);

%% INFO AND DISTORTION SWEEPS

%   (1) Run an info sweep (minimize distortion at each target info)
%   (2) Run a distortion sweep (maximize info at each target distortion)

nSteps = 11;                                % number of sweep steps
nCols  = ceil(sqrt(nSteps));                % grid columns (near-square)
nRows  = ceil(nSteps / nCols);              % grid rows

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) INFO SWEEP — minimize distortion subject to target info
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sweepAxis_info = 'info';

[LMSSweep_info, rgbLinSweep_info, ...
 LMSRenderedSweep_info, rgbLinRenderedSweep_info, ...
 TSweep_info, targetInfoNorm_info, targetDistNorm_info, ...
 infoNormAch_info, distNormAch_info] = computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType,nSteps, pathName, sweepAxis_info);

% Target vs achieved info (normalized)
fprintf('\nInfo sweep (n=%d): target vs achieved info\n', numel(targetInfoNorm_info));
disp([targetInfoNorm_info(:), infoNormAch_info(:)]);

% Image montage
triStack_info = zeros(m, n, 3, nSteps);
diStack_info  = zeros(m, n, 3, nSteps);

for i = 1:nSteps
    % Trichromat (optimized) image for step i
    rgbTri = rgbLin2RGB(rgbLinSweep_info{i}, Disp);
    imgTri = CalFormatToImage(rgbTri,m,n);
    imgTri = min(max(imgTri,0),1);       
    triStack_info(:,:,:,i) = imgTri;

    % Dichromat render of the optimized image for step i
    rgbDi = rgbLin2RGB(rgbLinRenderedSweep_info{i}, Disp);
    imgDi = CalFormatToImage(rgbDi,m,n);
    imgDi = min(max(imgDi,0),1);
    diStack_info(:,:,:,i) = imgDi;
end

figName_info = sprintf('%s — %s — %d-step INFO sweep (Tri vs Di)', imgType, dichromatType, nSteps);
figInfo = figure('Color','w','Name', figName_info);
tlInfo  = tiledlayout(figInfo, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% Left panel: Trichromat grid
axL_info = nexttile(tlInfo, 1);
montage(triStack_info, 'Size', [nRows nCols], 'Parent', axL_info, ...
    'BorderSize', 1, 'BackgroundColor', 'w');
title(axL_info, 'Trichromat (optimized)');

% Right panel: Dichromat grid
axR_info = nexttile(tlInfo, 2);
montage(diStack_info, 'Size', [nRows nCols], 'Parent', axR_info, ...
    'BorderSize', 1, 'BackgroundColor', 'w');
title(axR_info, sprintf('Dichromat render (%s)', dichromatType));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) DISTORTION SWEEP — maximize info subject to target distortion
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sweepAxis_dist = 'distortion';

[LMSSweep_dist, rgbLinSweep_dist, ...
 LMSRenderedSweep_dist, rgbLinRenderedSweep_dist, ...
 TSweep_dist, targetInfoNorm_dist, targetDistNorm_dist, ...
 infoNormAch_dist, distNormAch_dist] = computeSweep(theDaltonizer,LMSCalFormat, imgParams, dichromatType, nSteps, pathName, sweepAxis_dist);

% Target vs achieved distortion (normalized)
fprintf('\nDistortion sweep (n=%d): target vs achieved distortion\n', numel(targetDistNorm_dist));
disp([targetDistNorm_dist(:), distNormAch_dist(:)]);

% Image montage
triStack_dist = zeros(m, n, 3, nSteps);
diStack_dist  = zeros(m, n, 3, nSteps);

for i = 1:nSteps
    % Trichromat (optimized) image for step i
    rgbTri = rgbLin2RGB(rgbLinSweep_dist{i}, Disp);
    imgTri = CalFormatToImage(rgbTri,m,n);
    imgTri = min(max(imgTri,0),1);
    triStack_dist(:,:,:,i) = imgTri;

    % Dichromat render of the optimized image for step i
    rgbDi = rgbLin2RGB(rgbLinRenderedSweep_dist{i}, Disp);
    imgDi = CalFormatToImage(rgbDi,m,n);
    imgDi = min(max(imgDi,0),1);
    diStack_dist(:,:,:,i) = imgDi;
end

figName_dist = sprintf('%s — %s — %d-step DISTORTION sweep (Tri vs Di)', imgType, dichromatType, nSteps);
figDist = figure('Color','w','Name', figName_dist);
tlDist  = tiledlayout(figDist, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% Left panel: Trichromat grid
axL_dist = nexttile(tlDist, 1);
montage(triStack_dist, 'Size', [nRows nCols], 'Parent', axL_dist, ...
    'BorderSize', 1, 'BackgroundColor', 'w');
title(axL_dist, 'Trichromat (optimized)');

% Right panel: Dichromat grid
axR_dist = nexttile(tlDist, 2);
montage(diStack_dist, 'Size', [nRows nCols], 'Parent', axR_dist, ...
    'BorderSize', 1, 'BackgroundColor', 'w');
title(axR_dist, sprintf('Dichromat render (%s)', dichromatType));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlay achieved info vs achieved distortion (both sweeps)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Name','Achieved Info vs Distortion — INFO vs DIST sweeps');
plot(distNormAch_info, infoNormAch_info, 'o-','LineWidth',1.5, ...
    'DisplayName','Info sweep (minimize distortion)'); hold on;
plot(distNormAch_dist, infoNormAch_dist, 's-','LineWidth',1.5, ...
    'DisplayName','Distortion sweep (maximize info)');
grid on; axis square;
xlabel('Achieved Distortion (normalized)');
ylabel('Achieved Info (normalized)');
title(sprintf('%s — %s — %d steps (overlay)', imgType, dichromatType, nSteps));
legend('Location','southeast');
% xlim([0 1]); ylim([0 1]);

   

end

