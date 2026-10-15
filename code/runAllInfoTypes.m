function runAllInfoTypes()
% Sweep info across multiple info metrics (methods) for a single image
%
% Syntax:
%   run_all_infoTypes()
%
% Description:
%   For an input image and observer type, this script loops over several
%   info metrics (e.g., regression, Wade, L–M difference). For each
%   method it:
%     1) Optimization with that info function and parameters
%     2) Runs an info sweep with nSteps
%     3) Creates figure of achieved vs target info and Distortion vs Target Info
%     4) Creates image montage of results:
%             Row 1: Trichromat transforms for each step
%             Row 2: Dichromat renderings for each step


imageType     = 'flower1.png';
dichromatType = 'Deuteranopia';
m = 64; n = 64;
nSteps  = 11;
epsInfo = 1e-6;

methodDefs = { ...
    struct('name','Regress L,M,S from (L and S)', ...
    'infoFcn', @computeInfo_regress, ...
    'infoParams', struct('predictingWhat','L,M,S', 'predictingFromWhat','L and S')), ...
    struct('name','Wade (ours)', ...
    'infoFcn', @computeInfo_Wade, ...
    'infoParams', struct()), ...
    struct('name','L–M difference', ...
    'infoFcn', @computeInfo_LMdifference, ...
    'infoParams', struct()) ...
    };

distortionFcn    = @computeDistortion_squared;
distortionParams = struct();
renderFcn        = @DichromRenderLinear;
renderParams     = struct();

projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');

setType = 1; clearFlag = 0;

% Where to save the image
[LMSCalFormat, ~, ~, ~, Disp, imgParams, pathName] = colorCorrectionGenerateImages( ...
    imageType, setType, m, n, dichromatType, clearFlag);
% pathName looks like: 'Deuteranopia/flower1.png/s1_m30_n30'

%% Sweep over all the methods you want to check
for k = 1:numel(methodDefs)
    fprintf('Method %d/%d: %s\n', k, numel(methodDefs), methodDefs{k}.name);

    % Daltonization object for this metric
    d = daltonize( ...
        methodDefs{k}.infoFcn, methodDefs{k}.infoParams, ...
        distortionFcn, distortionParams, ...
        renderFcn, renderParams, ...
        Disp);

    % Info sweep
    [~, rgbLinTri_s, ~, rgbLinDi_s, ~, targetInfoVec, infoN_s, distN_s] = d.computeSweep( ...
        LMSCalFormat, imgParams, dichromatType, nSteps, pathName, 'info', 'epsInfo', epsInfo);

    % Figure out what the folder name is
    metricFolder = buildMetricFolderName(methodDefs{k}.infoFcn, methodDefs{k}.infoParams, distortionFcn);
    runFolder    = sprintf('%dsteps', nSteps);
    saveDir      = fullfile(outputDir, 'testImagesTransformed', pathName, metricFolder, runFolder);
    if ~exist(saveDir,'dir'), mkdir(saveDir); end


    % Plotting:
    %%%%%%%%%% Figure 1: Plotting achieved vs target (Info & Distortion) %%%%%%%%%%
    xTarget = targetInfoVec(:).';
    achInfo = cellfun(@double, infoN_s(:)).';
    achDist = cellfun(@double, distN_s(:)).';

    figName = sprintf('%s — %s — %s', imageType, dichromatType, methodDefs{k}.name);
    figA = figure('Color','w','Name', [figName ' — Metrics']);
    TL = tiledlayout(figA,1,2,'TileSpacing','compact','Padding','compact');
     
    ax1 = nexttile(TL,1);
    plot(ax1, xTarget, achInfo, '-o', 'LineWidth', 1.5); hold(ax1,'on');
    xmin = min(xTarget); xmax = max(xTarget);
    plot(ax1, [xmin xmax],[xmin xmax],'k--','LineWidth',1); hold(ax1,'off');
    grid(ax1,'on'); axis(ax1,'square');
    xlabel(ax1,'Target info (normalized)'); ylabel(ax1,'Achieved info (normalized)');
    title(ax1,'Achieved vs Target Info');

    ax2 = nexttile(TL,2);
    plot(ax2, xTarget, achDist, '-o', 'LineWidth', 1.5);
    grid(ax2,'on'); axis(ax2,'square');
    xlabel(ax2,'Target info (normalized)'); ylabel(ax2,'Achieved distortion (normalized)');
    title(ax2,'Distortion vs Target Info');

    sgtitle(TL, figName, 'FontWeight','bold');

    pngA = sprintf('FIG_%02d_%s.png', k, regexprep(methodDefs{k}.name,'[^A-Za-z0-9]+','_'));
    saveas(figA, fullfile(saveDir, pngA));

    % %%%%%%%%%% Figure 2: montage of all steps (Tri on top, Di on bottom) %%%%%%%%%%
    figB = figure('Color','w','Name', [figName ' — Two-Panel Grids']);
    set(figB,'InvertHardcopy','off');

    % Two panels side-by-side
    leftRect  = [0.02 0.08 0.47 0.88];   % [x y w h] normalized region for Tri
    rightRect = [0.51 0.08 0.47 0.88];   % region for Di

    % Near-square grid for thumbnails
    nCols = ceil(sqrt(nSteps));
    nRows = ceil(nSteps / nCols);

    % Margins
    lm=0.04; rm=0.02; tm=0.06; bm=0.05; hs=0.015; vs=0.025;
    cellW = (1 - lm - rm - (nCols-1)*hs) / nCols;
    cellH = (1 - tm - bm - (nRows-1)*vs) / nRows;

    % Panel titles drawn as annotations (no UI panels, so saveas works)
    annotation(figB,'textbox',[leftRect(1) leftRect(2)+leftRect(4) leftRect(3) 0.04], ...
        'String','Trichromat','LineStyle','none','HorizontalAlignment','left','FontWeight','bold');
    annotation(figB,'textbox',[rightRect(1) rightRect(2)+rightRect(4) rightRect(3) 0.04], ...
        'String','Dichromat','LineStyle','none','HorizontalAlignment','left','FontWeight','bold');

    for j = 1:nSteps
        % Grid position (row top->bottom)
        r = ceil(j/nCols);
        c = j - (r-1)*nCols;
        x = lm + (c-1)*(cellW+hs);
        y = 1 - tm - r*cellH - (r-1)*vs;

        % Map cell position into figure coords for each panel
        posL = [leftRect(1)  + x*leftRect(3),  leftRect(2)  + y*leftRect(4),  cellW*leftRect(3),  cellH*leftRect(4)];
        posR = [rightRect(1) + x*rightRect(3), rightRect(2) + y*rightRect(4), cellW*rightRect(3), cellH*rightRect(4)];

        % Left panel: trichromat image 
        axL = axes('Parent', figB, 'Units','normalized', 'Position',posL, ...
                   'Color','w','XColor','none','YColor','none','Box','off');
        RGBcal = rgbLin2RGB(rgbLinTri_s{j}, Disp);      
        RGBimg = CalFormatToImage(RGBcal, n, m);       
        image(axL, RGBimg); axis(axL,'image','off');
        title(axL, sprintf('Step %d\ninfo=%.3g', j, xTarget(j)), 'FontSize', 9);

        % Right panel: dichromat image 
        axR = axes('Parent', figB, 'Units','normalized', 'Position',posR, ...
                   'Color','w','XColor','none','YColor','none','Box','off');
        RGBcalD = rgbLin2RGB(rgbLinDi_s{j}, Disp);
        RGBimgD = CalFormatToImage(RGBcalD, n, m);
        image(axR, RGBimgD); axis(axR,'image','off');
        title(axR, sprintf('Step %d', j), 'FontSize', 9);
    end

    % Big title across the whole figure
    annotation(figB,'textbox',[0 0.93 1 0.06], 'String',[figName ' — Sweep Grids'], ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'LineStyle','none','FontWeight','bold','FontSize',12);

    % Save next to sweepOutputs.mat
    gridsDir = fullfile(saveDir, 'MontageFigures'); if ~exist(gridsDir,'dir'), mkdir(gridsDir); end
    pngB = sprintf('GRIDS_%02d_%s.png', k, regexprep(methodDefs{k}.name,'[^A-Za-z0-9]+','_'));
    saveas(figB, fullfile(gridsDir, pngB));

   
end

fprintf('\nSaved under:\n%s\n', ...
    fullfile(outputDir, 'testImagesTransformed', pathName));
end


