function saveTransformedOutputs(outputs, pathName, nSteps, infoFcn, distortionFcn, Disp)
% saveTransformedOutputs  Save outputs, per-step images, and sweep montages.
%
% Inputs:
%   outputs        - Struct (with sweep cell array)
%   pathName       - Relative image path (e.g., 'Deuteranopia/flower/s1_m32_n32')
%   nSteps         - Number of sweep steps
%   isSweep        - True if this is a sweep, false if single output
%   infoFcn        - Handle to the info function
%   distortionFcn  - Handle to the distortion function
%   Disp           - Display struct (for rendering)

% Output folder structure
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');
saveBase    = fullfile(outputDir, 'testImagesTransformed');

infoFcnName       = func2str(infoFcn);
distortionFcnName = func2str(distortionFcn);
metricFolder      = sprintf('%s_%s', infoFcnName, distortionFcnName);
stepFolder        = sprintf('%dsteps', nSteps);

saveSubdir    = fullfile(saveBase, pathName, metricFolder, stepFolder);
imageSaveDir  = fullfile(saveSubdir, 'IndividualImages');
figureSaveDir = fullfile(saveSubdir, 'MontageFigures');

if ~exist(saveSubdir, 'dir'),    mkdir(saveSubdir); end
if ~exist(imageSaveDir, 'dir'),  mkdir(imageSaveDir); end
if ~exist(figureSaveDir, 'dir'), mkdir(figureSaveDir); end

% Save mat file
save(fullfile(saveSubdir, 'sweepOutputs.mat'), 'outputs');
fprintf('[saveTransformedOutputs] Saved sweep outputs to: %s\n', saveSubdir);

% Save image of every step
for i = 1:nSteps
    stepOut = outputs{i};  % Assume outputs is nSteps×1 cell array of structs

    infoVal = stepOut.targetInfoNormalized;

    % trichromat
    triRGB = CalFormatToImage(rgbLin2RGB(stepOut.rgbLinDaltonizedCalFormat, Disp), ...
                              stepOut.imgParams.m, stepOut.imgParams.n);
    triName = sprintf('step%02d_Info%.2f_trichromat.png', i, infoVal);
    imwrite(triRGB, fullfile(imageSaveDir, triName));

    % dichromat
    diRGB = CalFormatToImage(rgbLin2RGB(stepOut.rgbLinDaltonizedRenderedCalFormat, Disp), ...
                             stepOut.imgParams.m, stepOut.imgParams.n);
    diName = sprintf('step%02d_Info%.2f_dichromat.png', i, infoVal);
    imwrite(diRGB, fullfile(imageSaveDir, diName));
end

% Save montage figures
% Trichromat montage
figure('Name', 'Trichromat Sweep', 'Visible', 'off');
nCols = ceil(sqrt(nSteps));
nRows = ceil(nSteps / nCols);
for i = 1:nSteps
    subplot(nRows, nCols, i);
    triRGB = CalFormatToImage(rgbLin2RGB(outputs{i}.rgbLinDaltonizedCalFormat, Disp), ...
                              outputs{i}.imgParams.m, outputs{i}.imgParams.n);
    imshow(triRGB);
    title(LiteralUnderscore(sprintf('Info %.2f', outputs{i}.targetInfoNormalized)), 'FontSize', 8);
end
% sgtitle('Trichromat Sweep');
sgtitle(sprintf('Trichromat Sweep – %s, %s', infoFcnName, distortionFcnName));
saveas(gcf, fullfile(figureSaveDir, 'trichromatSweep.png'));
close(gcf);

% Dichromat montage
figure('Name', 'Dichromat Sweep', 'Visible', 'off');
for i = 1:nSteps
    subplot(nRows, nCols, i);
    diRGB = CalFormatToImage(rgbLin2RGB(outputs{i}.rgbLinDaltonizedRenderedCalFormat, Disp), ...
                             outputs{i}.imgParams.m, outputs{i}.imgParams.n);
    imshow(diRGB);
    title(LiteralUnderscore(sprintf('Info %.2f', outputs{i}.targetInfoNormalized)), 'FontSize', 8);
end
% sgtitle('Dichromat Sweep');
sgtitle(LiteralUnderscore(sprintf('Dichromat Sweep – %s, %s', infoFcnName, distortionFcnName)));
saveas(gcf, fullfile(figureSaveDir, 'dichromatSweep.png'));
close(gcf);

end
