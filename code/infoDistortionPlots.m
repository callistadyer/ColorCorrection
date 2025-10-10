function [infoNorm, distortionNorm] = infoDistortionPlots(dichromatType, imageName, m, n, infoName, distortionName, nSteps, Live, Save, bPLOT)
% infoDistortionPlots  Visualize the tradeoff between info and distortion across sweep steps
%
% Syntax:
% Description:
%   Loads precomputed sweep outputs (from computeInfoSweep) for a given image,
%   dichromat type, and info/distortion metric combination. Extracts normalized information and
%   distortion values for each sweep step, and generates plots to visualize Info vs. distortion
%
% Inputs:
%   dichromatType   - String specifying dichromacy type:
%                       'Protanopia'
%                       'Deuteranopia'
%                       'Tritanopia'
%   imageName       - Image string (e.g., 'flower1.png')
%   m, n            - Size of the image (e.g., 60, 60)
%   infoName        - Name of information metric used in optimization
%   distortionName  - Name of distortion metric used in optimization
%   nSteps          - Number of sweep steps (e.g., 30)
%   Live            - Do you want the plot to plot live?
%                           true = leave figures open
%                           false = close after plotting
%   Save            - Do you want to save the figure as png? 
%                           true = save info-vs-distortion plot as PNG
%
% Example:
%{
   infoDistortionPlots('Deuteranopia', 'Gaugin.png', 128, 128, ...
       'Wade', 'squared', 30, true, true);
%}

% Get the directory where we can find the outputs files with info and
% distortion values 
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');
rootDir     = fullfile(outputDir, 'testImagesTransformed');

% Load the sweep outputs
outputs = loadSweepOutputs(rootDir, dichromatType, imageName, m, n, infoName, distortionName, nSteps);
n = numel(outputs);
distortionNorm = nan(1,n);  
infoNorm = nan(1,n);  

% Loop through each sweep step and grab the info and distortion
for i = 1:n
    s = outputs{i};
    distortionNorm(i) = s.distortionNormalized;
    infoNorm(i) = s.infoNormalized;
end

% Only do this plotting to visualize the lambda case
% figure();
% subplot(1,2,1)
% plot(linspace(0,1,30),infoNorm,'-ko',LineWidth=2,MarkerSize=10)
% title('Info (norm) as a function of Lambda')
% xlabel('Lambda','FontSize',12)
% ylabel('Info Normalized','FontSize',12)
% ylim([0 (max(infoNorm)+(max(infoNorm)/8))])
% axis square
% subplot(1,2,2)
% plot(linspace(0,1,30),distortionNorm,'-ko',LineWidth=2,MarkerSize=10)
% title('Distortion (norm) as a function of Lambda')
% xlabel('Lambda','FontSize',12)
% ylabel('Distortion Normalized','FontSize',12)
% ylim([0 (max(distortionNorm)+(max(distortionNorm)/8))])
% axis square

if bPLOT == 1
if Live
    figVis = 'on';
else
    figVis = 'off';
end

% Plot the distortion vs info 
figure('Color','w');
plot(distortionNorm, infoNorm, '-o', 'LineWidth',2, 'MarkerSize',8, 'MarkerFaceColor','w', 'Color','k');
grid on; axis square;
xlabel('distortionNormalized');
ylabel('infoNormalized');

end
% Save the plots if you want
if Save
    outPath = fullfile(rootDir,'info_vs_distortion.png');
    try
        exportgraphics(gca, outPath, 'Resolution',200,'BackgroundColor','w');
    catch
        set(f,'PaperPositionMode','auto');
        print(f, outPath, '-dpng','-r200');
    end
    fprintf('[infoDistortionPlots] wrote: %s\n', outPath);
end



end
