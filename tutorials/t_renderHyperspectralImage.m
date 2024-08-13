
function t_renderHyperspectralImage(image)

% function t_renderHyperspectralImage()
% Demonstrate how to read and then render a hyperspectral image.

% History
%   07/30/2024  dhb, cmd  Initial go.

%% Close all open figures
clear; close all

%% Load hyperspectral image data 
if isempty(image) || ~exist(image)
defaultImage = true;
end
if (defaultImage)
    % This image comes with ISETBio.
    scene = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(scene,'fov',2);
% else
%     % This will work if you are in the Brainard Lab and have the
%     % HyperspectralSceneTutorial folder on your lab dropbox path.
%     dropboxDirPath = localDropboxDir();
%     scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
%     load(fullfile(scenesDir, 'scene3.mat'), 'scene');
end

% Get the wavelength sampling and the actual hyperspectral image data in energy units.
wls = sceneGet(scene,'wave');
hyperspectralImage = sceneGet(scene,'energy');
clear scene

% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

% Get some monitor primaries
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

% Get cone responses for every pixel of the hyperspectral image
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Render lms image that a trichromat or dichromat sees
lmsImageCalFormat = T_cones*hyperspectralImageCalFormat;
[RGBImageCalFormat_trichromat,scaleFactor_tri] = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n);

% Image format
RGBImage_trichromat = CalFormatToImage(RGBImageCalFormat_trichromat,m,n);

%%% Testing different constant values of m-cone (or other missing cone)
% alpha        = linspace(0,1,20);
% mCones       = lmsImageCalFormat(2,:);
% contrast     = (mCones - mean(mCones)) ./ mean(mCones);
% contrast_new = alpha' .* contrast;
% m_new        = (contrast_new .* mean(mCones)) + mean(mCones);
l_cone = lmsImageCalFormat(1,:);
m_cone = lmsImageCalFormat(2,:);
s_cone = lmsImageCalFormat(3,:);

renderType = 'Protanopia';
dichromImageCalFormat = lmsImageCalFormat;
deuterMLScale = 0.65;
protoLMScale = 1/0.65;
tritanSMScale = 0.25;
switch (renderType)
    case 'Deuteranopia'
        dichromImageCalFormat(2,:)  =  deuterMLScale *l_cone; % replace M cones with L cone
    case 'Protanopia'
        dichromImageCalFormat(1,:)  = protoLMScale*m_cone; % replace L cones with M cone
    case 'Tritanopia'
        dichromImageCalFormat(3,:)  = tritanSMScale*(m_cone + deuterMLScale*l_cone)/2; % replace S cones with M cone
end

% Get dichromat image for looking at
[RGBImage_dichromatCalFormat,scaleFactor_di,rgbImage_dichromat] = LMS2RGBimg(dichromImageCalFormat,d,T_cones,P_monitor,m,n);
RGBImage_dichromat = CalFormatToImage(RGBImage_dichromatCalFormat,m,n);
% Show the trichromatic image and the dichromatic image
figure('position',[896         896        1152         363]); 
subplot(1,2,1);
imshow(RGBImage_trichromat);    % TRICHROMAT
title('Trichromat rendering');

subplot(1,2,2);
imshow(RGBImage_dichromat);     % DICHROMAT
title([renderType ' rendering']);

% Check by reversing RGB to LMS image
lmsImageDichromatFromrgb          = rgb2LMSimg(rgbImage_dichromat,T_cones,P_monitor,scaleFactor_di,m,n);
lmsImageDichromatFromrgbCalFormat = ImageToCalFormat(lmsImageDichromatFromrgb);
lmsImageDichromatFromRGB          = RGB2LMSimg(RGBImage_dichromat,d,T_cones,P_monitor,scaleFactor_di,m,n); 
lmsImageDichromatFromRGBCalFormat = ImageToCalFormat(lmsImageDichromatFromRGB);


%%% SCATTER PLOTS OF DESIRED VS RECOVERED LMS VALUES %%%
switch (renderType)
    case 'Deuteranopia'
        figure('position',[743         503        1200         1200]);
        subplot(2,2,1)
        scatter(l_cone,lmsImageDichromatFromRGBCalFormat(1,:),'red','Marker','.','LineWidth',2)
        xlabel('desired L'); ylabel('recovered from RGB L'); title('L from RGB');
        axis square;
        
        subplot(2,2,2)
        scatter(s_cone,lmsImageDichromatFromRGBCalFormat(3,:),'blue','Marker','.','LineWidth',2)
        xlabel('desired S'); ylabel('recovered from RGB S'); title('S from RGB');
        axis square;

        subplot(2,2,3)
        scatter(l_cone,lmsImageDichromatFromrgbCalFormat(1,:),'red','Marker','.','LineWidth',2)
        xlabel('desired L'); ylabel('recovered from rgb L'); title('L from rgb');
        axis square;
        
        subplot(2,2,4)
        scatter(s_cone,lmsImageDichromatFromrgbCalFormat(3,:),'blue','Marker','.','LineWidth',2)
        xlabel('desired S'); ylabel('recovered from rgb S'); title('S from rgb');
        axis square;

        sgtitle([renderType ' rendering and reversal: L and S values'])
    case 'Protanopia'
        figure('position',[743         503        1200         1200]);
        subplot(2,2,1)
        scatter(m_cone,lmsImageDichromatFromRGBCalFormat(2,:),'green','Marker','.','LineWidth',2)
        xlabel('desired M'); ylabel('recovered from RGB M'); title('M from RGB');
        axis square;
        
        subplot(2,2,2)
        scatter(s_cone,lmsImageDichromatFromRGBCalFormat(3,:),'blue','Marker','.','LineWidth',2)
        xlabel('desired S'); ylabel('recovered from RGB S'); title('S from RGB');
        axis square;

        subplot(2,2,3)
        scatter(m_cone,lmsImageDichromatFromrgbCalFormat(2,:),'green','Marker','.','LineWidth',2)
        xlabel('desired M'); ylabel('recovered from rgb M'); title('M from rgb');
        axis square;
        
        subplot(2,2,4)
        scatter(s_cone,lmsImageDichromatFromrgbCalFormat(3,:),'blue','Marker','.','LineWidth',2)
        xlabel('desired S'); ylabel('recovered from rgb S'); title('S from rgb');
        axis square;

        sgtitle([renderType ' rendering and reversal: M and S values'])
    case 'Tritanopia'
        figure('position',[743         503        1200         1200]);
        subplot(2,2,1)
        scatter(l_cone,lmsImageDichromatFromRGBCalFormat(1,:),'red','Marker','.','LineWidth',2)
        xlabel('desired L'); ylabel('recovered from RGB L'); title('L from RGB');
        axis square;
        
        subplot(2,2,2)
        scatter(m_cone,lmsImageDichromatFromRGBCalFormat(2,:),'green','Marker','.','LineWidth',2)
        xlabel('desired M'); ylabel('recovered from RGB M'); title('M from RGB');
        axis square;

        subplot(2,2,3)
        scatter(l_cone,lmsImageDichromatFromrgbCalFormat(1,:),'red','Marker','.','LineWidth',2)
        xlabel('desired L'); ylabel('recovered from rgb L'); title('L from rgb');
        axis square;
        
        subplot(2,2,4)
        scatter(m_cone,lmsImageDichromatFromrgbCalFormat(2,:),'green','Marker','.','LineWidth',2)
        xlabel('desired M'); ylabel('recovered from rgb M'); title('M from rgb');
        axis square;

        sgtitle([renderType ' rendering and reversal: L and M values'])
end


end

%%%%%%% GO FROM LMS IMAGE TO RGB IMAGE %%%%%%%
function [RGBImageCalFormat,scaleFactor,rgbImage] = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n)

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
rgbImageCalFormat = M_cones2rgb*lmsImageCalFormat;
rgbImage = CalFormatToImage(rgbImageCalFormat,m,n);

% For right now, normalize so that maximum value in rgb is 1
scaleFactor = max(rgbImage(:)); % save scale factor for later 
rgbImage = rgbImage/scaleFactor;

% Truncated version for gamma correction
rgbImageTruncate = rgbImage;
rgbImageTruncate(rgbImageTruncate < 0) = 0;

% Gamma correct
iGtable = displayGet(d,'inversegamma');
RGBImage = rgb2dac(rgbImageTruncate,iGtable)/(2^displayGet(d,'dacsize')-1);

% Transform to cal format
RGBImageCalFormat = ImageToCalFormat(RGBImage);

end

%%%%%%% GO FROM RGB IMAGE TO LMS IMAGE %%%%%%%
function lmsImage = RGB2LMSimg(RGBImage,d,T_cones,P_monitor,scaleFactor,m,n)

% Reverse the gamma correction
gammaTable = displayGet(d,'gammatable');
rgbImage = dac2rgb(RGBImage, gammaTable)*(2^displayGet(d,'dacsize')-1);

% Undo scaling and convert to LMS
lmsImage = rgb2LMSimg(rgbImage,T_cones,P_monitor,scaleFactor,m,n);

end

function lmsImage = rgb2LMSimg(rgbImage,T_cones,P_monitor,scaleFactor,m,n)

% Undo the scaling 
rgbImage = rgbImage * scaleFactor;

% Cal format
rgbImageCalFormat = ImageToCalFormat(rgbImage);

% lms image
M_rgb2cones = T_cones*P_monitor;
lmsImageCalFormat = M_rgb2cones * rgbImageCalFormat;

lmsImage = CalFormatToImage(lmsImageCalFormat,m,n);
end



