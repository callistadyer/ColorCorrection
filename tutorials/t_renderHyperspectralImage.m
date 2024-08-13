
function t_renderHyperspectralImage(image,renderType)

% function t_renderHyperspectralImage(image,renderType)
% Demonstrate how to read and then render a hyperspectral image.

% Example images to use
% t_renderHyperspectralImage('scene1.mat','Protanopia')
% t_renderHyperspectralImage('scene2.mat','Deuteranopia')
% t_renderHyperspectralImage('scene3.mat','Tritanopia')
% t_renderHyperspectralImage('scene4.mat','Deuteranopia')
% t_renderHyperspectralImage('scene5.mat','Deuteranopia')

% History
%   07/30/2024  dhb, cmd  Initial go.


%% Load hyperspectral image data 
if isempty(image)
    % This image comes with ISETBio.
    scene = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(scene,'fov',2);
% else
%     % This will work if you are in the Brainard Lab and have the
%     % HyperspectralSceneTutorial folder on your lab dropbox path.
%     dropboxDirPath = localDropboxDir();
%     scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
%     load(fullfile(scenesDir, 'scene3.mat'), 'scene');
else
    dropboxDirPath = localDropboxDir();
    scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
    load(fullfile(scenesDir, image), 'scene');
end

% call function that creates isochromatic plates
[RGB_modulated lms_ModuledCalFormat] = isochromaticPlates(image,renderType,.0005);

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

dichromImageCalFormat = lmsImageCalFormat;
deuterMLScale = 0.65;
protoLMScale = 1/0.65;
tritanSMScale = 0.25;

switch (renderType)
    case 'Deuteranopia'
        lms_ModuledCalFormat(2,:)       =  deuterMLScale *l_cone; % replace M cones with L cone PLATE
        dichromImageCalFormat(2,:)      =  deuterMLScale *l_cone; % replace M cones with L cone
    case 'Protanopia'
        lms_ModuledCalFormat(1,:)       = protoLMScale*m_cone; % replace L cones with M cone PLATE
        dichromImageCalFormat(1,:)      = protoLMScale*m_cone; % replace L cones with M cone
    case 'Tritanopia'
        lms_ModuledCalFormat(3,:)       = tritanSMScale*(m_cone + deuterMLScale*l_cone)/2; % replace S cones with M cone PLATE
        dichromImageCalFormat(3,:)      = tritanSMScale*(m_cone + deuterMLScale*l_cone)/2; % replace S cones with M cone
end

% Get dichromat image for looking at
[RGBImage_dichromatCalFormat,scaleFactor_di]                = LMS2RGBimg(dichromImageCalFormat     ,d,T_cones,P_monitor,m,n); % no modulation
[RGBModulatedPlate_dichromatCalFormat,scaleFactor_di_plate] = LMS2RGBimg(lms_ModuledCalFormat,d,T_cones,P_monitor,m,n); % plate 

% convert to image from cal format
RGBImage_dichromat          = CalFormatToImage(RGBImage_dichromatCalFormat         ,m,n); % no modulation
RGBModulatedPlate_dichromat = CalFormatToImage(RGBModulatedPlate_dichromatCalFormat,m,n); % plate

% Get dichromat image for comparing lms values
[rgbImage_dichromatCalFormat,scaleFactor_di] = LMS2rgbLinimg(dichromImageCalFormat,d,T_cones,P_monitor,m,n);
rgbImage_dichromat = CalFormatToImage(rgbImage_dichromatCalFormat,m,n);

% Show the trichromatic image and the dichromatic image
figure('position',[ 896         364        1231         883]); 
subplot(2,2,1);
imshow(RGBImage_trichromat);    % TRICHROMAT
title('Trichromat rendering - no modulation','FontSize',20);

subplot(2,2,3);
imshow(RGB_modulated);          % TRICHROMAT PLATE
title('Trichromat rendering - plate','FontSize',20);

subplot(2,2,2);
imshow(RGBImage_dichromat);     % DICHROMAT
title([renderType ' rendering - no modulation'],'FontSize',20);

subplot(2,2,4);
imshow(RGBModulatedPlate_dichromat); % DICHROMAT PLATE
title([renderType ' rendering - plate'],'FontSize',20);

% Check by reversing RGB to LMS image
lmsImageDichromatFromrgb          = rgbLin2LMSimg(rgbImage_dichromat,T_cones,P_monitor,scaleFactor_di,m,n);
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