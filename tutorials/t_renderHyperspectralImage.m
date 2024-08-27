
function t_renderHyperspectralImage(image,renderType,bPLOTscatter,bScale,bMinMod)

% function t_renderHyperspectralImage(image,renderType)
% Demonstrate how to read and then render a hyperspectral image.

% Example images to use
%{
t_renderHyperspectralImage([],'Deuteranopia'  ,0,1,0)          
t_renderHyperspectralImage('scene1.mat','Protanopia'  ,0,1,0)    
t_renderHyperspectralImage('scene2.mat','Deuteranopia',0,1,0)    
t_renderHyperspectralImage('scene3.mat','Tritanopia'  ,0,1,0)    
t_renderHyperspectralImage('scene4.mat','Deuteranopia',0,1,0)  
t_renderHyperspectralImage('scene5.mat','Deuteranopia',0,1,0)    
t_renderHyperspectralImage('gray','Deuteranopia',0,0,0)           
%}

% Inputs
% image:        image to be rendered
% renderType:   type of dichromat
%               Deuteranopia
%               Protanopia
%               Tritanopia
% bPLOTscatter: plot scatter plots of lms values (1 or 0)
% bScale:       scale the image values (1 or 0)
% bMinMod:      when modulating the image in L M or S cone dimension, do this separately for each pixel (0)
%               or take the minimum modulation for all pixels (1)

% History
%   07/30/2024  dhb, cmd  Initial go.


%% Load hyperspectral image data 
if isempty(image)
    % This image comes with ISETBio.
    scene = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(scene,'fov',2);
    hyperspectralImage = sceneGet(scene,'energy');
elseif strcmp(image,'gray')
    % Grab a scene to define primaries
    scene = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(scene,'fov',2);

    % Get some monitor primaries
    wls = sceneGet(scene,'wave');
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

    % Create gray hyperspectral image
    % 256 x 256 gray image
    [grayImgCalFormat,m,n] = ImageToCalFormat(ones(256,256));
    % Gray 0.5 rgb at each pixel in image 
    grayImgCalFormat       = (0.5.*(repmat(grayImgCalFormat,3,1)));
    grayImgrgb             = CalFormatToImage(grayImgCalFormat,m,n);
    % Make hyperspectral img by multiplying primaries * rgb values at each pixel 
    hyperspecGrayCalFormat = P_monitor * grayImgCalFormat;
    % Image format
    hyperspectralImage     = CalFormatToImage(hyperspecGrayCalFormat,m,n);
else
    % This will work if you are in the Brainard Lab and have the
    % HyperspectralSceneTutorial folder on your lab dropbox path.
    dropboxDirPath = localDropboxDir();
    scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
    load(fullfile(scenesDir, image), 'scene');
    hyperspectralImage = sceneGet(scene,'energy');
end

% Get the wavelength sampling and the actual hyperspectral image data in energy units.
wls = sceneGet(scene,'wave');
clear scene

% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

% Get some monitor primaries
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

% Get cone responses for every pixel of the hyperspectral image
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Render lms image that a trichromat sees
lmsImageCalFormat = T_cones*hyperspectralImageCalFormat;
[RGBImageCalFormat_trichromat] = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n,bScale);

% RGB Image format
RGBImage_trichromat = CalFormatToImage(RGBImageCalFormat_trichromat,m,n);

% Get original image into rgb so you can maximize gamut contrast
[rgbImageCalFormat,scaleFactor] = LMS2rgbLinimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n,bScale);

% Get modulation for isochromatic plate modulation
lmsModulationImgFormat = getModulation(rgbImageCalFormat,renderType,bMinMod,T_cones,P_monitor,scaleFactor,m,n,bScale);

% Create isochromatic plates
[RGBModulated lmsModuledCalFormat] = isochromaticPlates(image,renderType,lmsModulationImgFormat,bScale);

% Mean of absent cone in original image. Used for replacing that cone value in dichromat image 
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        cone_mean_orig = mean(lmsModuledCalFormat(2,:));
    case 'Protanopia'   % l cone deficiency
        cone_mean_orig = mean(lmsModuledCalFormat(1,:));
    case 'Tritanopia'   % s cone deficiency
        cone_mean_orig = mean(lmsModuledCalFormat(3,:));
end 

% Dichromat manipulation (push trichromat LMS image into this function to get out LMS of dichromat)  
lmsDichromImageCalFormat = tri2dichromatLMS(lmsImageCalFormat,renderType,cone_mean_orig); % gray
lmsModuledCalFormat      = tri2dichromatLMS(lmsModuledCalFormat,renderType,cone_mean_orig); %

% Dichromat LMS --> RGB
[RGBImage_dichromatCalFormat,scaleFactor_di]       = LMS2RGBimg(lmsDichromImageCalFormat,d,T_cones,P_monitor,m,n,bScale); % no modulation
[RGBPlate_dichromatCalFormat,scaleFactor_di_plate] = LMS2RGBimg(lmsModuledCalFormat, d,T_cones,P_monitor,m,n,bScale);  % isochromatic plate 

% Dichromat RGB cal format --> image format
RGBImage_dichromat          = CalFormatToImage(RGBImage_dichromatCalFormat,m,n); % no modulation
RGBPlate_dichromat          = CalFormatToImage(RGBPlate_dichromatCalFormat,m,n); % isochromatic plate

% plotting
% Show the trichromatic image, dichromatic image, and trichromatic plate, dichromatic plate
figure('position',[ 896         364        1231         883]); 
subplot(2,2,1);
imshow(RGBImage_trichromat);    % TRICHROMAT  
title('Trichromat rendering - no modulation','FontSize',20);

subplot(2,2,3);
imshow(RGBModulated);          % TRICHROMAT PLATE
title('Trichromat rendering - plate','FontSize',20);

subplot(2,2,2);
imshow(RGBImage_dichromat);     % DICHROMAT
title([renderType ' rendering - no modulation'],'FontSize',20);

subplot(2,2,4);
imshow(RGBPlate_dichromat);     % DICHROMAT PLATE 
title([renderType ' rendering - plate'],'FontSize',20);





if bPLOTscatter == 1

    % Get dichromat image for comparing lms values
    [rgbImage_dichromatCalFormat,scaleFactor_di] = LMS2rgbLinimg(lmsDichromImageCalFormat,d,T_cones,P_monitor,m,n,bScale);
    rgbImage_dichromat                           = CalFormatToImage(rgbImage_dichromatCalFormat,m,n);

    % Check by reversing RGB to LMS image
    lmsImageDichromatFromrgb          = rgbLin2LMSimg(rgbImage_dichromat,T_cones,P_monitor,scaleFactor_di,m,n,bScale);
    lmsImageDichromatFromrgbCalFormat = ImageToCalFormat(lmsImageDichromatFromrgb);
    lmsImageDichromatFromRGB          = RGB2LMSimg(RGBImage_dichromat,d,T_cones,P_monitor,scaleFactor_di,m,n,bScale);
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


end