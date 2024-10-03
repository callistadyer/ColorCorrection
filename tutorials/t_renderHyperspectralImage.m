function [lmsImageCalFormatTrichromat,lmsModuledCalFormatTrichromat,lmsDichromImageCalFormat,lmsDichromModuledCalFormat,cone_mean_orig] = t_renderHyperspectralImage(image,renderType,bPLOTscatter,bScale,bMinMod)
% Demonstrate how to read, add cone directed info, and render for tri- and dichromat
%
% Syntax:
%   t_renderHyperspectralImage(image,renderType,bPLOTscatter,bScale,bMinMod)
%
% Description:
%
% Inputs:
%   image:        - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.                       
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   bPLOTscatter  - Boolean. Plot scatter plots of lms values (1 or 0)
%   bScale:       - Boolean. Scale the image values into display range (1
%                   or 0).  A good idea except for 'gray'.
%   bMinMod:      - Boolean. When modulating the image in L M or S cone
%                   dimension, do this separately for each pixel (0),
%                   maximizing the modulation within the display gamut.  Or
%                   (1) take the minimum of the modulations that are within
%                   gamut for all pixels.
%
% Outputs:
%   None
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

% History
%   07/30/2024  dhb, cmd  Initial go.
%   08/27/2024  dhb, cmd  It's working, cleaning up.

% Examples:
%{
t_renderHyperspectralImage([],'Deuteranopia'  ,0,1,0)          
t_renderHyperspectralImage('scene1.mat','Protanopia'  ,0,1,0)    
[lmsImageCalFormat,lmsModuledCalFormat] = t_renderHyperspectralImage('scene2.mat','Deuteranopia',0,1,0);    
t_renderHyperspectralImage('scene3.mat','Tritanopia'  ,0,1,0)    
t_renderHyperspectralImage('scene4.mat','Deuteranopia',0,1,0)  
t_renderHyperspectralImage('scene5.mat','Deuteranopia',0,1,0)    
t_renderHyperspectralImage('gray','Deuteranopia',0,0,0)           
%}

%% Load hyperspectral image data
[hyperspectralImage wls d P_monitor] = loadImage(image);

% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

% Get some monitor primaries
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

% Get cone responses for every pixel of the hyperspectral image
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Render lms image that a trichromat sees.  Convert from cal format to
% image as well.
lmsImageCalFormatTrichromat = T_cones*hyperspectralImageCalFormat;
[RGBImageCalFormat_trichromat,rgbLinImageCalFormat,scaleFactor] = LMS2RGBCalFormat(lmsImageCalFormatTrichromat,d,T_cones,P_monitor,m,n,bScale);
RGBImage_trichromat = CalFormatToImage(RGBImageCalFormat_trichromat,m,n);

% Get original image into rgb so you can maximize gamut contrast.
% CAN MODIFY LMS2RGB image to return linear rgb image as well, since it is
% computed on the way.
[rgbLinImageCalFormat2,scaleFactor] = LMS2rgbLinCalFormat(lmsImageCalFormatTrichromat,d,T_cones,P_monitor,m,n,bScale);
% Get modulation for isochromatic plate modulation.
% This function is taking rgbLinImageCalFormat2 and using MaximizeGamutContrast to determine how much we can move 
% in a given cone direction (specifically, the direction of the missing cone) without going out of gamut
lmsModulationImgFormat = getDichromatConfusionModulation(rgbLinImageCalFormat2,renderType,bMinMod,T_cones,P_monitor,scaleFactor,m,n,bScale);

% Create isochromatic plates
[RGBModulated lmsModuledCalFormatTrichromat] = isochromaticPlates(image,renderType,lmsModulationImgFormat,bScale, ...
    'verbose',true);

% Mean of absent cone in original image. Used for replacing that cone value in dichromat image 
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        cone_mean_orig = mean(lmsImageCalFormatTrichromat(2,:));
    case 'Protanopia'   % l cone deficiency
        cone_mean_orig = mean(lmsImageCalFormatTrichromat(1,:));
    case 'Tritanopia'   % s cone deficiency
        cone_mean_orig = mean(lmsImageCalFormatTrichromat(3,:));
end 

% Dichromat manipulation (push trichromat LMS image into this function to get out LMS of dichromat)  
lmsDichromImageCalFormat        = tri2dichromatLMSCalFormat(lmsImageCalFormatTrichromat,renderType,cone_mean_orig); % gray
lmsDichromModuledCalFormat      = tri2dichromatLMSCalFormat(lmsModuledCalFormatTrichromat,renderType,cone_mean_orig); %

% Dichromat LMS --> RGB
[RGBImage_dichromatCalFormat,scaleFactor_di]       = LMS2RGBCalFormat(lmsDichromImageCalFormat,d,T_cones,P_monitor,m,n,bScale); % no modulation
[RGBPlate_dichromatCalFormat,scaleFactor_di_plate] = LMS2RGBCalFormat(lmsDichromModuledCalFormat, d,T_cones,P_monitor,m,n,bScale);  % isochromatic plate 

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
    [rgbImage_dichromatCalFormat,scaleFactor_di] = LMS2rgbLinCalFormat(lmsDichromImageCalFormat,d,T_cones,P_monitor,m,n,bScale);
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