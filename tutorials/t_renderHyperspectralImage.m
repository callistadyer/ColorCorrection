function [triLMScalFormat,triLMSCalFormat_plate,diLMScalFormat,diLMScalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(image,renderType,bPLOTscatter,bScale,nSquares,modType)
% Demonstrate how to read, add cone directed info, and render for tri- and dichromat
%
% Syntax:
%   [triLMScalFormat,triLMSCalFormat_plate,diLMScalFormat,diLMScalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(image,renderType,bPLOTscatter,bScale,nSquares,modType)
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
t_renderHyperspectralImage([],'Deuteranopia'  ,0,0,10,'rand')
t_renderHyperspectralImage('scene1.mat','Protanopia'  ,0,0,10)    
[lmsImageCalFormat,lmsModuledCalFormat] = t_renderHyperspectralImage('scene2.mat','Deuteranopia',0,1,0,10,'rand');    
t_renderHyperspectralImage('scene3.mat','Tritanopia'  ,0,0,10,'rand')    
t_renderHyperspectralImage('scene4.mat','Deuteranopia',0,0,10,'rand')  
t_renderHyperspectralImage('scene5.mat','Deuteranopia',0,0,10,'rand')    
t_renderHyperspectralImage('gray','Deuteranopia',0,0,10,'rand')           
%}

%% Load hyperspectral image data
[hyperspectralImage Disp] = loadImage(image);

% Get cone responses for every pixel of the hyperspectral image
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Render lms image that a trichromat sees.  Convert from cal format to image as well.
triLMScalFormat = Disp.T_cones*hyperspectralImageCalFormat;
[triRGBCalFormat,triRGBLinCalFormat,scaleFactor] = LMS2RGBCalFormat(triLMScalFormat,Disp,bScale);
triRGBImgFormat                                  = CalFormatToImage(triRGBCalFormat,m,n);


% Note: if you change modType to be "Deuteranopia" etc, then the following
% code will be overwritten inside of getDichromatConfusionModulation.m ...
% the vectors will be changed to the missing cone direction. otherwise,
% random color directions will be used

vectors = rand(nSquares, 3); % Uniform random values between 0 and 1
magnitudes = sqrt(sum(vectors.^2, 2)); % Compute the magnitudes
modDirection = vectors ./ magnitudes;  % Normalize to unit vectors

% Get modulation for isochromatic plate modulation.
% This function is taking rgbLinImageCalFormat2 and using MaximizeGamutContrast to determine how much we can move
% in a given cone direction (specifically, the direction of the missing cone) without going out of gamut
for i = 1:nSquares
    [lmsModulationImgFormat(:,:,:,i)] = getDichromatConfusionModulation(triRGBLinCalFormat,modDirection(i,:)',modType,Disp,scaleFactor,bScale);
end

% Create isochromatic plates
if ~strcmp(image,'74')
    [triRGBCalFormat_plate, triLMSCalFormat_plate] = isochromaticPlates(image,renderType,lmsModulationImgFormat,Disp,bScale,nSquares, ...
        'verbose',true);
    triRGBImgFormat_plate = CalFormatToImage(triRGBCalFormat_plate,Disp.m,Disp.n);
else
    % no need for plate when using plate 
    triLMSCalFormat_plate = triLMScalFormat;
    triRGBImgFormat_plate = triRGBImgFormat;
end
% Dichromat manipulation (push trichromat LMS image into this function to get out LMS of dichromat)  
diLMScalFormat        = tri2dichromatLMSCalFormat(triLMScalFormat,renderType,Disp,bScale); % gray
diLMScalFormat_plate  = tri2dichromatLMSCalFormat(triLMSCalFormat_plate,renderType,Disp,bScale); % modulated 

% Dichromat LMS --> RGB
[diRGBCalFormat,scaleFactor_di]             = LMS2RGBCalFormat(diLMScalFormat,Disp,bScale); % no modulation
[diRGBCalFormat_plate,scaleFactor_di_plate] = LMS2RGBCalFormat(diLMScalFormat_plate, Disp,bScale);  % isochromatic plate 

% Dichromat RGB cal format --> image format
diRGBImgFormat          = CalFormatToImage(diRGBCalFormat,m,n); % no modulation
diRGBImgFormat_plate    = CalFormatToImage(diRGBCalFormat_plate,m,n); % isochromatic plate

% Plotting
% Show the trichromatic image, dichromatic image, and trichromatic plate, dichromatic plate
% figure('position',[ 896         364        1231         883]); 
% subplot(2,2,1);
% imshow(triRGBImgFormat);    % TRICHROMAT  
% title('Trichromat rendering - no modulation','FontSize',20);
% 
% subplot(2,2,3);
% imshow(triRGBImgFormat_plate);          % TRICHROMAT PLATE
% title('Trichromat rendering - plate','FontSize',20);
% 
% subplot(2,2,2);
% imshow(diRGBImgFormat);     % DICHROMAT
% title([renderType ' rendering - no modulation'],'FontSize',20);
% 
% subplot(2,2,4);
% imshow(diRGBImgFormat_plate);     % DICHROMAT PLATE 
% title([renderType ' rendering - plate'],'FontSize',20);


if bPLOTscatter == 1

    % Get dichromat image for comparing lms values
    [rgbImage_dichromatCalFormat,scaleFactor_di] = LMS2rgbLinCalFormat(diLMScalFormat,d,T_cones,P_monitor,m,n,bScale);
    rgbImage_dichromat                           = CalFormatToImage(rgbImage_dichromatCalFormat,m,n);

    % Check by reversing RGB to LMS image
    lmsImageDichromatFromrgb          = rgbLin2LMSimg(rgbImage_dichromat,T_cones,P_monitor,scaleFactor_di,m,n,bScale);
    lmsImageDichromatFromrgbCalFormat = ImageToCalFormat(lmsImageDichromatFromrgb);
    lmsImageDichromatFromRGB          = RGB2LMSimg(diRGBImgFormat,d,T_cones,P_monitor,scaleFactor_di,m,n,bScale);
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