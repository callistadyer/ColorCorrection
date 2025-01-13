
function [triRGBImgFormatCorrected] = dichromatCorrection(img,renderType,bScale,bMinMod,nSquares)
% Uses PCA to move 3D trichromatic image into 2 dimensions in attempt to
% create an accessible image for a dichromat
%
% Syntax:
%   [RGBImage_dichromat] = dichromatCorrection(img,renderType,bScale,bMinMod)
%
% Description:
%
% Inputs:
%   img:          - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.                       
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   bScale:       - Boolean. Scale the image values into display range (1
%                   or 0).  A good idea except for 'gray'.
%   bMinMod:      - Boolean. When modulating the image in L M or S cone
%                   dimension, do this separately for each pixel (0),
%                   maximizing the modulation within the display gamut.  Or
%                   (1) take the minimum of the modulations that are within
%                   gamut for all pixels.
%
% Outputs:
%   RGBImage_dichromat:  - Transformed RGB image after PCA and scaling. Also replaced missing cone as done in other code 
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

% History
%   09/05/2024  cmd  Initial go.
%
% Examples:
%{
[RGBImage_dichromat] = dichromatCorrection('gray','Deuteranopia',0,0,10)
[RGBImage_dichromat] = dichromatCorrection('scene2.mat','Deuteranopia',1,0,10)
%}

% Close out any stray figures
close all;

% Get trichromatic (LMS) image
[triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(img,renderType,0,bScale,bMinMod,nSquares);    

% Color correction to aid dichromacy
[triLMScalFormatCorrected       triLMSmeans]                   = colorCorrection(triLMSCalFormat,renderType,Disp,bScale);   % Original image
[triLMScalFormatCorrected_plate triLMSmeans_plate]             = colorCorrection(triLMSCalFormat_plate,renderType,Disp,bScale); % Image with plate
% correctedLMS = K_opt_plate * D_mnew + T_mean_plate;


% Create RGB image from LMS  
%%%%% ORIGINAL IMAGE AND PLATE %%%%%
% Dichromat simulation of original image
[diRGBCalFormatOrig]        = LMS2RGBCalFormat(diLMSCalFormat, Disp,bScale);
[diRGBCalFormatOrig_plate]  = LMS2RGBCalFormat(diLMSCalFormat_plate, Disp,bScale);         % isochromatic plate 

% Trichromat simulation of original image
[triRGBcalFormatOrig]       = LMS2RGBCalFormat(triLMSCalFormat, Disp,bScale);
[triRGBcalFormatOrig_plate] = LMS2RGBCalFormat(triLMSCalFormat_plate, Disp,bScale);  % isochromatic plate 

%%%%% CORRECTED IMAGE AND PLATE %%%%%
% Corrected trichromat image
[triRGBcalFormatCorrected]        = LMS2RGBCalFormat(triLMScalFormatCorrected, Disp,bScale);
[triRGBcalFormatCorrected_plate]  = LMS2RGBCalFormat(triLMScalFormatCorrected_plate, Disp,bScale);       % isochromatic plate 

% Corrected dichromat image
diLMSCalFormatCorrected                  = tri2dichromatLMSCalFormat(triLMScalFormatCorrected,renderType,Disp,bScale); 
diRGBCalFormatCorrected                  = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp,bScale);
diLMSCalFormatCorrected_plate            = tri2dichromatLMSCalFormat(triLMScalFormatCorrected_plate,renderType,Disp,bScale);      % isochromatic plate 
diRGBCalFormatCorrected_plate            = LMS2RGBCalFormat(diLMSCalFormatCorrected_plate, Disp,bScale); % isochromatic plate 


% Transform from cal format to image for viewing
% original trichromat
triRGBImgFormatOrig              = CalFormatToImage(triRGBcalFormatOrig,Disp.m,Disp.n); % no modulation
triRGBImgFormatOrig_plate        = CalFormatToImage(triRGBcalFormatOrig_plate,Disp.m,Disp.n); % modulation

% corrected trichromat
triRGBImgFormatCorrected         = CalFormatToImage(triRGBcalFormatCorrected,Disp.m,Disp.n); % no modulation
triRGBImgFormatCorrected_plate   = CalFormatToImage(triRGBcalFormatCorrected_plate,Disp.m,Disp.n); % modulation

% original dichromat
diRGBImgFormatOrig               = CalFormatToImage(diRGBCalFormatOrig,Disp.m,Disp.n); % no modulation
diRGBImgFormatOrig_plate         = CalFormatToImage(diRGBCalFormatOrig_plate,Disp.m,Disp.n); % modulation

% corrected dichromat
diRGBImgFormatCorrected          = CalFormatToImage(diRGBCalFormatCorrected,Disp.m,Disp.n); % no modulation
diRGBImgFormatCorrected_plate    = CalFormatToImage(diRGBCalFormatCorrected_plate,Disp.m,Disp.n); % plate


figure('Position',[214         261         501        1076]);

% Create a tiled layout
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add images to the tiles
nexttile
imshow(triRGBImgFormatOrig);
title('trichromat');

nexttile
imshow(diRGBImgFormatOrig);
title('dichromat');

nexttile
imshow(triRGBImgFormatCorrected);
title('trichromat corrected');

nexttile
imshow(diRGBImgFormatCorrected);
title('dichromat corrected');

nexttile
imshow(triRGBImgFormatOrig_plate);
title('trichromat - plate');

nexttile
imshow(diRGBImgFormatOrig_plate);
title('dichromat - plate');

nexttile
imshow(triRGBImgFormatCorrected_plate);
title('trichromat corrected - plate');

nexttile
imshow(diRGBImgFormatCorrected_plate);
title('dichromat corrected - plate');

sgtitle('Left: Trichromat, Right: Dichromat')
end
