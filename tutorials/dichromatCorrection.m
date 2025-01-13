
function [RGBImage_dichromat] = dichromatCorrection(img,renderType,bScale,bMinMod,nSquares)
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
[lmsImageCalFormatTri,lmsModuledCalFormatTri,lmsDichromImageCalFormat,lmsDichromModuledCalFormat,cone_mean_orig,Disp,modDirection] = t_renderHyperspectralImage(img,renderType,0,bScale,bMinMod,nSquares);    

% Apply pca correction to aid dichromacy
[correctedLMS T_mean]                         = colorCorrectionPCA(img,lmsImageCalFormatTri,renderType,cone_mean_orig,Disp,bScale);   % Original image
[correctedLMS_plate T_mean_plate]             = colorCorrectionPCA(img,lmsModuledCalFormatTri,renderType,cone_mean_orig,Disp,bScale); % Image with plate
% correctedLMS = K_opt_plate * D_mnew + T_mean_plate;

rgbcheck = LMS2rgbLinCalFormat(lmsImageCalFormatTri,Disp,bScale);

% Scale corrected LMS values to be as close to possible to original LMS
% LMS_new       = correctedLMSadjust(correctedLMS,lmsImageCalFormatTri);
% LMS_new_plate = correctedLMSadjust(correctedLMS_plate,lmsModuledCalFormatTri);
% 
% correctedLMS       = LMS_new;
% correctedLMS_plate = LMS_new_plate;

% Create RGB image from LMS  
% Dichromat simulation of original image
[RGBImage_dichromatCalFormat_orig]        = LMS2RGBCalFormat(lmsDichromImageCalFormat, Disp,bScale);
[RGBImage_dichromatCalFormat_plate_orig]  = LMS2RGBCalFormat(lmsDichromModuledCalFormat, Disp,bScale);         % isochromatic plate 

% Trichromat simulation of original image
[RGBImage_trichromatCalFormat,scaleFactor_tri_plate] = LMS2RGBCalFormat(lmsImageCalFormatTri, Disp,bScale);
[RGBImage_trichromatCalFormat_plate,scaleFactor_tri] = LMS2RGBCalFormat(lmsModuledCalFormatTri, Disp,bScale);  % isochromatic plate 

% Corrected trichromat image via pca LMS values
[RGBImage_dichromatCalFormat,scaleFactor_di_plate]  = LMS2RGBCalFormat(correctedLMS, Disp,bScale);
[RGBImage_dichromatCalFormat_plate,scaleFactor_di]  = LMS2RGBCalFormat(correctedLMS_plate, Disp,bScale);       % isochromatic plate 

% Corrected dichromat image via pca LMS values
cone_mean_processed = mean(correctedLMS,2);
LMSfixedDichromat_plate                  = tri2dichromatLMSCalFormat(correctedLMS_plate,renderType,cone_mean_processed(2),Disp,bScale);      % isochromatic plate 
[RGBImage_fixedDichromatCalFormat_plate] = LMS2RGBCalFormat(LMSfixedDichromat_plate, Disp,bScale); % isochromatic plate 
LMSfixedDichromat                        = tri2dichromatLMSCalFormat(correctedLMS,renderType,cone_mean_processed(2),Disp,bScale); 
[RGBImage_fixedDichromatCalFormat]       = LMS2RGBCalFormat(LMSfixedDichromat, Disp,bScale);


% Transform from cal format to image for viewing
RGBImage_dichromat_orig          = CalFormatToImage(RGBImage_dichromatCalFormat_orig,Disp.m,Disp.n); % no modulation
RGBImage_dichromat_plate_orig    = CalFormatToImage(RGBImage_dichromatCalFormat_plate_orig,Disp.m,Disp.n); % modulation

RGBImage_dichromat               = CalFormatToImage(RGBImage_dichromatCalFormat,Disp.m,Disp.n); % no modulation
RGBImage_dichromat_plate         = CalFormatToImage(RGBImage_dichromatCalFormat_plate,Disp.m,Disp.n); % modulation

RGBImage_trichromat              = CalFormatToImage(RGBImage_trichromatCalFormat,Disp.m,Disp.n); % no modulation
RGBImage_trichromat_plate        = CalFormatToImage(RGBImage_trichromatCalFormat_plate,Disp.m,Disp.n); % modulation

RGBImage_FixedDichromat_plate    = CalFormatToImage(RGBImage_fixedDichromatCalFormat_plate,Disp.m,Disp.n); % plate
RGBImage_FixedDichromat          = CalFormatToImage(RGBImage_fixedDichromatCalFormat,Disp.m,Disp.n); % no modulation


figure('Position',[214         261         501        1076]);

% Create a tiled layout
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add images to the tiles
nexttile
imshow(RGBImage_trichromat);
title('trichromat');

nexttile
imshow(RGBImage_dichromat_orig);
title('dichromat');

nexttile
imshow(RGBImage_dichromat);
title('trichromat corrected - no plate');

nexttile
imshow(RGBImage_FixedDichromat);
title('dichromat corrected - no plate');

nexttile
imshow(RGBImage_trichromat_plate);
title('trichromat - plate');

nexttile
imshow(RGBImage_dichromat_plate_orig);
title('dichromat - plate');

nexttile
imshow(RGBImage_dichromat_plate);
title('trichromat corrected - plate');

nexttile
imshow(RGBImage_FixedDichromat_plate);
title('dichromat corrected - plate');

sgtitle('Left: Trichromat, Right: Dichromat')
end
