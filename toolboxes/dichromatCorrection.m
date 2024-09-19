
function [RGBImage_dichromat] = dichromatCorrection(img,renderType,bScale,bMinMod)
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
[RGBImage_dichromat] = dichromatCorrection('gray','Deuteranopia',0,0)
[RGBImage_dichromat] = dichromatCorrection('scene2.mat','Deuteranopia',1,0)
%}

% Get trichromatic (LMS) image
[lmsImageCalFormatTri,lmsModuledCalFormatTri,lmsDichromImageCalFormat,lmsDichromModuledCalFormat,cone_mean_orig] = t_renderHyperspectralImage(img,renderType,0,bScale,bMinMod);    

% Apply pca correction to aid dichromacy
[correctedLMS K_opt D_mnew T_mean]                         = colorCorrectionPCA(img,lmsImageCalFormatTri,renderType,cone_mean_orig);   % Original image
[correctedLMS_plate K_opt_plate D_mnew_plate T_mean_plate] = colorCorrectionPCA(img,lmsModuledCalFormatTri,renderType,cone_mean_orig); % Image with plate
correctedLMS = K_opt_plate * D_mnew + T_mean_plate;

% Plot new image 
[hyperspectralImage wls d P_monitor] = loadImage(img);
% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);
% Do this to get m n dimensions
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);
P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

% Create RGB image from LMS  

% Dichromat simulation of original image
[RGBImage_dichromatCalFormat_orig]        = LMS2RGBimg(lmsDichromImageCalFormat, d,T_cones,P_monitor,m,n,bScale);
[RGBImage_dichromatCalFormat_plate_orig]  = LMS2RGBimg(lmsDichromModuledCalFormat, d,T_cones,P_monitor,m,n,bScale);         % isochromatic plate 

% Trichromat simulation of original image
[RGBImage_trichromatCalFormat,scaleFactor_tri_plate] = LMS2RGBimg(lmsImageCalFormatTri, d,T_cones,P_monitor,m,n,bScale);
[RGBImage_trichromatCalFormat_plate,scaleFactor_tri] = LMS2RGBimg(lmsModuledCalFormatTri, d,T_cones,P_monitor,m,n,bScale);  % isochromatic plate 

% Corrected trichromat image via pca LMS values
[RGBImage_dichromatCalFormat,scaleFactor_di_plate]  = LMS2RGBimg(correctedLMS, d,T_cones,P_monitor,m,n,bScale);
[RGBImage_dichromatCalFormat_plate,scaleFactor_di]  = LMS2RGBimg(correctedLMS_plate, d,T_cones,P_monitor,m,n,bScale);       % isochromatic plate 

% Corrected dichromat image via pca LMS values
cone_mean_processed = mean(correctedLMS,2);
LMSfixedDichromat_plate                  = tri2dichromatLMS(correctedLMS_plate,renderType,cone_mean_processed(2));      % isochromatic plate 
[RGBImage_fixedDichromatCalFormat_plate] = LMS2RGBimg(LMSfixedDichromat_plate, d,T_cones,P_monitor,m,n,bScale); % isochromatic plate 
LMSfixedDichromat                        = tri2dichromatLMS(correctedLMS,renderType,cone_mean_processed(2)); 
[RGBImage_fixedDichromatCalFormat]       = LMS2RGBimg(LMSfixedDichromat, d,T_cones,P_monitor,m,n,bScale);


% Transform from cal format to image for viewing
RGBImage_dichromat_orig          = CalFormatToImage(RGBImage_dichromatCalFormat_orig,m,n); % no modulation
RGBImage_dichromat_plate_orig    = CalFormatToImage(RGBImage_dichromatCalFormat_plate_orig,m,n); % modulation

RGBImage_dichromat               = CalFormatToImage(RGBImage_dichromatCalFormat,m,n); % no modulation
RGBImage_dichromat_plate         = CalFormatToImage(RGBImage_dichromatCalFormat_plate,m,n); % modulation

RGBImage_trichromat              = CalFormatToImage(RGBImage_trichromatCalFormat,m,n); % no modulation
RGBImage_trichromat_plate        = CalFormatToImage(RGBImage_trichromatCalFormat_plate,m,n); % modulation

RGBImage_FixedDichromat_plate    = CalFormatToImage(RGBImage_fixedDichromatCalFormat_plate,m,n); % plate
RGBImage_FixedDichromat          = CalFormatToImage(RGBImage_fixedDichromatCalFormat,m,n); % no modulation


figure();
subplot(4,2,1)
imshow(RGBImage_trichromat);  
title('trichromat')

subplot(4,2,2)
imshow(RGBImage_trichromat_plate);     
title('trichromat plate')

subplot(4,2,3)
imshow(RGBImage_dichromat_orig);     % DICHROMAT
title('dichromat')

subplot(4,2,4)
imshow(RGBImage_dichromat_plate_orig);
title('dichromat - plate')

subplot(4,2,5)
imshow(RGBImage_dichromat);     % DICHROMAT
title('trichromat corrected')

subplot(4,2,6)
imshow(RGBImage_dichromat_plate);
title('trichromat corrected - plate')

subplot(4,2,7)
imshow(RGBImage_FixedDichromat);
title('dichromat corrected')

subplot(4,2,8)
imshow(RGBImage_FixedDichromat_plate);
title('dichromat corrected - plate')

end
