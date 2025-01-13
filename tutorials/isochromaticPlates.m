function [RGBmodulatedCalFormat lmsModuledCalFormat] = isochromaticPlates(img,renderType,lmsModulationImgFormat,Disp,bScale,nSquares,options)

% function create isochromatic plates for testing dichromacy
%
% Syntax:
%   [RGBmodulatedCalFormat lmsModuledCalFormat] = isochromaticPlates(img,renderType,lmsModulationImgFormat,bScale,options)
%
% Description:
%
% Inputs:
%       img:                     Input image - which to create plate with?
%       renderType:              Type of dichromacy
%       lmsModulationImgFormat:  LMS plate modulation in img format
%       bScale:                  Scale rgb values or not?
%
% Outputs:
%       RGBmodulatedCalFormat    Gamma corrected RGB image with cone modulation included  
%       lmsModuledCalFormat      LMS values with cone modulation added 
%
% Optional key/value pairs:
%   None

% Examples:
%{
    [RGB_modulated lms_ModuledCalFormat] = isochromaticPlates('gray','Deuteranopia',.0005,0);
%}

%% Pick up optional arguments
arguments
    img
    renderType
    lmsModulationImgFormat
    Disp struct
    bScale
    nSquares
    options.verbose (1,1) logical = false;
end

if (options.verbose)
    fprintf('Starting execution of isochromaticPlates\n');
end

disp('Callista come back to this - make it so this func takes in LMS image')
% Load hyperspectral image
[hyperspectralImage Disp] = loadImage(img);
% Get LMS values
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);
lmsImageCalFormat = Disp.T_cones*hyperspectralImageCalFormat;
lmsImage          = CalFormatToImage(lmsImageCalFormat,Disp.m,Disp.n);

% Get original RGB image
[RGB_CalFormat rgbLinImageCalFormat]  = LMS2RGBCalFormat(lmsImageCalFormat,Disp,bScale);
RGB_img                               = CalFormatToImage(RGB_CalFormat,Disp.m,Disp.n);

% CREATE SQUARE MODULATIONS
delta_lms = plateSquare(size(lmsImage),lmsModulationImgFormat,nSquares);

% Add the delta to the L M S values to modulate cones (original LMS + modulation) 
lmsImage_mod = lmsImage + delta_lms;

% CHECK IF MODULATED LMS IS IN GAMUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmsImage_modCalFormat = ImageToCalFormat(lmsImage_mod);
inGamut = checkGamut(lmsImage_modCalFormat,Disp,bScale);
if inGamut == 0
    error(['isochromaticPlates: WARNING! rgb values are out of gamut... lmsImage_mod values outside of the range [0 1]']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% turn lmsImage to cal format so you can convert to RGB
lmsModuledCalFormat = ImageToCalFormat(lmsImage_mod);
% convert to RGB
RGB_modulatedCalFormat = LMS2RGBCalFormat(lmsModuledCalFormat,Disp,bScale);
% convert to image for viewing
RGBmodulatedCalFormat = CalFormatToImage(RGB_modulatedCalFormat,m,n);

figure('position',[927         886        1245         367]);
subplot(1,2,1)
imshow(RGB_img)
title('original image')

subplot(1,2,2)
imshow(RGBmodulatedCalFormat)
title([renderType ' testing plate'])

sgtitle('Isochromatic Plates')


end