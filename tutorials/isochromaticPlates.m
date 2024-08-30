function [RGBmodulatedCalFormat lmsModuledCalFormat] = isochromaticPlates(img,renderType,lmsModulationImgFormat,bScale,options)

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
    bScale
    options.verbose (1,1) logical = false;
end

if (options.verbose)
    fprintf('Starting execution of isochromaticPlates\n');
end

% Load hyperspectral image
[hyperspectralImage wls d P_monitor] = loadImage(img);

% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

% Get LMS values
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);
lmsImageCalFormat = T_cones*hyperspectralImageCalFormat;
lmsImage          = CalFormatToImage(lmsImageCalFormat,m,n);

% get original image
RGB_CalFormat     = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n,bScale);
RGB_img           = CalFormatToImage(RGB_CalFormat,m,n);

switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        idxLMS = 2;
        slice_lms = lmsImage(:,:,idxLMS);
    case 'Protanopia'   % l cone deficiency
        idxLMS = 1;
        slice_lms = lmsImage(:,:,idxLMS);
    case 'Tritanopia'   % s cone deficiency
        idxLMS = 3;
        slice_lms = lmsImage(:,:,idxLMS);
end

% carve out an image where there is 0s everywhere except the delta 
% delta = plateO(size(slice_lms),deltaModulation,100);
% no longer just for one slice...
delta_lms = plateSquare(size(lmsImage),lmsModulationImgFormat,100);

% add the delta to the L M or S values to modulate one cone type
lmsImage_mod = lmsImage + delta_lms;
% redefine that L M or S cone values
% lmsImage(:,:,idxLMS) = new_slice;

% turn lmsImage to cal format so you can convert to RGB
lmsModuledCalFormat = ImageToCalFormat(lmsImage_mod);
% convert to RGB
RGB_modulatedCalFormat = LMS2RGBimg(lmsModuledCalFormat,d,T_cones,P_monitor,m,n,bScale);
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