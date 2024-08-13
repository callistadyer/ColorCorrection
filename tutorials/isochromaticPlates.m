function [RGB_modulated lms_ModuledCalFormat] = isochromaticPlates(img,renderType,deltaModulation)

% isochromaticPlates([],'Deuteranopia',.0005)
% scene = 0.5 * ones(256, 256);
% scene = sceneSet(scene,'fov',2);
% hyperspectralGrayImage = sceneGet(scene,'energy');

if isempty(img)
    sceneImg = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(sceneImg,'fov',2);
else
    dropboxDirPath = localDropboxDir();
    scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
    load(fullfile(scenesDir, img), 'scene');
end

wls = sceneGet(scene,'wave');
hyperspectralImage = sceneGet(scene,'energy');

% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

% Get some monitor primaries
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

% Get LMS values
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);
lmsImageCalFormat = T_cones*hyperspectralImageCalFormat;
lmsImage          = CalFormatToImage(lmsImageCalFormat,m,n);

% get original image
RGB_CalFormat     = LMS2RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n);
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
delta = plateO(size(slice_lms),deltaModulation,200);
% add the delta to the L M or S values to modulate one cone type
new_slice = slice_lms + delta;
% redefine that L M or S cone values
lmsImage(:,:,idxLMS) = new_slice;

% turn lmsImage to cal format so you can convert to RGB
lms_ModuledCalFormat = ImageToCalFormat(lmsImage);
% convert to RGB
RGB_modulatedCalFormat = LMS2RGBimg(lms_ModuledCalFormat,d,T_cones,P_monitor,m,n);
% convert to image for viewing
RGB_modulated = CalFormatToImage(RGB_modulatedCalFormat,m,n);

figure('position',[927         886        1245         367]);
subplot(1,2,1)
imshow(RGB_img)
title('original image')

subplot(1,2,2)
imshow(RGB_modulated)
title([renderType ' testing plate'])

sgtitle('Isochromatic Plates')


end