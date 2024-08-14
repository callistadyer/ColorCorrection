function [RGB_modulated lms_ModuledCalFormat] = isochromaticPlates(img,renderType,deltaModulation)

% function create isochromatic plates for testing dichromacy
% Example call: [RGB_modulated lms_ModuledCalFormat] = isochromaticPlates([],'Deuteranopia',.0005)

if isempty(img)
    % Default scenes
    sceneImg = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(sceneImg,'fov',2);
    hyperspectralImage = sceneGet(scene,'energy');

    % Get wavelengths from scene
    wls = sceneGet(scene,'wave');
    % Get some monitor primaries
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);
elseif strcmp(img,'gray')
    % Need to call this scene to grab its primaries
    sceneImg = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(sceneImg,'fov',2);
    % Get wavelengths for creating primaries
    wls = sceneGet(scene,'wave');
    % Get some monitor primaries
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

    % Create gray hyperspectral image
    % 256 x 256 gray image
    [grayImgCalFormat,m,n] = ImageToCalFormat(ones(256,256));
    grayImgCalFormat = (0.5.*(repmat(grayImgCalFormat,3,1)));
    grayImgrgb          = CalFormatToImage(grayImgCalFormat,m,n);
    hyperspecGrayCalFormat = P_monitor * grayImgCalFormat;
    hyperspectralImage = CalFormatToImage(hyperspecGrayCalFormat,m,n);
else
    % Choose scene from manchester database of hyperspectral scenes
    dropboxDirPath = localDropboxDir();
    scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
    load(fullfile(scenesDir, img), 'scene');
    hyperspectralImage = sceneGet(scene,'energy');

    % Get wavelengths from scene
    wls = sceneGet(scene,'wave');
    % Get some monitor primaries
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);
end


% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

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
delta = plateO(size(slice_lms),deltaModulation,100);
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