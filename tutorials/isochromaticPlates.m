function [RGBmodulatedCalFormat lmsModuledCalFormat] = isochromaticPlates(img,renderType,lmsModulationImgFormat,bScale,options)

% function create isochromatic plates for testing dichromacy
% 
% inputs
% img:                    input image - which to create plate with?
% renderType:             type of dichromacy
% lmsModulationImgFormat: lms plate modulation in img format
% bScale:                 scale rgb values or not?

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