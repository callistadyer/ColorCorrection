% t_renderHyperspectralImage
%
% Demonstrate how to read and then render a hyperspectral image.

% History
%   07/30/2024  dhb, cmd  Initial go.

%% Close all open figures
clear; close all

%% Load hyperspectral image data 
defaultImage = false;
if (defaultImage)
    % This image comes with ISETBio.
    scene = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(scene,'fov',2);
else
    % This will work if you are in the Brainard Lab and have the
    % HyperspectralSceneTutorial folder on your lab dropbox path.
    dropboxDirPath = localDropboxDir();
    scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
    load(fullfile(scenesDir, 'scene3.mat'), 'scene');
end

% Get the wavelength sampling and the actual hyperspectral image data in
% energy units.
wls = sceneGet(scene,'wave');
hyperspectralImage = sceneGet(scene,'energy');
clear scene

% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

% Get some monitor primaries.  For right now we'll
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

% Get cone responses for every pixel of the hyperspectral image
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Render lms image that a trichromat or dichromat sees
lmsImageCalFormat = T_cones*hyperspectralImageCalFormat;

renderType = 'Tritanopia';
if strcmp(renderType,'trichromat')
elseif strcmp(renderType, 'Deuteranopia')
lmsImageCalFormat(2,:) = mean([lmsImageCalFormat(1,:),lmsImageCalFormat(3,:)],"all");
elseif strcmp(renderType, 'Protanopia')
lmsImageCalFormat(1,:) = mean([lmsImageCalFormat(2,:),lmsImageCalFormat(3,:)],"all");    
elseif strcmp(renderType, 'Tritanopia')
lmsImageCalFormat(3,:) = mean([lmsImageCalFormat(1,:),lmsImageCalFormat(2,:)],"all");    
end

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Get linear RGB from LMS
rgbImageCalFormat = M_cones2rgb*lmsImageCalFormat;
rgbImage = CalFormatToImage(rgbImageCalFormat,m,n);

% For right now, normalize so that maximum value in rgb is 1
rgbImage = rgbImage/max(rgbImage(:));
rgbImage(rgbImage < 0) = 0;

% Gamma correct
gammaTable = displayGet(d,'gammatable');
iGtable = displayGet(d,'inversegamma');
RGBImage = rgb2dac(rgbImage,iGtable)/(2^displayGet(d,'dacsize')-1);

% Show the imagew
figure; imshow(RGBImage);
title([renderType ' rendering']);
