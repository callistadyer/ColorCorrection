function [hyperspectralImage wls d P_monitor] = loadImage(image)

% Loads in images and makes them hyperspectral images for ColorCorrection
% project (dichromat rendering)
%
% Syntax:
%   [hyperspectralImage wls d P_monitor] = loadImage(image)
%
% Description:
%
% Inputs:
%   lmsImageCalFormat     - String. Image to be loaded. Options include:
%                               'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                               'gray'       - Gray spatially uniform field.  
% Outputs:
%   hyperspectralImage    - hyperspectral image of the input image
%   wls                   - wavelengths
%   d                     - display info
%   P_monitor             - monitor primaries
%
% Optional key/value pairs:
%   None
%
%% Load hyperspectral image data 
if isempty(image)
    % This image comes with ISETBio.
    scene = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
    scene = sceneSet(scene,'fov',2);
    wls = sceneGet(scene,'wave');
    hyperspectralImage = sceneGet(scene,'energy');

    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);
elseif strcmp(image,'gray')
    % Get some monitor primaries
    wls = (400:10:700)';
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);

    % Create gray hyperspectral image
    %
    % 256 x 256 gray image
    % [grayImgCalFormat,m,n] = ImageToCalFormat(ones(256,256));
    [grayImgCalFormat,m,n] = ImageToCalFormat(ones(128,128));

    % Gray 0.5 rgb at each pixel in image 
    grayImgCalFormat       = (0.5.*(repmat(grayImgCalFormat,3,1)));
    grayImgrgb             = CalFormatToImage(grayImgCalFormat,m,n);

    % Make hyperspectral img by multiplying primaries * rgb values at each pixel 
    hyperspecGrayCalFormat = P_monitor * grayImgCalFormat;

    % Image format
    hyperspectralImage     = CalFormatToImage(hyperspecGrayCalFormat,m,n);
else
    % This will work if you are in the Brainard Lab and have the
    % HyperspectralSceneTutorial folder on your lab dropbox path.
    % Get some monitor primaries
    dropboxDirPath = localDropboxDir();
    scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
    load(fullfile(scenesDir, image), 'scene');
    hyperspectralImage = sceneGet(scene,'energy');
    wls = sceneGet(scene,'wave');
    
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);
end

end