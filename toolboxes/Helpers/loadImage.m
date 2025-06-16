function [hyperspectralImage Disp] = loadImage(image)
error('WARNING: THIS IS OUT OF DATE... WE DONT GET THIS INFO FROM THIS FUNCTION ANYMORE')
end

% % Loads in images and makes them hyperspectral images for ColorCorrection
% % project (dichromat rendering)
% %
% % Syntax:
% %   [hyperspectralImage wls d P_monitor] = loadImage(image)
% %
% % Description:
% %
% % Inputs:
% %   lmsImageCalFormat     - String. Image to be loaded. Options include:
% %                               'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
% %                               'gray'       - Gray spatially uniform field.  
% % Outputs:
% %   hyperspectralImage    - hyperspectral image of the input image
% %   wls                   - wavelengths
% %   d                     - display info
% %   P_monitor             - monitor primaries
% %
% % Optional key/value pairs:
% %   None
% %
% %% Load hyperspectral image data 
% Disp = struct();
% 
% if isempty(image)
%     % This image comes with ISETBio.
%     scene  = sceneFromFile('StuffedAnimals_tungsten-hdrs','multispectral');
%     scene  = sceneSet(scene,'fov',2);
%     % imgXYZ = sceneGet(scene,'xyz');
%     % whiteXYZ = sceneGet(scene,'illuminant xyz');
%     wls    = sceneGet(scene,'wave');
%     hyperspectralImage = sceneGet(scene,'energy');
% 
%     d = displayCreate('LCD-Apple');
%     P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);
% elseif strcmp(image,'gray')
%     % Get some monitor primaries
%     wls = (400:10:700)';
%     % Create gray hyperspectral image
%     % 256 x 256 gray image
%     % [grayImgCalFormat,m,n] = ImageToCalFormat(ones(256,256));
%     [grayImgCalFormat,m,n] = ImageToCalFormat(ones(32,32));
% 
%     % Gray 0.5 rgb at each pixel in image 
%     grayImgCalFormat       = (0.5.*(repmat(grayImgCalFormat,3,1)));
%     grayImgrgb             = CalFormatToImage(grayImgCalFormat,m,n);
% 
%     d = displayCreate('LCD-Apple');
%     P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);
%     % Make hyperspectral img by multiplying primaries * rgb values at each pixel 
%     % This is a weighted sum of primaries 
%     hyperspecGrayCalFormat = P_monitor * grayImgCalFormat;
% 
%     % Image format
%     hyperspectralImage     = CalFormatToImage(hyperspecGrayCalFormat,m,n);
% 
%     % would be nice to create a scene so we have access to this stuff:
%     scene = sceneFromFile(grayImgrgb, 'rgb', [], d, wls);
% 
%     imgXYZ = sceneGet(scene,'xyz');
%     whiteXYZ = sceneGet(scene,'illuminant xyz');
%     Disp.m         = m;
%     Disp.n         = n;
%     Disp.imgXYZ    = imgXYZ;
%     Disp.whiteXYZ       = whiteXYZ;
% 
% elseif strcmp(image,'74')
%     % Ishihara test plate
%     rgbImage = imread('74.jpg');
%     % Resize and make sure it's appropriate data class
%     rgbImage = imresize(rgbImage, [128, 128]);
%     rgbImage = double(rgbImage)/255;
%     % Some numbers were slightly out of gamut... fix that
%     rgbImage(rgbImage>1) = .99;
% 
%     % Display
%     % Wavelengths for display
%     wls = (400:10:700)';
%     d = displayCreate('LCD-Apple');
%     P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);
% 
%     % Make scene from rgb values
%     scene = sceneFromFile(rgbImage, 'rgb', [], d, wls);
%     % Get XYZ values (needed for dichromat simulation algorithm from Brettel)  
%     imgXYZ = sceneGet(scene,'xyz');
%     % Hyperspectral image (needed for calculating cone responses)
%     hyperspectralImage = double(sceneGet(scene,'energy'));
%     % Extra parameters
%     Disp.m         = 128;
%     Disp.n         = 128;
%     Disp.imgXYZ    = imgXYZ;
% else
%     % This will work if you are in the Brainard Lab and have the
%     % HyperspectralSceneTutorial folder on your lab dropbox path.
%     % Get some monitor primaries
%     dropboxDirPath = localDropboxDir();
%     scenesDir = fullfile(dropboxDirPath, 'HyperspectralSceneTutorial', 'resources', 'manchester_database', '2004');
%     load(fullfile(scenesDir, image), 'scene');
%     hyperspectralImage = sceneGet(scene,'energy');
%     wls = sceneGet(scene,'wave');
% 
%     d = displayCreate('LCD-Apple');
%     P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),wls);
% end
% 
% % Get cone spectral sensitivities
% load T_cones_ss2;
% T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);
% % Converting into XYZ space for Brettel code
% load T_xyz1931.mat
% T_xyz = SplineCmf(S_xyz1931,T_xyz1931,wls);
% 
% % Save display parameters for easy calling later
% Disp.T_cones   = T_cones;
% Disp.d         = d;
% Disp.P_monitor = P_monitor;
% Disp.wls       = wls;
% Disp.T_xyz = T_xyz;
% 
% 
% end