% t_renderHyperspectralImage
%
% Demonstrate how to read and then render a hyperspectral image.

% History
%   07/30/2024  dhb, cmd  Initial go.

%% Close all open figures
clear; close all

%% Load hyperspectral image data 
defaultImage = true;
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
RGBImage_trichromat = RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n);


% Testing different constant values of m-cone (or other missing cone)
alpha        = linspace(0,1,20);
mCones       = lmsImageCalFormat(2,:);
contrast     = (mCones - mean(mCones)) ./ mean(mCones);
contrast_new = alpha' .* contrast;
m_new        = (contrast_new .* mean(mCones)) + mean(mCones);

l_cone = lmsImageCalFormat(1,:);
m_cone = lmsImageCalFormat(2,:);
s_cone = lmsImageCalFormat(3,:);

renderType = 'Protanopia';
if strcmp(renderType, 'Deuteranopia')

    % % use this to test different values for missing cone
    % for i = 1:size(m_new,1)
    % lmsImageCalFormat(2,:) = m_new(i);
    % RGBImage_dichromat(:,:,:,i) = RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n);
    % end

    % lmsImageCalFormat(2,:) = mean([lmsImageCalFormat(1,:),lmsImageCalFormat(3,:)],"all");

lmsImageCalFormat(2,:)  = l_cone; % replace m cones with L cone mean
RGBImage_dichromat      = RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n);

elseif strcmp(renderType, 'Protanopia')
lmsImageCalFormat(1,:) = m_cone;   
RGBImage_dichromat      = RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n);

elseif strcmp(renderType, 'Tritanopia')
lmsImageCalFormat(3,:) = mean([lmsImageCalFormat(1,:),lmsImageCalFormat(2,:)],"all");    
end


% for i = 1:size(m_new,1) 
% Show the image
figure('position',[896         896        1152         363]); 
subplot(1,2,1);
imshow(RGBImage_trichromat);
title('Trichromat rendering');

subplot(1,2,2);
imshow(RGBImage_dichromat); hold on;
title([renderType ' rendering']);
% title([renderType ' rendering: alpha = ' num2str(round(alpha(i),3)) ' constant = ' num2str(round(m_new(i),4))]);
% end


% check by reversing RGB to LMS image
LMSimage_trichromat = double(T_cones*hyperspectralImageCalFormat);
LMSimage_dichromat  = double(RGB2LMSimg(RGBImage_dichromat,d,T_cones,P_monitor,m,n));
sameLcone = isequal(LMSimage_trichromat(1,:),LMSimage_dichromat(1,:))
sameMcone = isequal(LMSimage_trichromat(2,:),LMSimage_dichromat(2,:))
sameScone = isequal(LMSimage_trichromat(3,:),LMSimage_dichromat(3,:))

function RGBImage = RGBimg(lmsImageCalFormat,d,T_cones,P_monitor,m,n)

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

end

function lmsImageCalFormat = RGB2LMSimg(RGBImage,d,T_cones,P_monitor,m,n)

% reverse the gamma correction
gammaTable = displayGet(d,'gammatable');
rgbImg = dac2rgb(RGBImage, gammaTable);

rgbImageCalFormat = ImageToCalFormat(rgbImg);

M_rgb2cones = T_cones*P_monitor;
lmsImageCalFormat = M_rgb2cones * rgbImageCalFormat;


end


