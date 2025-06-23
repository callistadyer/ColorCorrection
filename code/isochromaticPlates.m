function [RGBCalFormat_plate LMSCalFormat_plate] = isochromaticPlates(LMSImage,LMSImageModulation,Disp,nSquares,options)

% function create isochromatic plates for testing dichromacy
%
% Syntax:
%   [RGBCalFormat_plate LMSCalFormat_plate] = isochromaticPlates(img,renderType,LMSImageModulation,Disp,nSquares,options)
%
% Description:
%
% Inputs:
%       LMSImageModulation:      LMS plate modulation in img format
%       Disp:                    Display parameters
%       nSquares:                Number of squares for modulation
%
% Outputs:
%       RGBCalFormat_plate      Gamma corrected RGB image with cone modulation included  
%       LMSCalFormat_plate      LMS values with cone modulation added 
%
% Optional key/value pairs:
%   None
%
% See t_renderHyperspectralImage.m to see this used.


%% Pick up optional arguments
arguments
    LMSImage
    LMSImageModulation
    Disp struct
    nSquares
    options.verbose (1,1) logical = false;
end

if (options.verbose)
    fprintf('Starting execution of isochromaticPlates\n');
end

disp('Callista come back to this - make it so this func takes in LMS image')
% Load hyperspectral image
% [hyperspectralImage] = loadImage(img);

% Get LMS values
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(LMSImage);
LMSCalFormat = Disp.T_cones*hyperspectralImageCalFormat;
LMSImage     = CalFormatToImage(LMSCalFormat,Disp.m,Disp.n);

% Get original RGB image
[RGBCalFormat rgbLinCalFormat]  = LMS2RGBCalFormat(LMSCalFormat,Disp);
RGBimage                        = CalFormatToImage(RGBCalFormat,Disp.m,Disp.n);

% Create square modulations
deltaLMS = plateSquare(size(LMSImage),LMSImageModulation,nSquares);

% Add the delta to the L M S values to modulate cones (original LMS + modulation) 
lmsImage_mod = LMSImage + deltaLMS;

% CHECK IF MODULATED LMS IS IN GAMUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmsImage_modCalFormat = ImageToCalFormat(lmsImage_mod);
inGamut = checkGamut(lmsImage_modCalFormat,Disp);
if inGamut == 0
    error(['isochromaticPlates: WARNING! rgb values are out of gamut... lmsImage_mod values outside of the range [0 1]']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% turn lmsImage to cal format so you can convert to RGB
LMSCalFormat_plate = ImageToCalFormat(lmsImage_mod);

% Build matrix that goes from cones to monitor linear rgb
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

RGBCalFormat_plate = M_cones2rgb * LMSCalFormat_plate;
% convert to RGB
% RGBCalFormat_plate2 = LMS2RGBCalFormat(LMSCalFormat_plate,Disp);

% convert to image for viewing
RGBimageModulated = CalFormatToImage(RGBCalFormat_plate,m,n);

% figure('position',[927         886        1245         367]);
% subplot(1,2,1)
% imshow(RGBimage)
% title('original image')
% 
% subplot(1,2,2)
% imshow(RGBimageModulated)
% title([renderType ' testing plate'])
% 
% sgtitle('Isochromatic Plates')


end