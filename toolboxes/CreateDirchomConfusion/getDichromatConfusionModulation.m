function [modulationLMS] = getDichromatConfusionModulation(rgbImageCalFormat,modulationDirection_LMS,renderType,Disp,scaleFactor,bScale)

% function that calculates the modulation in the L M or S cone based on the
% input image and the gamut limitations
%
% inputs
% rgbImageCalFormat: rgb in cal format
% renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
% takeMin:           take the minimum modulation of all the pixels?
%                    0 -> dont take min
%                    1 -> take min
% Disp includes:
%   Disp.T_cones:           spectral sensitivities
%   Disp.P_monitor:         monitor primaries
%   Disp.m:                 dimension 1 of img
%   Disp.n:                 dimension 2 of img
% scaleFactor:       how much to scale the rgbImage values (or undo it)
% bScale:            do you wanna scale or not?

 
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% THIS OVERWRITES MODULATION INPUT TO MAKE ONLY INVISIBLE SQUARES... FOR
disp('CALLISTA!! Make sure you can toggle this to modulate just one cone or random')
% Which cone do you want to modulate
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        modulationDirection_LMS = [0 1 0]';
    case 'Protanopia'   % l cone deficiency
        modulationDirection_LMS = [1 0 0]';
    case 'Tritanopia'   % s cone deficiency
        modulationDirection_LMS = [0 0 1]';
end

% Convert modulation direction from LMS space to rgb space
modulationDirection_rgb = M_cones2rgb*modulationDirection_LMS;

if round(norm(modulationDirection_LMS),4) ~= 1
    error('modulation direction not a unit vector')
end

% NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
% is different from scaleFactor that goes into rgb->LMS conversions (and
% vice versa) which scales the image values to a sensible number
% Background is the value of each pixel: this determines cone contrast separately for each pixel
for i = 1:size(rgbImageCalFormat,2)
scaleFactor_rgb(:,i) = MaximizeGamutContrast(modulationDirection_rgb,rgbImageCalFormat(:,i)); % bg is in rgb cal format
end

% Stay away from the very edge
toleranceFactor = 0.9;
scaleFactor_rgb = toleranceFactor*scaleFactor_rgb;

% Scale modulation direction by scale factor to get modulation
modulation_rgb = scaleFactor_rgb.*modulationDirection_rgb;

% If there is more than 1 value for the modulation (ie. different for each pixel) 
if size(modulation_rgb,2) > 1 
    modulation_rgb_img = CalFormatToImage(modulation_rgb,Disp.m,Disp.n);
end

% Convert modulation from rgb to LMS format for output
modulationLMS = rgbLin2LMSimg(modulation_rgb_img,Disp,scaleFactor,bScale);

end