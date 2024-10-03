function modulationLMS = getDichromatConfusionModulation(rgbImageCalFormat,renderType,takeMin,T_cones,P_monitor,scaleFactor,m,n,bScale)

% function that calculates the modulation in the L M or S cone based on the
% input image and the gamut limitations
%
% inputs
% rgbImageCalFormat: rgb in cal format
% renderType:        type of dichromacy
% takeMin:           take the minimum modulation of all the pixels?
%                    0 -> dont take min
%                    1 -> take min
% T_cones:           spectral sensitivities
% P_monitor:         monitor primaries
% scaleFactor:       how much to scale the rgbImage values (or undo it)
% m:                 dimension 1 of img
% n:                 dimension 2 of img
% bScale:            do you wanna scale or not?

% Grab 
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Which cone do you want to modulate
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        modulationDirection = [0 1 0]';
    case 'Protanopia'   % l cone deficiency
        modulationDirection = [1 0 0]';
    case 'Tritanopia'   % s cone deficiency
        modulationDirection = [0 0 1]';
end

% Convert modulation direction from LMS space to rgb space
modulationDirection_rgb = M_cones2rgb*modulationDirection;

% NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
% is different from scaleFactor that goes into rgb->LMS conversions (and
% vice versa) which scales the image values to a sensible number
% Background is the value of each pixel: this determines cone contrast separately for each pixel
for i = 1:size(rgbImageCalFormat,2)
scaleFactor_rgb(:,i) = MaximizeGamutContrast(modulationDirection_rgb,rgbImageCalFormat(:,i)); % bg is in rgb cal format
end

if takeMin == 1
% Use these two lines to use the same modulation for all pixels... will need to take the minimum mod of all 
scaleFactor_rgb = min(scaleFactor_rgb(:));        % minimum modulation of all pixels
end

% Scale modulation direction by scale factor to get modulation
modulation_rgb = scaleFactor_rgb.*modulationDirection_rgb;

% If there is more than 1 value for the modulation (ie. different for each pixel) 
if size(modulation_rgb,2) > 1 
    modulation_rgb_img = CalFormatToImage(modulation_rgb,m,n);
end

% Convert modulation from rgb to LMS format for output
modulationLMS = rgbLin2LMSimg(modulation_rgb_img,T_cones,P_monitor,scaleFactor,m,n,bScale);

end