function [modulationLMSimage] = getDichromatConfusionModulation(rgbImageCalFormat,modulationDirection_LMS,modType,Disp,scaleFactor)

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

 
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Which cone do you want to modulate
switch (modType)
    case 'rand'
        % This code case allows the squares to be random colors
        modulationDirection_LMS = modulationDirection_LMS;
    case 'M' % m cone deficiency
        modulationDirection_LMS = [0 1 0]';
    case 'L'   % l cone deficiency
        modulationDirection_LMS = [1 0 0]';
    case 'S'   % s cone deficiency
        modulationDirection_LMS = [0 0 1]';
end

% Convert modulation direction from LMS space to rgb space
modulationDirection_rgb = M_cones2rgb*modulationDirection_LMS;


% NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
% is different from scaleFactor that goes into rgb->LMS conversions (and
% vice versa) which scales the image values to a sensible number
% Background is the value of each pixel: this determines cone contrast separately for each pixel
for i = 1:size(rgbImageCalFormat,2)
scaleFactor_rgb(:,i) = MaximizeGamutContrast(modulationDirection_rgb,rgbImageCalFormat(:,i)); % bg is in rgb cal format
end

% Stay away from the very edge
toleranceFactor = 0.8;
scaleFactor_rgb = toleranceFactor*scaleFactor_rgb;

% Scale modulation direction by scale factor to get modulation
modulation_rgb = scaleFactor_rgb.*modulationDirection_rgb;

% If there is more than 1 value for the modulation (ie. different for each pixel) 
% if size(modulation_rgb,2) > 1 
%     modulation_rgb_img = CalFormatToImage(modulation_rgb,Disp.m,Disp.n);
% end

% Convert modulation from rgb to LMS format for output
modulationLMSCalFormat = rgbLin2LMSCalFormat(modulation_rgb,Disp,scaleFactor);
modulationLMSimage = CalFormatToImage(modulationLMSCalFormat,Disp.m,Disp.n);
end