function [insideColors, outsideColors] = chooseIshiharaColors(renderType,plateType,Disp)
% Chooses inside and outside dot colors for Ishihara-style color vision test plates
%
% Syntax:
%   [insideColors, outsideColors] = chooseIshiharaColors(renderType,plateType,Disp)
%
% Inputs:
%   renderType:     String. Type of dichromat vision to simulate or correct for. Options are:
%                       'Deuteranopia'  (M-cone missing)
%                       'Protanopia'    (L-cone missing)
%                       'Tritanopia'    (S-cone missing)
%
%   plateType:      Integer. Determines style of plate used:
%                       1   - Gray background with inside dots modulated along missing cone axis
%                       2   - Background is a set of predefined reddish/orange colors
%                       3   - All backgrounds start from gray and are modulated in different LMS directions
%                       4   - Randomly generated RGB background colors
%
%   Disp:           Struct. Display parameters. Must contain the following fields:
%                       M_cones2rgb     - Inverse of RGB-to-cone matrix
%                       grayRGB         - RGB triplet for neutral gray (for use in plateType 3)
%
% Outputs:
%   insideColors:   3×3 matrix of RGB values for "inside" dots (i.e., the number in the plate)
%   outsideColors:  3×3 matrix of RGB values for "outside" dots (i.e., the background dots)
%
%
% Examples:
%{
renderType = 'Deuteranopia';
plateType = 1;
Disp = loadDisplay('ishihara')
[insideColors, outsideColors] = chooseIshiharaColors(renderType, plateType, Disp);
%%%% Visualize the plate %%%%
ishiharaRGB = generateIshiharaPlate('74', insideColors,outsideColors,Disp.m);
figure();imagesc(ishiharaRGB)
axis square;
%}

%
% History:
%   06/19/2025  cmd  comments

if plateType == 1 % gray with missing cone mod
    outsideColors = [
        0.50    0.50    0.50;
        0.50    0.50    0.50;
        0.50    0.50    0.50;
        ];

    switch (renderType)
        case 'Deuteranopia' % m cone deficiency
            modulationDirection_LMS = [0 1 0]';
        case 'Protanopia'   % l cone deficiency
            modulationDirection_LMS = [1 0 0]';
        case 'Tritanopia'   % s cone deficiency
            modulationDirection_LMS = [0 0 1]';
    end
    
    % Direction of color (e.g., move in M cone direction to make a plate
    % that deuteranopes cannot see
    modulationDirection_rgb      = Disp.M_cones2rgb*modulationDirection_LMS;

    % NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
    % is different from scaleFactor that goes into rgb->LMS conversions (and
    % vice versa) which scales the image values to a sensible number
    % Background is the value of each pixel: this determines cone contrast separately for each pixel
    for i = 1:size(outsideColors,1)
        scaleFactor_rgb(i)      = MaximizeGamutContrast(modulationDirection_rgb,outsideColors(i,:)'); % bg is in rgb cal format
        % Stay away from the very edge
        % toleranceFactor = 0.9;
        % Scale modulation direction by scale factor to get modulation=
        modulation_rgb(i,:)      = scaleFactor_rgb(i).*modulationDirection_rgb;
    end

    % New colors. Outside colors stay the same. Inside colors simply add
    % missing cone information to existing background colors.
    % insideColorsMod  = insideColors + modulation_rgb.*[.9 .7 .5]';
    insideColors  = outsideColors + modulation_rgb;

elseif plateType == 2 % background is random colors

    % outsideColors = rand(3,3);
    % outsideColors = 0.1 + 0.8 * rand(3,3);
    outsideColors = [
        0.91, 0.53, 0.36;  % light reddish-orange
        0.95, 0.69, 0.45;  % peach
        0.86, 0.41, 0.29;  % deeper orange-red
        ];
    %     outsideColors = [
    %     0.51, 0.49, 0.54;
    %     0.47, 0.5,  0.45;
    %     0.5, 0.53, 0.46;
    % ];

    switch (renderType)
        case 'Deuteranopia' % m cone deficiency
            modulationDirection_LMS = [0 1 0]';
        case 'Protanopia'   % l cone deficiency
            modulationDirection_LMS = [1 0 0]';
        case 'Tritanopia'   % s cone deficiency
            modulationDirection_LMS = [0 0 1]';
    end

    % Direction of color (e.g., move in M cone direction to make a plate
    % that deuteranopes cannot see
    modulationDirection_rgb      = Disp.M_cones2rgb*modulationDirection_LMS;

    % NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
    % is different from scaleFactor that goes into rgb->LMS conversions (and
    % vice versa) which scales the image values to a sensible number
    % Background is the value of each pixel: this determines cone contrast separately for each pixel
    for i = 1:size(outsideColors,1)
        scaleFactor_rgb(i)      = MaximizeGamutContrast(modulationDirection_rgb,outsideColors(i,:)'); % bg is in rgb cal format
        % Stay away from the very edge
        % toleranceFactor = 0.9;
        % Scale modulation direction by scale factor to get modulation=
        modulation_rgb(i,:)      = scaleFactor_rgb(i).*modulationDirection_rgb;
    end

    % New colors. Outside colors stay the same. Inside colors simply add
    % missing cone information to existing background colors.
    % insideColorsMod  = insideColors + modulation_rgb.*[.9 .7 .5]';
    insideColors  = outsideColors + modulation_rgb;


elseif plateType == 3

    % CHECK THIS???? CONFUSED, THE DICHROMAT SHOULD BE ABLE TO SEE THIS.
    % PERHAPS YOU ARE ADDING INCORRECTLY. 
    % CHECK TO SEE IF 74 IS REVEAL WHEN YOU LEAVE 

    LS_directions = [
        .2  0  0;
        .2  0 .2;
        0  0  .2;
        ]';

    % Normalize each direction to unit length in LMS contrast space
    LS_directions = LS_directions ./ vecnorm(LS_directions);

    outsideColors = zeros(3, 3);

    for i = 1:3
        modulation_LMS = LS_directions(:,i);

        % convert to RGB
        modulation_RGB = Disp.M_cones2rgb * modulation_LMS;

        % Maximize the amount of that RGB direction we can add to gray
        scaleFactor = MaximizeGamutContrast(modulation_RGB, Disp.grayRGB);
        scaleFactor = scaleFactor/2; % dont go all the way to the edge of gamut
        % Final color: gray + scaled RGB direction
        rgbModulated = Disp.grayRGB + scaleFactor * modulation_RGB;

        outsideColors(i, :) = rgbModulated';
    end

    switch (renderType)
        case 'Deuteranopia' % m cone deficiency
            modulationDirection_LMS = [0 1 0]';
        case 'Protanopia'   % l cone deficiency
            modulationDirection_LMS = [1 0 0]';
        case 'Tritanopia'   % s cone deficiency
            modulationDirection_LMS = [0 0 1]';
    end

    % Direction of color (e.g., move in M cone direction to make a plate
    % that deuteranopes cannot see
    modulationDirection_rgb      = Disp.M_cones2rgb*modulationDirection_LMS;

    % NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
    % is different from scaleFactor that goes into rgb->LMS conversions (and
    % vice versa) which scales the image values to a sensible number
    % Background is the value of each pixel: this determines cone contrast separately for each pixel
    for i = 1:size(outsideColors,1)
        scale = [.5 .7 .9];
        scaleFactor_rgb(i)      = MaximizeGamutContrast(modulationDirection_rgb,Disp.grayRGB); % bg is in rgb cal format
        % Scale modulation direction by scale factor to get modulation=
        modulation_rgb(i,:)      = scaleFactor_rgb(i).*modulationDirection_rgb .* scale(i);
    end

    % New colors. Outside colors stay the same. Inside colors simply add
    % missing cone information to existing background colors.
    % insideColorsMod  = insideColors + modulation_rgb.*[.9 .7 .5]';
    insideColors  = Disp.grayRGB + modulation_rgb;

elseif plateType == 4 % background is random colors
    % Generate 3 colors, each with RGB values between 0.4 and 0.6
    rng(1);
    outsideColors = 0.3 + (0.7 - 0.3) * rand(3, 3);

    % 1. Create 15 random colors
    % baseColors = 0.4 + (0.6 - 0.4) * rand(15, 3);
    % 
    % % 2. Create a brightness scaling factor (e.g., 0.8 for slightly dimmer)
    % brightnessScale = 0.8;
    % 
    % % 3. Create the brightness-scaled versions
    % scaledColors = baseColors * brightnessScale;
    % 
    % % 4. Concatenate them
    % outsideColors = [baseColors; scaledColors];

    switch (renderType)
        case 'Deuteranopia' % m cone deficiency
            modulationDirection_LMS = [0 1 0]';
        case 'Protanopia'   % l cone deficiency
            modulationDirection_LMS = [1 0 0]';
        case 'Tritanopia'   % s cone deficiency
            modulationDirection_LMS = [0 0 1]';
    end

    % Direction of color (e.g., move in M cone direction to make a plate
    % that deuteranopes cannot see
    modulationDirection_rgb      = Disp.M_cones2rgb*modulationDirection_LMS;

    % NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
    % is different from scaleFactor that goes into rgb->LMS conversions (and
    % vice versa) which scales the image values to a sensible number
    % Background is the value of each pixel: this determines cone contrast separately for each pixel
    for i = 1:size(outsideColors,1)
        scaleFactor_rgb(i)      = MaximizeGamutContrast(modulationDirection_rgb,outsideColors(i,:)'); % bg is in rgb cal format
        % Stay away from the very edge
        scale = [.8 .9 1];
        scale_repeated = repmat(scale, 1, 10);
        scale = scale_repeated;
        % toleranceFactor = 0.9;
        % Scale modulation direction by scale factor to get modulation=
        modulation_rgb(i,:)      = scaleFactor_rgb(i).*modulationDirection_rgb.* scale(i);
    end

    % New colors. Outside colors stay the same. Inside colors simply add
    % missing cone information to existing background colors.
    % insideColorsMod  = insideColors + modulation_rgb.*[.9 .7 .5]';
    insideColors  = outsideColors + modulation_rgb;


else
    error('undefined plateType')
end

insideColors(insideColors<0) = 0;
insideColors(insideColors>1) = 1;

outsideColors(outsideColors<0) = 0;
outsideColors(outsideColors>1) = 1;



end