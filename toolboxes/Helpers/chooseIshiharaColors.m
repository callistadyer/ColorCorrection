function [insideColors, outsideColors] = chooseIshiharaColors(renderType,plateType,Disp)

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
        scaleFactor = scaleFactor/3; % dont go all the way to the edge of gamut
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

    insideColors(insideColors<0) = 0;
    insideColors(insideColors>1) = 1;
elseif plateType == 4 % background is random colors
    % Generate 3 colors, each with RGB values between 0.4 and 0.6
    rng(1);
    outsideColors = 0.3 + (0.7 - 0.3) * rand(3, 3);

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
    insideColors(insideColors<0) = 0;


else
    error('undefined plateType')
end



end