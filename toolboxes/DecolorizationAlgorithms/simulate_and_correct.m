
function simulate_and_correct(image, cbTypes)
% Simulates and applies color correction for specified color blindness types.
%
% Syntax:
%   simulate_and_correct(image, cbTypes)
%
% Inputs:
%   image: The input RGB image to process.
%   cbTypes: Array of integers specifying the types of color blindness to simulate:
%            1 = Protanopia, 2 = Deuteranopia, 3 = Tritanopia.
%
% Description:
%   This function simulates dichromatic vision for the given types of color
%   blindness using Brettel's method and then applies color correction
%   followed by re-simulation on the corrected image.
%
% Outputs:
%   Displays the simulated and corrected images.
%
% History:
%   01/10/2025 - Modularized simulation and correction process.
%
% Examples:
%   Simulate and correct for Protanopia and Deuteranopia:
%   simulate_and_correct([], [1, 2]);


    cbTypeNames = {'Protanopia', 'Deuteranopia', 'Tritanopia'};
    
    % Simulate and display initial dichromatic vision
    % figure('Name', 'Simulated Dichromatic Vision', 'NumberTitle', 'off');
    for idx = 1:length(cbTypes)
        cbType = cbTypes(idx);
        [srgb_dichromat{idx}, lmsDichromat{idx}, lmsTrichromat{idx}] = Brettel(image, cbType);
        srgb_dichromat = srgb_dichromat{idx};
        subplot(1, length(cbTypes), idx);
        imshow(srgb_dichromat{:});
        title(['Simulated ' cbTypeNames{cbType}]);
    end
    correctedLMS = optLinTransform(lmsTrichromat{2});
    correctedLMS = CalFormatToImage(correctedLMS',size(lmsTrichromat{2},1),size(lmsTrichromat{2},2));
    % Correct colors and re-simulate
    % correctedLMS = daltonization(lmsDichromat{2}, lmsTrichromat{2});
    correctedXYZ = imageLinearTransform(correctedLMS, colorTransformMatrix('lms2xyz'));
    correctedRGBImage = imageLinearTransform(correctedXYZ, colorTransformMatrix('xyz2srgb'));

    % VERY ODD THIS NORMALIZION THING IS WHAT DOES SOMETHING NEAT
    % correctedRGBImage = correctedRGBImage ./ max(correctedRGBImage); % Normalize
    correctedRGBImage = uint8(correctedRGBImage * 255); % Scale to 0-255 and cast to uint8

    correctedSimulatedImage = Brettel(correctedRGBImage, 2);
    imshow(correctedSimulatedImage{:});
end

