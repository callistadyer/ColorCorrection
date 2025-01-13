function correctedLMS = daltonization(lmsDichromat, lmsTrichromat)
% Applies color correction for dichromatic viewers in LMS space.
%
% Syntax:
%   correctedImage = color_correction(image, lmsDichromat, lmsTrichromat)
%
% Inputs:
%   image: The input RGB image.
%   lmsDichromat: LMS data of the dichromatic simulation.
%   lmsTrichromat: LMS data of the original trichromatic image.
%
% Description:
%   This function calculates the LMS error between dichromatic and
%   trichromatic vision, applies a weighted correction, and normalizes
%   the LMS values before converting them back to RGB format.
%
% Outputs:
%   correctedImage: The color-corrected RGB image.
%
% History:
%   01/10/2025 - Added normalization for corrected LMS values.
%
% Examples:
%   correctedImage = color_correction(lmsDichromat, lmsTrichromat);

    errorLMS = lmsTrichromat - lmsDichromat{:};
    weight = 2.5; % Weight for correction
    correctedLMS = lmsTrichromat + weight .* errorLMS;
    % correctedLMS = correctedLMS ./ max(correctedLMS(:)); % Normalize

end
