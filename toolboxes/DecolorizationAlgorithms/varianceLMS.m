function variance = varianceLMS(varianceType,renderType,LMS_old,LMS_new)
% Computes variance metric between original and transformed LMS image
%
% Syntax:
%
%
% Inputs:
%   varianceType:      type of variance metric
%                           "newConeVar" -> maximize variance on the
%                            available cones of the transformed LMS image
%   renderType:        type of dichromacy
%   LMS_old:           original LMS values
%   LMS_new:           transformed LMS values
%   LMSc_new           transformed LMS values, contrast
%
% Outputs:
%   triLMSCalFormatOpt: Transformed LMS values
%
% Constraints:
%   - RGB values must be between 0 and 1
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        index = [1 3];
        missingIdx = 2;
    case 'Protanopia'   % l cone deficiency
        index = [2 3];
        missingIdx = 1;
    case 'Tritanopia'   % s cone deficiency
        index = [1 2];
        missingIdx = 3;
end
varianceType = 'LMdifferenceContrast';
% Variance term
switch (varianceType)
    case 'newConeVar'
        variance = (var(LMS_new(index(1), :)) + var(LMS_new(index(2), :)));
    case 'delta'
        delta1 = LMS_old(index(1),:) - LMS_new(index(1),:);
        delta2 = LMS_old(index(2),:) - LMS_new(index(2),:);

        missingCone = LMS_old(missingIdx,:); % old m

        alpha = 4;  %
        % weight = abs(missingCone).^alpha;
        % weight = missingCone;

        % Relative to background
        M_background = mean(missingCone);
        relative_M = missingCone - M_background;
        alpha = 2;
        weight = max(relative_M, 0).^alpha;
        weight = weight / max(weight);

        % variance = norm([delta1 .* missingCone; delta2 .* missingCone])^2;
        variance = norm([delta1 .* weight; delta2 .* weight])^2;
    case 'DoGweightedContrast'
        % Reshape missing cone contrast to 2D image
        sz = sqrt(size(LMS_old, 2));
        if mod(sz,1) ~= 0
            error('Image must be square or provide known reshaping strategy.');
        end
        missingConeImg = reshape(LMS_old(missingIdx, :), sz, sz);

        % Difference of Gaussians to find local M contrast bumps
        sigma1 = 0.5; sigma2 = 2;
        fine   = imgaussfilt(missingConeImg, sigma1);
        coarse = imgaussfilt(missingConeImg, sigma2);
        dogMap = fine - coarse;

        % Only focus on extra positive M contrast
        dogMap = max(dogMap, 0);

        weight = dogMap(:).^2;  % you can tune power

        % Compute deltas
        delta1 = LMS_old(index(1), :) - LMS_new(index(1), :);
        delta2 = LMS_old(index(2), :) - LMS_new(index(2), :);

        % Final weighted contrast loss
        variance = norm([delta1 .* weight'; delta2 .* weight'])^2;
    case 'LMdifferenceContrastDOG'
        % Reshape LMS contrast to 2D images
        sz = sqrt(size(LMS_old, 2));
        if mod(sz,1) ~= 0
            error('Image must be square or provide known reshaping strategy.');
        end

        L_img = reshape(LMS_old(1, :), sz, sz);
        M_img = reshape(LMS_old(2, :), sz, sz);

        % Compute L-M opponency map (absolute difference)
        LM_diff = abs(L_img - M_img);

        % Optionally emphasize spatial structure with DoG
        sigma1 = 0.5; sigma2 = 2;
        fine   = imgaussfilt(LM_diff, sigma1);
        coarse = imgaussfilt(LM_diff, sigma2);
        dogMap = fine - coarse;

        % Positive contrast emphasis
        dogMap = max(dogMap, 0);
        weight = dogMap(:).^2;

        % Compute deltas for available cones
        delta1 = LMS_old(index(1), :) - LMS_new(index(1), :);
        delta2 = LMS_old(index(2), :) - LMS_new(index(2), :);

        % Final weighted contrast loss
        % weight = reshape(weight, 1, []);  % Ensure it's a row
        variance = norm([delta1 .* weight'; delta2 .* weight'])^2;
    case 'LMdifferenceContrast'

        % Extract L and M contrast images
        L = LMS_old(1, :);
        M = LMS_old(2, :);

        % Compute absolute difference (opponency cue)
        LM_diff = abs(L - M);

        % Apply power to emphasize large differences
        alpha = 2;
        weight = LM_diff .^ alpha;
        weight = reshape(weight, 1, []);  % Ensure row vector

        % Compute deltas in available cone planes
        delta1 = LMS_old(index(1), :) - LMS_new(index(1), :);
        delta2 = LMS_old(index(2), :) - LMS_new(index(2), :);

        % Final weighted contrast loss
        variance = norm([delta1 .* weight; delta2 .* weight])^2;
    otherwise
        error('Unknown variance type specified');
end


end
