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
%                           "delta"      -> weight remaining cone contrast
%                            by where missing cone contrast is greatest
%                           "LMdifferenceContrast" -> weight is now
%                            determined by the difference in L and M cone
%                            activation instead of just missing cone
%                            contrast
%   renderType:        type of dichromacy
%   LMS_old:           original LMS values
%   LMS_new:           transformed LMS values
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
% Variance term
switch (varianceType)
    case 'newConeVar'
        variance = (var(LMS_new(index(1), :)) + var(LMS_new(index(2), :)));
    case 'delta'
        delta1 = LMS_old(index(1),:) - LMS_new(index(1),:);
        delta2 = LMS_old(index(2),:) - LMS_new(index(2),:);

        missingCone = LMS_old(missingIdx,:); % old m
        variance = norm([delta1 .* missingCone; delta2 .* missingCone])^2;

        % weight = abs(missingCone).^alpha;
        % weight = missingCone;

        % Relative to background
        % M_background = mean(missingCone);
        % relative_M = missingCone - M_background;
        % alpha = 2;
        % weight = max(relative_M, 0).^alpha;

    case 'LMdifferenceContrast'

        % NOTE: needs to be adapted for all color blindness types...

        % Extract L and M contrast images
        L = LMS_old(1, :);
        M = LMS_old(2, :);

        % Compute absolute difference
        LM_diff = abs(L - M);

        % Apply power to emphasize large differences
        alpha = 2;
        weight = LM_diff .^ alpha;
        weight = reshape(weight, 1, []);  % Ensure row vector

        % Compute deltas in available cone planes
        delta1 = LMS_old(index(1), :) - LMS_new(index(1), :);
        delta2 = LMS_old(index(2), :) - LMS_new(index(2), :);

        % Weighted contrast loss
        variance = norm([delta1 .* weight; delta2 .* weight])^2;

    case 'regress'
        % Extract relevant cone channels
        availableCones_old = LMS_old(index, :);        % 2 x n
        missingCone_old    = LMS_old(missingIdx, :);   % 1 x n

        % Add bias term for regression 
        X = [availableCones_old', ones(size(availableCones_old, 2), 1)];  % n x 3
        y = missingCone_old';  % n x 1

        % Linear regression to predict missing cone from the available ones
        beta = X \ y;          % Least squares fit
        y_hat = X * beta;      % Predicted missing cone

        % Residual: part not explained by other two cones
        residual = y - y_hat;  % n x 1

        % Weighting: use squared residuals to emphasize areas poorly explained
        weight = residual'.^2;  % 1 x n

        % Compute deltas in available cone planes
        delta1 = LMS_old(index(1), :) - LMS_new(index(1), :);
        delta2 = LMS_old(index(2), :) - LMS_new(index(2), :);

        % Weighted contrast loss
        variance = norm([delta1 .* weight; delta2 .* weight])^2;



    otherwise
        error('Unknown variance type specified');
end


end
