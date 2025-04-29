function variance = varianceLMS(varianceType,renderType,LMS_old,LMS_new,Disp)
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

        % Add bias term???
        X = availableCones_old';  % n x 3
        y = missingCone_old';  % n x 1

        % Linear regression to predict missing cone from the available ones
        % XB = y --> B = X\y
        %
        % Deal with numerical issues
        smallestInteresingContrast = 1e-5;
        X(X < smallestInteresingContrast) = 0;
        if (all(X == 0))
            residual = y;
        else
            beta = X \ y;          % Least squares fit
            y_hat = X * beta;      % Predicted missing cone

            % Residual: part not explained by other two cones
            residual = y - y_hat;
        end

        % Weighting: use squared residuals to emphasize areas poorly explained
        weight = residual'.^2;  

        % Compute deltas in available cone planes
        delta1 = LMS_old(index(1), :) - LMS_new(index(1), :);
        delta2 = LMS_old(index(2), :) - LMS_new(index(2), :);

        % Weighted contrast loss
        variance = norm([delta1 .* weight; delta2 .* weight])^2;


    case 'detail'
        % Convert transformed LMS to dichromat simulation
        constraintWL = 585;
        LMS_CVD = DichromSimulateLinear(LMS_new, Disp.grayLMS, constraintWL, renderType, Disp);

        % Only analyze one cone channel (e.g., L) or average over channels
        % For better perceptual relevance, we can use L+M (luminance-like)
        LMS_old_lum = LMS_old(1,:) + LMS_old(2,:);    % 1 x n
        LMS_CVD_lum = LMS_CVD(1,:) + LMS_CVD(2,:);    % 1 x n

        % Convert to image dimensions (Disp.m = width, Disp.n = height)
        img_orig = reshape(LMS_old_lum, [Disp.n, Disp.m]);  % H x W
        img_cvd  = reshape(LMS_CVD_lum,  [Disp.n, Disp.m]);

        % Define local patch size (odd number)
        patchSize = 11;
        halfSize = floor(patchSize / 2);

        % Pad images to avoid border issues
        img_orig_padded = padarray(img_orig, [halfSize halfSize], 'symmetric');
        img_cvd_padded  = padarray(img_cvd,  [halfSize halfSize], 'symmetric');

        % Preallocate local variance maps
        var_orig = zeros(size(img_orig));
        var_cvd  = zeros(size(img_cvd));

        % Compute local variance in sliding window
        for i = 1:Disp.n
            for j = 1:Disp.m
                patch_orig = img_orig_padded(i:i+patchSize-1, j:j+patchSize-1);
                patch_cvd  = img_cvd_padded(i:i+patchSize-1, j:j+patchSize-1);

                mu_orig = mean(patch_orig(:));
                mu_cvd  = mean(patch_cvd(:));

                var_orig(i,j) = mean((patch_orig(:) - mu_orig).^2);
                var_cvd(i,j)  = mean((patch_cvd(:)  - mu_cvd).^2);
            end
        end

        % Compute squared difference between variances
        % varianceMap = (var_orig - var_cvd).^2;
        varianceMap = var_orig .* (var_orig - var_cvd).^2;

        % Total detail variance loss (sum or mean — your choice)
        variance = sum(varianceMap(:));  % could use mean(...) too
    case 'missingDetailGradient'
        % Convert LMS_new to dichromat-simulated LMS image
        constraintWL = 585;  % use whatever your standard is
        LMS_CVD = DichromSimulateLinear(LMS_new, Disp.grayLMS, constraintWL, renderType, Disp);

        % Compute luminance-like signals (L+M) for original and simulated images
        LMS_old_lum = LMS_old(1,:) + LMS_old(2,:);   % original luminance
        LMS_CVD_lum = LMS_CVD(1,:) + LMS_CVD(2,:);   % simulated luminance

        % Reshape to images for spatial filtering
        img_orig = reshape(LMS_old_lum, [Disp.n, Disp.m]);  % H x W
        img_cvd  = reshape(LMS_CVD_lum,  [Disp.n, Disp.m]);

        % Compute image gradients (magnitude of spatial change)
        [Gx_orig, Gy_orig] = gradient(img_orig);
        [Gx_cvd,  Gy_cvd]  = gradient(img_cvd);

        grad_orig = sqrt(Gx_orig.^2 + Gy_orig.^2);  % edge strength
        grad_cvd  = sqrt(Gx_cvd.^2  + Gy_cvd.^2);

        % Compute loss of gradient energy: what disappears after simulation
        missing_detail = max(grad_orig - grad_cvd, 0);
        weight = missing_detail .^ 2;  % square to emphasize large losses

        % Reshape to vector
        weight = reshape(weight, 1, []);

        % Compute deltas in available cone channels
        delta1 = LMS_old(index(1), :) - LMS_new(index(1), :);
        delta2 = LMS_old(index(2), :) - LMS_new(index(2), :);

        % Weighted contrast loss focused on disappearing detail
        variance = norm([delta1 .* weight; delta2 .* weight])^2;
        
    otherwise
        error('Unknown variance type specified');
end


end
