function similarity = similarityLMS(similarityType,LMS_old,LMS_new)
% Computes similarity metric between original and transformed LMS image
%
% Syntax:
%
%
% Inputs:
%   similarityType:    type of similarity metric
%                           "squared"   -> sum of squared error
%                           "luminance" -> keep chromaticity similar only
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

switch (similarityType)
    % case 'angle'
    %     similarity = (LMS_new(:)'*LMS_old(:))/(norm(LMS_new(:))*norm(LMS_old(:)));
    % case 'distance'
    %     similarity = -norm(LMS_new(:)-LMS_old(:))/norm(LMS_old(:));
    case "squared"
        % is big positive when bad
        similarity = sum((LMS_new(:)-LMS_old(:)).^2);
    case "luminance"
        L_old = LMS_old(:,1);
        M_old = LMS_old(:,2);
        S_old = LMS_old(:,3);

        L_new = LMS_new(:,1);
        M_new = LMS_new(:,2);
        S_new = LMS_new(:,3);

        sum_old = L_old + M_old + S_old;
        sum_new = L_new + M_new + S_new;

        sum_old(sum_old == 0) = eps;
        sum_new(sum_new == 0) = eps;

        % Chromaticity (l = L / (L+M+S), m = M / (L+M+S))
        l_old = L_old ./ sum_old;
        m_old = M_old ./ sum_old;

        l_new = L_new ./ sum_new;
        m_new = M_new ./ sum_new;
        % Squared chromaticity difference (preserve hue!)
        chroma_diff = (l_new - l_old).^2 + (m_new - m_old).^2;

        % --- NEW LUMINANCE DEFINITION ---
        % Weighted luminance from Stockman & Sharpe 10°
        % k_L = 0.69283932;
        % k_M = 0.34967567;
        % k_S = 0;    % (ignored)
        % 
        % lum_old = k_L * L_old + k_M * M_old;
        % lum_new = k_L * L_new + k_M * M_new;
        % 
        % % Squared luminance difference
        % % luminance difference (L + M)
        % % lum_old = L_old + M_old;
        % % lum_new = L_new + M_new;
        % lum_diff = (lum_new - lum_old).^2;
        % 
        % % weights for keeping it similar
        % w_chroma = 1.0;
        % w_luminance = 0.01;

        similarity = sum(w_chroma * chroma_diff + w_luminance * lum_diff);

    otherwise
        error('Unknown similarity type specified');
end

end
