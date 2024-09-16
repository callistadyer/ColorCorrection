function obj = T_EstObjectiveFunction(K, D_mnew, T_mean, d, T_cones, P_monitor, m, n, bScale)

% Make sure K hasn't entered the twilight zone
if (any(isnan(K(:))))
    obj = Inf;
    return;
end

% Cone values scaled by some K, then add back in T_mean
% K is what you are optimizing
T_est = K * D_mnew + T_mean;

M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);
rgbLinImageCalFormat = M_cones2rgb*T_est;

if min(rgbLinImageCalFormat(:)) < 0 || max(rgbLinImageCalFormat(:)) > 1
    obj = Inf;

else

    % Compute rgb values from the new cone estimates
    % T_est_rgbImg = LMS2rgbLinimg(T_est, d, T_cones, P_monitor, m, n, bScale);
    [~,T_est_rgbImg] = LMS2RGBimg(T_est,d,T_cones,P_monitor,m,n,bScale);

    % Compute min and max of rgb values
    minrgb = min(T_est_rgbImg(:));
    maxrgb = max(T_est_rgbImg(:));
    if (minrgb < 0 || maxrgb > 1)
        error('Something is not consistent.');
    end

    % Goal: resulting rgb needs to be between 0 and 1
    % Calculate Pneg and Ppos
    % if minrgb >= 0
    %     Pneg = 0;
    % elseif minrgb < 0
    %     Pneg = 10e50*minrgb^2;
    % end
    %
    % if maxrgb <= 1
    %     Ppos = 0;
    % elseif maxrgb > 1
    %     Ppos = 10e50*maxrgb^2;
    % end

    if minrgb < 0 || maxrgb > 1
        obj = Inf;
        return
    end

    % Compute the penalty
    % penalty = Ppos + Pneg;

    % object to minimize
    % want to maximize k
    % largest value will be the smallest (most negative) value
    % obj = -sum(K(:).^2) + penalty;
    obj = -sum(K(:).^2);

end

end