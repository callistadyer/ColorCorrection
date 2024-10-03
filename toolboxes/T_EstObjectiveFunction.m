function obj = T_EstObjectiveFunction(kVec, D_mnew, T_mean, d, T_cones, P_monitor, m, n, bScale)

% Make sure K hasn't entered the twilight zone
if (any(isnan(kVec(:))))
    obj = Inf;
    return;
end

% Convert kVec to diagonal matrix K
K = diag(kVec);

% Cone values scaled by some K, then add back in T_mean
% K is what you are optimizing
T_est = K * D_mnew;
% T_est = K * D_mnew + T_mean;

% Get rgb values
M_rgb2cones = T_cones*P_monitor;
M_cones2rgb = inv(M_rgb2cones);
rgbLinImageCalFormat = M_cones2rgb*T_est;

% Ensure rgb is not <0 or >1
if min(rgbLinImageCalFormat(:)) < 0 || max(rgbLinImageCalFormat(:)) > 1
    obj = Inf;

else

    % Compute rgb values from the new cone estimates
    % T_est_rgbImg = LMS2rgbLinCalFormat(T_est, d, T_cones, P_monitor, m, n, bScale);
    [~,T_est_rgbImg] = LMS2RGBCalFormat(T_est,d,T_cones,P_monitor,m,n,bScale);

    % Compute min and max of rgb values
    minrgb = min(T_est_rgbImg(:));
    maxrgb = max(T_est_rgbImg(:));
    if (minrgb < 0 || maxrgb > 1)
        error('Something is not consistent.');
    end

    if minrgb < 0 || maxrgb > 1
        obj = Inf;
        return
    end

    obj = -sqrt(sum(K(:).^2));

end

end