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
rgbLinImageCalFormat
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

    if minrgb < 0 || maxrgb > 1
        obj = Inf;
        return
    end

    obj = -sum(K(:).^2);

end

end