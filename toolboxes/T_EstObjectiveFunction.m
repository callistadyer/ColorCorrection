function obj = T_EstObjectiveFunction(K, D_mnew, T_mean, d, T_cones, P_monitor, m, n, bScale)
    % Reshape K into a diagonal matrix
    % K = diag(K);
    T_est = K * D_mnew + T_mean;

    % Compute T_est_rgbImg
    T_est_rgbImg = LMS2rgbLinimg(T_est, d, T_cones, P_monitor, m, n, bScale);
    
    % Compute minrgb and maxrgb
    minrgb = min(T_est_rgbImg(:));
    maxrgb = max(T_est_rgbImg(:));
    
    % Calculate Pneg and Ppos
    if minrgb >= 0
        Pneg = 0;
    elseif minrgb < 0
        Pneg = minrgb^2;
    end

    if maxrgb <= 1
        Ppos = 0;
    elseif maxrgb > 1
        Ppos = maxrgb^2;
    end

    % Compute the penalty
    penalty = Ppos + Pneg;
    % object to minimize

    % want to maximize k
    % largest value will be the smallest (most negative) value
    obj = -sum(K(:).^2) + penalty;
end