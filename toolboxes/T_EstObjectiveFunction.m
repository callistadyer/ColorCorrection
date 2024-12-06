function obj = T_EstObjectiveFunction(kVec, D_mnew, T_mean, Disp, bScale)
% Objective function for choosing how to scale values after PCA to be in
% similar range as LMS values
%
% Syntax:
%   obj = T_EstObjectiveFunction(kVec, D_mnew, T_mean, Disp, bScale)
%
% Description:
%
% Inputs:
%   kVec:         - [1x3] vector. Initial guess for scalars on PCA outputs                    
%   D_mnew:       - mean subtracted LMS values mapped onto two principle
%                   components (and one is copied into the row of the
%                   missing cone)
%   T_mean        - mean values of original LMS
%   d             - Struct.  Contains display information, displayCreate('LCD-Apple'); 
%   T_cones       - [3xnWl]. Cone spectral sensitivities
%   P_monitor     - [nWlx3]. Display primaries
%   m             - Scalar.  Row dimension of image     
%   n             - Scalar.  Column dimension of image    
%   bScale:       - Boolean. Scale the image values into display range (1
%                   or 0).  A good idea except for 'gray'.
%
% Outputs:
%   obj:          - objective to minimize in fmincon
%                   see colorCorrectionPCA.m
%
% Optional key/value pairs:
%   None
%


% Make sure K hasn't entered the twilight zone
if (any(isnan(kVec(:))))
    obj = Inf;
    return;
end

% Convert kVec to diagonal matrix K
K = diag(kVec);

% Cone values scaled by some K, then add back in T_mean
% K is what you are optimizing
% T_est = K * D_mnew;
T_est = K * D_mnew + T_mean;

% Get rgb values
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);
rgbLinImageCalFormat = M_cones2rgb*T_est;

% Ensure rgb is not <0 or >1
if min(rgbLinImageCalFormat(:)) < 0 || max(rgbLinImageCalFormat(:)) > 1
    obj = Inf;

else

    % Compute rgb values from the new cone estimates
    % T_est_rgbImg = LMS2rgbLinCalFormat(T_est, d, T_cones, P_monitor, m, n, bScale);
    [~,T_est_rgbImg] = LMS2RGBCalFormat(T_est,Disp,bScale);

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