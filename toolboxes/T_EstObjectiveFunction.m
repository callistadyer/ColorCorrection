function obj = T_EstObjectiveFunction(kVec, triLMSCalFormatOpt, triLMSmeans, renderType, Disp, bScale)
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
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        cbType = 2;
    case 'Protanopia'   % l cone deficiency
        cbType = 1;
    case 'Tritanopia'   % s cone deficiency
        cbType = 3;
end

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
triLMSCalFormatEst = K * triLMSCalFormatOpt + triLMSmeans;

% Get rgb values
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% dichromat rendering
% LMS --> Linear RGB (so we can go from RGB --> XYZ)
triRGBlinCalFormat = LMS2rgbLinCalFormat(triLMSCalFormatEst,Disp,0);

% Matrix to convert from rgb to xyz
% Note: this matrix must be applied on the LEFT!!
M_rgb2xyz = Disp.T_xyz*Disp.P_monitor;
M_xyz2rgb = inv(M_rgb2xyz);

% Linear RGB --> XYZ (Brettel takes in XYZ)
triXYZCalFormat = M_rgb2xyz * triRGBlinCalFormat;
% Cal Format --> Image Format
triXYZImgFormat = CalFormatToImage(triXYZCalFormat,Disp.m,Disp.n);

%%%%%%%%% Dichromat Simulation (Brettel) %%%%%%%%%%%%%%%%%%%%%%
[diXYZ] = DichromatSimulateBrettel(triXYZImgFormat, cbType, []);

% Image Format --> CalFormat
diXYZCalFormat = ImageToCalFormat(diXYZ);
% XYZ --> Linear RGB
diRGBLinCalFormat = M_xyz2rgb * diXYZCalFormat;
diRGBLinImgFormat = CalFormatToImage(diRGBLinCalFormat,Disp.m,Disp.n);
% Quick snipping if the vals are only over by a small amount
if (max(diRGBLinImgFormat(:)) > 1) && (max(diRGBLinImgFormat(:)) < 1.01)
    diRGBLinImgFormat(diRGBLinImgFormat>1) = .99;
end
% Linear RGB --> LMS 
diLMSImgFormat = rgbLin2LMSimg(diRGBLinImgFormat,Disp,1,0);
% Image --> Cal Format
diLMSCalFormat = ImageToCalFormat(diLMSImgFormat);



% DichromatSimulateBrettel()
% dichromatLMSCalFormat = DichromatSimulateBrettel(triLMSCalFormatEst,renderType,Disp);
% % Get rgb values
% M_rgb2cones = Disp.T_cones*Disp.P_monitor;
% M_cones2rgb = inv(M_rgb2cones);
% di_rgbLinImageCalFormat = M_cones2rgb*dichromatLMSCalFormat;

% Ensure rgb is not <0 or >1
if min(triRGBlinCalFormat(:)) < 0 || max(triRGBlinCalFormat(:)) > 1 ...
        || min(diRGBLinImgFormat(:)) < 0 || max(diRGBLinImgFormat(:)) > 1 
    obj = Inf;

else

    % Compute rgb values from the new cone estimates
    % Trichromat rendering
    [~,T_est_rgbImg]          = LMS2RGBCalFormat(triLMSCalFormatEst,Disp,bScale);
    % Dichromat rendering
    [~,dichromatRGBCalFormat] = LMS2RGBCalFormat(diLMSCalFormat,Disp,bScale);

    allrgb = [T_est_rgbImg(:), dichromatRGBCalFormat(:)];
    % Compute min and max of rgb values
    % minrgb = min(T_est_rgbImg(:));
    % maxrgb = max(T_est_rgbImg(:));
    minrgb = min(allrgb(:));
    maxrgb = max(allrgb(:));
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