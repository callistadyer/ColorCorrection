function [triLMSCalFormatOpt, triRGBCalFormat_T, info, infoNormalized, transformRGBmatrix_opt, targetInfoVals] = ...
    computeInfoSweep(obj, LMSCalFormat, dichromatType, imgParams, nSteps)
% computeInfoSweep  Sweep through target info values and optimize color correction.
%
% Syntax:
%   [triLMSCalFormatOpt, triRGBCalFormat_T, info, infoNormalized, transformRGBmatrix_opt, targetInfoVals] = ...
%       obj.computeInfoSweep(LMSCalFormat, dichromatType, imgParams, nSteps)
%
% Description:
%   This method sweeps between the info values produced by lambda=0 and lambda=1,
%   then optimizes the transformation for each intermediate targetInfo.
%
% Inputs:
%   LMSCalFormat   - 3 x N LMS values of the input image
%   dichromatType  - type of dichromacy
%                            'Protaniopia'
%                            'Deuteranopia'
%                            'Tritanopia'
%   imgParams      - Struct with image-related params
%   nSteps         - Number of target info steps to interpolate (default: 10).
%
% Outputs:
%   triLMSCalFormatOpt      - Cell array of optimized LMS images
%   triRGBCalFormat_T       - Cell array of optimized RGB images
%   info                    - Cell array of info values at each step
%   infoNormalized          - Cell array of normalized info values
%   transformRGBmatrix_opt  - Cell array of 3x3 transformation matrices
%   targetInfoVals          - Vector of target info values used

if nargin < 5
    nSteps = 10;
end

% Step 1: Prepare simulation
Disp = obj.Disp;
LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

[calFormatLMS_prot, ~, ~] = DichromSimulateLinear(LMSCalFormat, 'Protanopia', Disp);
[calFormatLMS_deut, ~, ~] = DichromSimulateLinear(LMSCalFormat, 'Deuteranopia', Disp);
[calFormatLMS_trit, ~, ~] = DichromSimulateLinear(LMSCalFormat, 'Tritanopia', Disp);
LMSCalFormat_new = [calFormatLMS_prot(1,:); calFormatLMS_deut(2,:); calFormatLMS_trit(3,:)];
LMSContrastCalFormat_new = (LMSCalFormat_new - Disp.grayLMS) ./ Disp.grayLMS;

% Step 2: Normalize distortion and info
normalizerValueToGetRawValue = 1;
imgParams.infoNorm       = normalizerValueToGetRawValue;
imgParams.distortionNorm = normalizerValueToGetRawValue;

infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, dichromatType, Disp, imgParams, obj.infoParams);
distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, obj.distortionParams);

imgParams.infoNorm       = infoNormalizer;
imgParams.distortionNorm = distortionNormalizer;

T_prev = eye(3);

% Step 3: Get info values for lambda = 0 and lambda = 1
[~,~,info_0,~,~] = colorCorrectionOptimize("lambda", 0, LMSCalFormat, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, T_prev, Disp, imgParams);
[~,~,info_1,~,~] = colorCorrectionOptimize("lambda", 1, LMSCalFormat, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, T_prev, Disp, imgParams);

targetInfoVals = linspace(info_0, info_1, nSteps);

% Step 4: Preallocate outputs
triLMSCalFormatOpt     = cell(1, nSteps);
triRGBCalFormat_T      = cell(1, nSteps);
info                   = cell(1, nSteps);
infoNormalized         = cell(1, nSteps);
transformRGBmatrix_opt = cell(1, nSteps);

    %% SWEEP THROUGH TARGET INFOS
    T_I = eye(3);

    for i = 1:nSteps
        thisTargetInfo = targetInfoVals(i);

        % Optimize from identity
        [LMS_TI, RGB_TI, info_TI, normInfo_TI, T_TI] = colorCorrectionOptimize( ...
            "targetInfo", thisTargetInfo, LMSCalFormat, ...
            dichromatType, obj.infoFcn, obj.distortionFcn, ...
            T_I, Disp, imgParams);

        % Optimize from T_prev
        [LMS_Tprev, RGB_Tprev, info_Tprev, normInfo_Tprev, T_Tprev] = colorCorrectionOptimize( ...
            "targetInfo", thisTargetInfo, LMSCalFormat, ...
            dichromatType, obj.infoFcn, obj.distortionFcn, ...
            T_prev, Disp, imgParams);

        % Compare losses
        loss_TI    = lossFunction("targetInfo", thisTargetInfo, T_TI(:), LMSCalFormat, dichromatType, obj.infoFcn, obj.distortionFcn, Disp, imgParams);
        loss_Tprev = lossFunction("targetInfo", thisTargetInfo, T_Tprev(:), LMSCalFormat, dichromatType, obj.infoFcn, obj.distortionFcn, Disp, imgParams);

        % Select better of the two
        if loss_TI <= loss_Tprev
            triLMSCalFormatOpt{i}     = LMS_TI;
            triRGBCalFormat_T{i}      = RGB_TI;
            info{i}                   = info_TI;
            infoNormalized{i}         = normInfo_TI;
            transformRGBmatrix_opt{i} = T_TI;
            T_prev = T_TI;  % UPDATE FOR NEXT SWEEP STEP
        else
            triLMSCalFormatOpt{i}     = LMS_Tprev;
            triRGBCalFormat_T{i}      = RGB_Tprev;
            info{i}                   = info_Tprev;
            infoNormalized{i}         = normInfo_Tprev;
            transformRGBmatrix_opt{i} = T_Tprev;
            T_prev = T_Tprev;  % UPDATE FOR NEXT SWEEP STEP
        end
    end
end
