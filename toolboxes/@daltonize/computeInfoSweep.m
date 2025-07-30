function [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
          LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
          transformRGBmatrixSweep, targetInfoNormalized, infoNormalized, distortionNormalized] = computeInfoSweep(obj,...
          LMSCalFormat, imgParams, dichromatType, nSteps)
% computeInfoSweep  Sweep through target info values and optimize color correction.
%
% Syntax:
%    [LMSDaltonizedCalFormatSweep, rgbLinDaltonizedCalFormatSweep,...
%           LMSDaltonizedRenderedCalFormatSweep,rgbLinDaltonizedRenderedCalFormatSweep,...
%           transformRGBmatrixSweep, infoNormalized, distortionNormalized, targetInfoNormalized] = computeInfoSweep(obj,...
%           LMSCalFormat, imgParams, dichromatType, nSteps)
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
%   LMSDaltonizedCalFormatSweep            - Cell array of optimized LMS images (this is
%                                            the size of 1xnSteps)
%   rgbLinDaltonizedCalFormatSweep         - Cell array of optimized RGB images
%   LMSDaltonizedRenderedCalFormatSweep    - Cell array of
%                                            LMSDaltonizedCalFormatSweep rendered for a dichromat 
%   rgbLinDaltonizedRenderedCalFormatSweep - Cell array of
%                                            rgbLinDaltonizedCalFormatSweep rendered for a dichromat 
%   transformRGBmatrixSweep - Cell array of 3x3 transformation matrices
%   targetInfoNormalized    - Vector of target info values used
%   infoNormalized          - Cell array of normalized info values
%   distortionNormalized    - Cell array of normalized distortion values

if isempty(nSteps)
    nSteps = 10;
end

Disp = obj.Disp;
LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

[calFormatLMS_prot, ~, ~] = DichromSimulateLinear(LMSCalFormat, 'Protanopia', Disp);
[calFormatLMS_deut, ~, ~] = DichromSimulateLinear(LMSCalFormat, 'Deuteranopia', Disp);
[calFormatLMS_trit, ~, ~] = DichromSimulateLinear(LMSCalFormat, 'Tritanopia', Disp);
LMSCalFormat_new = [calFormatLMS_prot(1,:); calFormatLMS_deut(2,:); calFormatLMS_trit(3,:)];
LMSContrastCalFormat_new = (LMSCalFormat_new - Disp.grayLMS) ./ Disp.grayLMS;

% Normalize distortion and info
normalizerValueToGetRawValue = 1;

infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, normalizerValueToGetRawValue, obj.distortionParams);

% Get info values for lambda = 0 and lambda = 1
[~,~,~,info_0,infoNormalized_0,distortion0, distortionNormalized0] = colorCorrectionOptimize("lambda", 0, LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

[~,~,~,info_1,infoNormalized_1,distortion1, distortionNormalized1] = colorCorrectionOptimize("lambda", 1, LMSCalFormat, imgParams, dichromatType, ...
    obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

% Get target info values interpolated between lambdas 0 and 1 
targetInfoNormalized = linspace(infoNormalized_0, infoNormalized_1, nSteps);

% Preallocate outputs
LMSDaltonizedCalFormatSweep     = cell(1, nSteps);
rgbLinDaltonizedCalFormatSweep  = cell(1, nSteps);
LMSDaltonizedRenderedCalFormatSweep      = cell(1, nSteps);
rgbLinDaltonizedRenderedCalFormatSweep   = cell(1, nSteps);
infoNormalized           = cell(1, nSteps);
distortionNormalized     = cell(1, nSteps);
transformRGBmatrixSweep  = cell(1, nSteps);

    %% Sweep through target infos 
    T_I = eye(3);
    T_prev = eye(3);
    for i = 1:nSteps
        thisTargetInfo = targetInfoNormalized(i);

        % Optimize from identity starting point
        [LMS_TI, RGB_TI, T_TI, info_TI,normInfo_TI,distortion_TI, normDistortion_TI] = colorCorrectionOptimize( ...
            "targetInfo", thisTargetInfo, LMSCalFormat, imgParams, ...
            dichromatType, obj.infoFcn, obj.distortionFcn, ...
            infoNormalizer, distortionNormalizer,...
            Disp,'T_init',T_I);

        % Optimize from T_prev starting point
        [LMS_Tprev, RGB_Tprev, T_Tprev, info_Tprev, normInfo_Tprev,distortion_Tprev, normDistortion_Tprev] = colorCorrectionOptimize( ...
            "targetInfo", thisTargetInfo, LMSCalFormat, imgParams, ...
            dichromatType, obj.infoFcn, obj.distortionFcn, ...
            infoNormalizer, distortionNormalizer,...
            Disp,'T_init',T_prev);

        % Compare losses
        loss_TI    = lossFunction("targetInfo", thisTargetInfo, T_TI(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);
        loss_Tprev = lossFunction("targetInfo", thisTargetInfo, T_Tprev(:), LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

        % Select better of the two
        if loss_TI <= loss_Tprev
            LMSDaltonizedCalFormatSweep{i}     = LMS_TI;
            rgbLinDaltonizedCalFormatSweep{i}  = RGB_TI;
            infoNormalized{i}          = normInfo_TI;
            distortionNormalized{i}    = normDistortion_TI;
            transformRGBmatrixSweep{i} = T_TI;
            T_prev = T_TI;  % Update T_prev based on previous step 
        else
            LMSDaltonizedCalFormatSweep{i}     = LMS_Tprev;
            rgbLinDaltonizedCalFormatSweep{i}  = RGB_Tprev;
            infoNormalized{i}          = normInfo_Tprev;
            distortionNormalized{i}    = normDistortion_Tprev;
            transformRGBmatrixSweep{i} = T_Tprev;
            T_prev = T_Tprev;  % Update T_prev based on previous step 
        end
         
        [LMSDaltonizedRenderedCalFormatSweep{i},rgbLinDaltonizedRenderedCalFormatSweep{i},~] = DichromSimulateLinear(LMSDaltonizedCalFormatSweep{i},dichromatType,Disp);
    end
end
