function [loss, info, infoNormalized, distortion, distortionNormalized] = lossFunction(useLambdaOrTarget, lambdaOrTarget,... 
    t_vec, LMSCalFormat, imgParams, dichromatType, infoFcn, distortionFcn, infoNormalizer, distortionNormalizer, Disp)
% lossFunction  Objective for color correction optimization
%
% Syntax:
 % [loss, info, infoNormalized] = lossFunction(useLambdaOrTargetInfo, lambdaOrTargetInfo,... 
 %    t_vec, LMSCalFormat, imgParams, dichromatType, infoFcn, distortionFcn, infoNormalizer, distortionNormalizer, Disp)
% Description:
%   Computes the loss used in optimization by applying a 3x3 transformation
%   to RGB contrast values, then converting to LMS, computing information gain
%   and distortion, and combining them based on the loss formulation.
%
% Inputs:
%   useLambdaOrTarget        
%   lambdaOrTarget            Scalar lambda (0–1) or target info value
%   t_vec                     9×1 vectorized transformation matrix
%   LMSCalFormat              3×N LMS values in PTB CalFormat
%   imgParams                 Struct with image and normalization info
%   dichromatType             String: 'Protanopia', 'Deuteranopia', or 'Tritanopia'
%   infoFcn                   Function handle to information metric
%   distortionFcn             Function handle to distortion metric
%   infoNormalizer
%   distortionNormalizer
%   Disp                      Display structure
%
% Outputs:
%   loss                      Scalar loss to minimize
%   info                      Raw info value
%   infoNormalized            Normalized info value
%   distortion                Raw distortion value
%   distortionNormalized      Normalized distortion value

    % Reshape transformation matrix
    T = reshape(t_vec, 3, 3);

    % Original RGB cal format. Also convert to contrast image. 
    RGBCalFormat = Disp.M_cones2rgb * LMSCalFormat;
    RGBContrastCalFormat = (RGBCalFormat - Disp.grayRGB) ./ Disp.grayRGB;
    LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

    % Apply transformation to RGB contrast image
    newRGBContrastCalFormat = T * RGBContrastCalFormat;

    % Add back in the gray (AKA turn it back to regular RGB image)
    newRGBCalFormat = newRGBContrastCalFormat .* Disp.grayRGB + Disp.grayRGB;

    % Now get the transformed LMS values
    newLMSCalFormat = inv(Disp.M_cones2rgb) * newRGBCalFormat;

    % Get the contrast LMS value
    newLMSContrastCalFormat = (newLMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

    %%%%%%%%%%%%%%%%%%% Calculate info in the image %%%%%%%%%%%%%%%%%%%
    % To calculate info, we use contrast LMS. We compare the old and new:
    LMSold = LMSContrastCalFormat;     % Original LMS contrast
    LMSnew = newLMSContrastCalFormat;  % Transformed LMS contrast 

    paramsStruct = struct();
    % Info function
    [info, infoNormalized] = infoFcn(LMSold, LMSnew, imgParams, dichromatType, infoNormalizer, Disp, paramsStruct);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%% Calculate distortion in the image %%%%%%%%%%%%%%%%%
    % To calculate info, we use regular LMS excitations. We compare the old and new:
    LMSold = LMSCalFormat;
    LMSnew = newLMSCalFormat;

    % Distortion function
    [distortion, distortionNormalized] = distortionFcn(LMSold, LMSnew, imgParams, distortionNormalizer, paramsStruct);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loss functions
    fminconFactor = 1e8;

    % if strcmp(lower(useLambdaOrTarget), 'info')
    %     % Squared difference from the desired info value
    %     info_diff = infoNormalized - lambdaOrTarget;
    %     info_scalar = 1e20;
    %     loss = fminconFactor * (info_scalar * info_diff^2);
    % 
    %     % want to bake in distortion here a little bit
    %     % nonlinear constraint, take parameter, compute info 
    %     % minimize distortion, subject to the constraint that the info is
    %     % what i just found here 
    %     % nonlinear constraint function that computes the info from
    %     % whatever parameters its passed, takes abs(target-actualinfo) and
    %     % have some tolerance to deviate from that targetInfo
    %     % - abs(target-actualinfo) < tolerance
    % 
    % 
    %     % CHANGE THE useLambdaOrTargetInfo TO useLambdaOrTarget
    % elseif strcmp(lower(useLambdaOrTarget), 'distortion')
    %     distortion_diff = distortionNormalized - lambdaOrTarget;
    %     distortion_scalar = 1e20;
    %     loss = fminconFactor * (distortion_scalar * distortion_diff^2);

    if strcmp(lower(useLambdaOrTarget), 'lambda')
        % Lambda-weighted loss function
        infoWeighted       = lambdaOrTarget   * infoNormalized;
        distortionWeighted = (1 - lambdaOrTarget) * distortionNormalized;
        
        % Want info to be big. So minimize negative of info. Distortion you
        % want to be small, so minimize the regular positive distortion.
        loss = (fminconFactor) * (-infoWeighted + distortionWeighted);

    % Use this distortion one when you do an unconstrained info search:
    elseif strcmp(lower(useLambdaOrTarget), 'distortion')
        
        % The distortion sweep is having a hard time. Try getting the
        % transformation T that is able to achieve the correct distortion,
        % then start your search at that T value.
        distDiff = (distortionNormalized - lambdaOrTarget).^2;

        loss = fminconFactor * distDiff;
    else
        % error('Invalid option for useLambdaOrTargetInfo: must be ''lambda'' or ''targetInfo''');
        error('gotta use lambda?')
    end
end
