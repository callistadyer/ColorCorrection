function [loss, info, infoNormalized, distortion, distortionNormalized, LMSCalFormat_pred, regressionBeta] = lossFunction(useLambdaOrTarget, lambdaOrTarget,...
    t_vec, LMSCalFormat, imgParams, dichromatType, infoFcn, distortionFcn, infoNormalizer, distortionNormalizer, Disp, paramsStruct)
% lossFunction  Objective for color correction optimization
%
% Syntax:
%   [loss, info, infoNormalized, distortion, distortionNormalized, LMSCalFormat_pred] = ...
%       lossFunction(useLambdaOrTarget, lambdaOrTarget, t_vec, LMSCalFormat, imgParams, dichromatType, ...
%                    infoFcn, distortionFcn, infoNormalizer, distortionNormalizer, Disp, paramsStruct)
%
% Description:
%   Applies a 3x3 transformation to the image's RGB contrast values,
%   converts the result to LMS, computes the information and distortion of
%   the transformed image, and combines them into a scalar loss for fmincon.
%
%   The loss is a lambda-weighted tradeoff:
%       loss = fminconFactor * ( -lambda*infoNormalized + (1-lambda)*distortionNormalized )
%   Info and distortion are always computed and returned, so callers that
%   only need those values (constraints, diagnostics) call it and ignore
%   the loss output.
%
% HOW THE REST OF THE CODE USES IT
%   lambda takes three kinds of values:
%       lambda = 0     ->  loss = +distortion   (minimize distortion only)
%       lambda = 1     ->  loss = -info         (maximize info only)
%       0 < lambda < 1 ->  weighted tradeoff    (lambda sweep only)
%
%   The info and distortion sweeps do NOT hit their targets through this
%   loss. A target is enforced as a separate band constraint
%   (constraintBand in colorCorrectionOptimize.m), and constraintBand also
%   calls this function, only to read off info or distortion:
%
%   Info sweep (computeSweep, sweepAxis = 'info'), via colorCorrectionOptimize('info', ...):
%       objective  = lossFunction with lambda = 0  -> minimize distortion
%       constraint = constraintBand: |info - targetInfo| <= bandTol
%       => "least distortion that reaches this info level"
%
%   Distortion sweep (computeSweep, sweepAxis = 'distortion'), via colorCorrectionOptimize('distortion', ...):
%       objective  = lossFunction with lambda = 1  -> maximize info
%       constraint = constraintBand: |distortion - targetDist| <= bandTol
%       => "most info at this distortion level"
%
%   Sweep endpoints (computeEndpoints), via colorCorrectionOptimize('lambda', ...):
%       lambda = 0 -> forced to the identity (fmincon skipped); gives info0
%       lambda = 1 -> maximize info with no constraint; gives info1
%       The sweep targets are spaced evenly between these two.
%
% Inputs:
%   useLambdaOrTarget         Must be 'lambda'. (Kept as an argument so existing
%                             calls work; older target-based modes were removed.)
%   lambdaOrTarget            Lambda, 0 to 1 (see HOW THE REST OF THE CODE USES IT)
%   t_vec                     9x1 vectorized 3x3 transformation (acts on RGB contrast)
%   LMSCalFormat              3xN LMS values of the original image (PTB CalFormat)
%   imgParams                 Struct with image parameters (m, n, ...)
%   dichromatType             'Protanopia', 'Deuteranopia', or 'Tritanopia'
%   infoFcn                   Function handle to the information metric
%   distortionFcn             Function handle to the distortion metric
%   infoNormalizer            Scalar; infoNormalized = info / infoNormalizer
%   distortionNormalizer      Scalar; distortionNormalized = distortion / distortionNormalizer
%   Disp                      Display structure
%   paramsStruct              Struct of metric parameters (e.g. predictingWhat,
%                             predictingFromWhat for computeInfo_regressCIE)
%
% Outputs:
%   loss                      Scalar loss to minimize
%   info                      Raw info value
%   infoNormalized            Normalized info value
%   distortion                Raw distortion value
%   distortionNormalized      Normalized distortion value
%   LMSCalFormat_pred         Reconstructed LMS image from the info metric, if the
%                             info function returns one (computeInfo_regressCIE
%                             does); otherwise []
%   regressionBeta            The info metric's regression weights (rows =
%                             predictors, columns = targets), if the info
%                             function returns them (computeInfo_regressCIE
%                             does); otherwise []. Used for diagnostics.

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
newLMSCalFormat = Disp.M_rgb2cones * newRGBCalFormat;

% Get the contrast LMS value
newLMSContrastCalFormat = (newLMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

%%%%%%%%%%%%%%%%%%% Calculate info in the image %%%%%%%%%%%%%%%%%%%
% To calculate info, we use contrast LMS. We compare the old and new:
LMSold = LMSContrastCalFormat;     % Original LMS contrast
LMSnew = newLMSContrastCalFormat;  % Transformed LMS contrast

% Info function. Some info functions (e.g. computeInfo_regressCIE) also
% return the reconstructed LMS image as a third output; ask for it only
% when the function provides one.
% Regression weights are only fetched when the caller asks for them.
regressionBeta = [];
if nargout(infoFcn) >= 4 && nargout >= 7
    [info, infoNormalized, LMSCalFormat_pred, regressionBeta] = infoFcn(LMSold, LMSnew, imgParams, dichromatType, infoNormalizer, Disp, paramsStruct);
elseif nargout(infoFcn) >= 3
    [info, infoNormalized, LMSCalFormat_pred] = infoFcn(LMSold, LMSnew, imgParams, dichromatType, infoNormalizer, Disp, paramsStruct);
else
    LMSCalFormat_pred = [];
    [info, infoNormalized] = infoFcn(LMSold, LMSnew, imgParams, dichromatType, infoNormalizer, Disp, paramsStruct);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Calculate distortion in the image %%%%%%%%%%%%%%%%%
% To calculate distortion, we use regular LMS excitations (not contrast).
% We compare the old and new:
LMSold = LMSCalFormat;
LMSnew = newLMSCalFormat;

% Distortion function
[distortion, distortionNormalized] = distortionFcn(LMSold, LMSnew, imgParams, distortionNormalizer, Disp, paramsStruct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale factor applied to the loss.
% Tested fminconFactor = 1 vs 1e8 (40x40 flower1 + Gaugin, L,M,S, 10/28):
% factor 1 removed most stalls (exitflag 2) but caused "no feasible
% point" failures (exitflag -2), took ~12% longer, and gave the same
% final answers. Kept at 1e8.
fminconFactor = 1e8;

% Combine info and distortion into the loss. Targets (info sweep, distortion
% sweep) are NOT handled here; they are band constraints in
% colorCorrectionOptimize (see HOW THE REST OF THE CODE USES IT above).
if ~strcmp(useLambdaOrTarget, 'lambda')
    error('lossFunction: useLambdaOrTarget must be ''lambda'' (got ''%s'').', char(string(useLambdaOrTarget)));
end

% Lambda-weighted loss function
infoWeighted       = lambdaOrTarget   * infoNormalized;
distortionWeighted = (1 - lambdaOrTarget) * distortionNormalized;

% Want info to be big. So minimize negative of info. Distortion you
% want to be small, so minimize the regular positive distortion.
loss = (fminconFactor) * (-infoWeighted + distortionWeighted);
end
