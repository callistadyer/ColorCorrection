function [LMSDaltonizedCalFormat, rgbLinDaltonizedCalFormat, transformRGBmatrix, info, infoNormalized, distortion, distortionNormalized] = ...
    colorCorrectionOptimize(target, lambda, targetInfo, targetDist, ...
    triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, varargin)
% Optimizes a linear transformation to enhance color contrast for dichromats
%
% % Syntax:
% [triLMSCalFormatOpt, trirgbLinCalFormatOpt, transformRGBmatrix, info, infoNormalized, distortion, distortionNormalized] = ...
%     colorCorrectionOptimize(useLambdaOrTargetInfo, lambdaOrTargetInfo, ...
%     triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, varargin)
% Inputs:
%   useLambdaOrTargetInfo: String. Optimization mode:
%                           'lambda'       – vary tradeoff weight lambda (0 to 1)
%                           'targetInfo'   – vary based on target contrast info  
%
%   lambdaOrTargetInfo:   Either the lambda or the var value, depending on useLambdaOrTargetInfo 
%                                 If 'lambdaOrTargetInfo' = 'lambda': (0 < lambda < 1). 
%                                    Tradeoff weight balancing contrast maximization
%                                    and similarity to the original image.
%                                 If 'lambdaOrTargetInfo' = 'targetinfo', this specifies
%                                    a contrast info target
%
%   triLMSCalFormat:    Nx3 matrix of original LMS-calibrated image data to be transformed.
%
%   imgParams:          image parameters
%
%   dichromatType:      String. Type of color vision deficiency:
%                           'Deuteranopia' (M-cone missing)
%                           'Protanopia'   (L-cone missing)
%                           'Tritanopia'   (S-cone missing)
%
%   infoFnc:            Function handle used to define how much info is in the image.
%
%   distortionFcn:      Function handle used to define similarity metric to preserve naturalness.
%                           "squared"   -> sum of squared error
%                           "luminance" -> keep chromaticity similar only
%   infoNormalizer:     normalizing constant
%   distortionNormalizer:  normalizing constant
%
%   Disp:               Struct with display parameters. Must include:
%                           .M_cones2rgb  – Matrix to convert LMS to RGB
%                           .T_cones      – Cone fundamentals
%                           .P_monitor    – Monitor spectral power distribution
%                           .wls          – Wavelength sampling
%
% Outputs:
%   LMSDaltonizedCalFormat:    daltonized LMS rendering for trichromat
%   rgbLinDaltonizedCalFormat: daltonized rgb Linear rendering for trichromat
%   transformRGBmatrix:        transformation matrix
%   info:                      info from infoFcn
%   infoNormalized:            info from infoFcn, normalized
%   distortion:                distortion from distortionFcn
%   distortionNormalized:      distortion from distortionFcn, normalized 
%
% Constraints:
%   - RGB values must be between 0 and 1
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

% History
%   09/14/2024  cmd  Initial go.
%
% Examples:
%{
%}

% Sample data to see how code workswe
if isempty(triLMSCalFormat)
    % Pretend data to visualize what the function does
    N = 100;                         % Number of points
    x = randn(1, N);                 % Random data for the first dimension
    y = 2 * x + randn(1, N) * 0.2;   % Strong correlation with the first dimension
    z = 0.5 * x + 0.3 * randn(1, N); % Weaker correlation with the first dimension
    triLMSCalFormat = [x; y; z];
end

%% Key value pairs
parser = inputParser;
parser.addParameter('T_init', eye(3), @(x) isnumeric(x) && isequal(size(x),[3 3]));
parser.parse(varargin{:});
T_init  = parser.Results.T_init;


switch lower(target)
    case 'lambda'
        if isempty(lambda), error('Mode ''lambda'' requires a scalar lambda.'); end
    case 'info'
        if isempty(targetInfo), error('Mode ''targetinfo'' requires targetInfo.'); end
    case 'distortion'
        if isempty(targetDist)
            error('Mode ''targetdist'' requires both targetInfo and targetDist.');
        end
    otherwise
        error('useMode must be ''lambda'', ''targetinfo'', or ''targetdist''.');
end

rng(1);

% Get linear constraints
[A_total, b_total, M_triRGBc2diRGBc] = buildGamutConstraints(triLMSCalFormat, dichromatType, Disp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Just reached optimization')

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...  
    'ConstraintTolerance',  1e-10, ...  
    'StepTolerance',        1e-10, ... 
    'Display','iter', ...
    'MaxIterations',200);


%%% NEW! Attempt to add in nonlinear constraint where we keep info a constant
%%% value and search over distortion values (want to keep info as the
%%% target info and then get the transformation that minimizes distortion
switch lower(target)

    case 'lambda'
        % Trade off info vs distortion via lambda
        fun     = @(t_vec) lossFunction('lambda', lambda, t_vec, ...
                        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp);
        nonlcon = [];

    case 'info'
        % Minimize distortion subject to info constraint (=targetInfo)
        epsInfo = 1e-6;

        % Minimize distortion only (lambda=0)
        fun = @(t_vec) lossFunction('lambda', 0.0, t_vec, ...
                        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp);
        % Nonlinear constraint on info
        nonlcon = @(t_vec) constraintBand( ...
                        t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp, ...
                        'info', targetInfo, epsInfo);

    case 'distortion'
        % Maximize info subject to distortion constraint (=targetDist)
        % epsDist = 1e-8;
        epsDist = max(1e-6, 0.02*abs(targetDist));


        % Maximize info only (lambda=1)
        fun = @(t_vec) lossFunction('lambda', 1.0, t_vec, ...
                        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp);
        % Nonlinear constraint on distortion
        nonlcon = @(t_vec) constraintBand( ...
                        t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp, ...
                        'distortion', targetDist, epsDist);
        % options = optimoptions(options,'Algorithm','sqp');  % NEW

end

% Now do the minimization with that 
[transformRGB_opt, fval] = fmincon(fun, T_init(:), ...            
    A_total, b_total, ...           % Linear inequality constraints (gamut)
    [], [], [], [], ...             % Aeq, beq, lb, ub 
    nonlcon, ...                    % Nonlinear constraint
    options);

% Check that the constraint worked
% if ~isempty(nonlcon)
%     [cchk, ~] = nonlcon(transformRGB_opt);
%     if any(cchk > 1e-8)
%         warning('Info constraint not fully satisfied (max c = %.3g).', max(cchk));
%     end
% end

%%% OLD CODE FOR OPTIMIZATION:
% % Optimization - start with identity transformation matrix
% options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-10, 'Display', 'iter','MaxIterations',200);
% % fmincon
% [transformRGB_opt, fval] = fmincon(@(transformRGB) lossFunction(useLambdaOrTargetInfo,lambdaOrTargetInfo,transformRGB, triLMSCalFormat,imgParams,dichromatType,infoFnc,distortionFcn,infoNormalizer, distortionNormalizer, Disp), ...
%     T_init(:), A_total, b_total, [], [], [], [], [], options);
% % Test loss function with final transformation matrix values
% [loss, info, infoNormalized, distortion, distortionNormalized] = lossFunction(useLambdaOrTargetInfo,lambdaOrTargetInfo,transformRGB_opt, triLMSCalFormat,imgParams,dichromatType,infoFnc,distortionFcn,infoNormalizer, distortionNormalizer, Disp);

disp('Just finished optimization')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the info and distortion values for the transformation matrix that
% ultimately was chosen after the optimization
[~, info, infoNormalized, distortion, distortionNormalized] = ...
    lossFunction('lambda', 0, transformRGB_opt, ...
        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
        infoNormalizer, distortionNormalizer, Disp);

% See if constraint worked
bCheck = A_total*transformRGB_opt(:);
if (any(bCheck > b_total))
    fprintf('Failed to satisfy constraint\n');
end  
 
% Reshape optimal solution into matrix
transformRGBmatrix = reshape(transformRGB_opt, 3, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 7: TRANSFORM CONTRAST IMAGE W OPTIMAL %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert LMS to RGB
triRGBCalFormat = Disp.M_cones2rgb * triLMSCalFormat;

% Convert RGB to contrast RGB
triRGBContrastCalFormat = (triRGBCalFormat - Disp.grayRGB) ./ Disp.grayRGB;

% Transform RGB contrast image
triRGBContrastCalFormat_T = transformRGBmatrix * triRGBContrastCalFormat;
diRGBContrastCalFormat_T  = M_triRGBc2diRGBc * triRGBContrastCalFormat_T;

% Add back in gray before outputting the image
% min(min(diRGBCalFormat_T(:)))
% max(max(diRGBCalFormat_T(:)))
diRGBCalFormatOpt = (diRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;
rgbLinDaltonizedCalFormat = (triRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;

% Cut off values outside of gamut, when there is some weird numerical out
% of bounds 
if (max(rgbLinDaltonizedCalFormat(:))>1)% && max(trirgbLinCalFormat_T(:))<1+1e-2)
    rgbLinDaltonizedCalFormat(rgbLinDaltonizedCalFormat>1)=1;
end

if (min(rgbLinDaltonizedCalFormat(:))<0)% && min(trirgbLinCalFormat_T(:))>0-1e-2)
    rgbLinDaltonizedCalFormat(rgbLinDaltonizedCalFormat<0)=0;
end

triRGBCalFormatOpt = rgbLin2RGB(rgbLinDaltonizedCalFormat,Disp);

% Get LMS values to output
LMSDaltonizedCalFormat = Disp.M_rgb2cones * rgbLinDaltonizedCalFormat;

end

% function [c,ceq] = infoBandConstraint(t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
%     infoNormalizer, distortionNormalizer, Disp, targetInfo, epsInfo)
% 
% [~, ~, infoNorm_here] = lossFunction('lambda', 0.0, t_vec, ...
%     triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
%     infoNormalizer, distortionNormalizer, Disp);
% 
% c   = (infoNorm_here - targetInfo).^2 - (epsInfo.^2);
% ceq = [];
% end

% function [c,ceq] = infoBandConstraint(t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
%     infoNormalizer, distortionNormalizer, Disp, targetInfo, epsInfo)
% 
% [~, ~, infoNorm_here,distortionNorm_here] = lossFunction('lambda', 0.0, t_vec, ...
%     triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
%     infoNormalizer, distortionNormalizer, Disp);
% 
% c   = (infoNorm_here - targetInfo).^2 - (epsInfo.^2);
% ceq = [];
% end

function [c, ceq] = constraintBand(t_vec, triLMSCalFormat, imgParams, dichromatType,infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, whichMetric, targetVal, epsBand)
% Inputs:
%   whichMetric : 'info' or 'distortion' (case-insensitive)
%   targetVal   : desired normalized metric value
%   epsBand     : half-width of the allowed band around targetVal
%
% Output for fmincon:
%   c   : inequality vector, must be <= 0
%   ceq : equality vector 

    % Doesn't really matter what the loss is here, I just want to grab the
    % info and distortion norm vals from the loss function
    [~, ~, infoNorm, distNorm] = lossFunction('lambda', 0.0, t_vec, ...
        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
        infoNormalizer, distortionNormalizer, Disp);

    switch lower(whichMetric)
        case {'info','information'}
            metricVal = infoNorm;
        case {'dist','distortion'}
            metricVal = distNorm;
        otherwise
            error('constraintBand: whichMetric must be ''info'' or ''distortion''.');
    end

    % Band inequality: (metric - target)^2 - eps^2 <= 0
    c   = (metricVal - targetVal).^2 - (epsBand.^2);
    ceq = [];
end
