function [LMSDaltonizedCalFormat, rgbLinDaltonizedCalFormat, transformRGBmatrix, info, infoNormalized, distortion, distortionNormalized] = ...
    colorCorrectionOptimize(target, lambda, targetInfo, targetDist, ...
    triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, varargin)
% Optimizes a linear transformation to enhance color contrast for dichromats
%
% % Syntax:
% [triLMSCalFormatOpt, trirgbLinCalFormatOpt, transformRGBmatrix, info, infoNormalized, distortion, distortionNormalized] = ...
%     colorCorrectionOptimize(useLambdaOrTargetInfo, lambdaOrTargetInfo, ...
%     triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, varargin)
%
% Inputs:
%   target                     Are you looking for a specific 'lambda','info', or 'distortion'?
%   lambda                     Tradeoff weight in [0,1] when target='lambda'.
%                                   lambda=1 means maximize information
%                                   lambda=0 means minimize distortion
%                              
%   targetInfo                 Desired normalized information when target='info'.
%   targetDist                 Desired normalized distortion value when target='distortion'.
%
%   triLMSCalFormat            Trichromatic LMS data in CalFormat
%   imgParams                  Image/metric params used by infoFnc and distortionFcn.
%   dichromatType              Dichromat type
%                                   'Deuteranopia'
%                                   'Protanopia'
%                                   'Tritanopia'.
%   infoFnc                    Handle for the information metric.
%   distortionFcn              Handle for similarity metric to preserve naturalness.
%                                   'squared' (sum of squared error),
%   infoNormalizer             Scalar to normalize the info metric.
%   distortionNormalizer       Scalar to normalize the distortion metric.
%   Disp:                      Struct with display parameters. Must include:
%                                   .M_cones2rgb  – Matrix to convert LMS to RGB
%                                   .T_cones      – Cone fundamentals
%                                   .P_monitor    – Monitor spectral power distribution
%                                   .wls          – Wavelength sampling
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
%   11/03?/2025  cmd  Initial go.
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
        if isempty(lambda), error('target ''lambda'' requires a lambda.'); end
    case 'info'
        if isempty(targetInfo), error('target ''info'' requires target info.'); end
    case 'distortion'
        if isempty(targetDist)
            error('target ''distortion'' requires target distortion');
        end
    otherwise
        error('target must be ''lambda'', ''targetinfo'', or ''targetdist''.');
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


% Attempt to add in nonlinear constraint where we keep info a constant
% value and search over distortion values (want to keep info as the
% target info and then get the transformation that minimizes distortion, or
% keep distortion fixed and maximize info)
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
        epsDist = 1e-4;

        % Maximize info only (lambda=1)
        fun = @(t_vec) lossFunction('lambda', 1.0, t_vec, ...
                        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp);
        % Nonlinear constraint on distortion
        nonlcon = @(t_vec) constraintBand( ...
                        t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp, ...
                        'distortion', targetDist, epsDist);

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
 
if strcmpi(target, 'distortion')
    % Report target vs achieved normalized distortion
    fprintf('[DISTORTION] targetDist = %.6f, achieved = %.6f, diff = %.6f\n', ...
        targetDist, distortionNormalized, distortionNormalized - targetDist);
end
% the idea here was that maybe the constraint for distortion worked but
% keeping it in gamut didn't? or maybe the attempt to keep it in gamut is
% what's making it fail the search? 

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
diRGBCalFormatOpt = (diRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;
rgbLinDaltonizedCalFormat = (triRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;

% Cut off values outside of gamut, when there is some weird numerical out
% of bounds 
% if (max(rgbLinDaltonizedCalFormat(:))>1)% && max(trirgbLinCalFormat_T(:))<1+1e-2)
%     rgbLinDaltonizedCalFormat(rgbLinDaltonizedCalFormat>1)=1;
% end
% 
% if (min(rgbLinDaltonizedCalFormat(:))<0)% && min(trirgbLinCalFormat_T(:))>0-1e-2)
%     rgbLinDaltonizedCalFormat(rgbLinDaltonizedCalFormat<0)=0;
% end

triRGBCalFormatOpt = rgbLin2RGB(rgbLinDaltonizedCalFormat,Disp);

% Get LMS values to output
LMSDaltonizedCalFormat = Disp.M_rgb2cones * rgbLinDaltonizedCalFormat;

end

function [c, ceq] = constraintBand(t_vec, triLMSCalFormat, imgParams, dichromatType,infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, whichMetric, targetVal, epsBand)

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

% Try to keep the distortion or info value as close as possible to the
% target value (epsBand gives some wiggle room)
% Band inequality: (metric - target)^2 - eps^2 <= 0
c   = (metricVal - targetVal).^2 - (epsBand.^2);
ceq = [];
end
