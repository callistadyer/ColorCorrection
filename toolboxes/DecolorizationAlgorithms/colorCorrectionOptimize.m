function [triLMSCalFormatOpt, trirgbLinCalFormat_T, info, infoNormalized, transformRGBmatrix_opt] = ...
    colorCorrectionOptimize(useLambdaOrTargetInfo, lambdaOrTargetInfo, ...
    triLMSCalFormat, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, imgParams, varargin)
% Optimizes a linear transformation to enhance color contrast for dichromats
%
% % Syntax:
% [triLMSCalFormatOpt, trirgbLinCalFormat_T, info, infoNormalized, transformRGBmatrix_opt] = ...
%     colorCorrectionOptimize(useLambdaOrTargetInfo, lambdaOrTargetInfo, ...
%     triLMSCalFormat, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, imgParams, varargin)

% Inputs:
%   useLambdaOrTargetInfo: String. Optimization mode:
%                           'lambda'       – vary tradeoff weight lambda (0 to 1)
%                           'targetInfo'   – vary based on target contrast info  
%
%   lambdaOrTargetInfo:   Either the lambda or the var value, depending on useLambdaOrTargetInfo 
%                                 If 'lambdaOrTargetInfo' = 'lambda': (0 ≤ lambda ≤ 1). 
%                                    Tradeoff weight balancing contrast maximization
%                                    and similarity to the original image.
%                                 If 'lambdaOrTargetInfo' = 'targetInfo', this specifies
%                                    a contrast info target
%
%   triLMSCalFormat:    Nx3 matrix of original LMS-calibrated image data to be transformed.
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
%   constraintWL:       Scalar. Wavelength (in nm) used to define the confusion plane
%                       for dichromat projection (e.g., 585 for Deuteranopia).
%
%   T_prev:             3×3 matrix. Initial RGB transformation matrix. Usually start with eye(3).
%                       Useful for warm-starting optimization across a lambda sweep.
%
%   Disp:               Struct with display parameters. Must include:
%                           .M_cones2rgb  – Matrix to convert LMS to RGB
%                           .T_cones      – Cone fundamentals
%                           .P_monitor    – Monitor spectral power distribution
%                           .wls          – Wavelength sampling
%
%   imgParams:          Struct with image parameters
%
% Outputs:
%   triLMSCalFormatOpt: Transformed LMS-calibrated image optimized for colorblind viewing.
%   s_raw:              Similarity of original image to itself (baseline).
%   v_raw:              Variance of original image (baseline).
%   s_bal:              Similarity of optimized image to original.
%   v_bal:              Variance of optimized image.
%   transformRGBmatrix_opt:
%                       3×3 transformation matrix that maps original RGB values to corrected RGB.
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
Disp = loadDisplay()
% Load LMS values for this image
imgParams = buildSetParameters('gray',1,128,128);
[triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = loadLMSvalues('gray','Deuteranopia',Disp,imgParams);
T_prev = eye(3);
[triLMSOpt, s0, v0, s1, v1, T_opt] = colorCorrectionOptimize('lambda', 0.5, triLMSCalFormat,'Deuteranopia', 'LMdifferenceContrast', 'squared', T_prev, Disp,imgParams,[],[]);
%}

% Sample data to see how code works
if isempty(triLMSCalFormat)
    % Pretend data to visualize what the function does
    N = 100;                         % Number of points
    x = randn(1, N);                 % Random data for the first dimension
    y = 2 * x + randn(1, N) * 0.2;   % Strong correlation with the first dimension
    z = 0.5 * x + 0.3 * randn(1, N); % Weaker correlation with the first dimension
    triLMSCalFormat = [x; y; z];
end

%% Find optimal transformation matrix

% -------------------------------------------------------------------------
% PARSE OPTIONAL PARAMETERS
% -------------------------------------------------------------------------
parser = inputParser;
parser.addParameter('T_init', eye(3), @(x) isnumeric(x) && isequal(size(x), [3 3]));
parser.parse(varargin{:});
T_init = parser.Results.T_init;

rng(1);

% Get linear constraints
[A_total, b_total, M_triRGBc2diRGBc] = buildGamutConstraints(triLMSCalFormat, dichromatType, Disp);

% Initial guess at transformation matrix - start with identity
T_I    = eye(3, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Just reached optimization')
% Optimization - start with identity transformation matrix
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'StepTolerance', 1e-10, 'Display', 'iter','MaxIterations',200);
% fmincon
[transformRGB_opt_TI, fval] = fmincon(@(transformRGB) lossFunction(useLambdaOrTargetInfo,lambdaOrTargetInfo,transformRGB, triLMSCalFormat, dichromatType,infoFnc,distortionFcn,infoNormalizer, distortionNormalizer, Disp,imgParams), ...
    T_I(:), A_total, b_total, [], [], [], [], [], options);
% Test loss function with final transformation matrix values
[fValOpt_TI, info, infoNormalized] = lossFunction(useLambdaOrTargetInfo,lambdaOrTargetInfo,transformRGB_opt_TI, triLMSCalFormat,dichromatType,infoFnc,distortionFcn,infoNormalizer, distortionNormalizer, Disp,imgParams);
transformRGB_opt = transformRGB_opt_TI;

disp('Just finished optimization')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See if constraint worked
bCheck = A_total*transformRGB_opt(:);
if (any(bCheck > b_total))
    fprintf('Failed to satisfy constraint\n');
end  
 
% Reshape optimal solution into matrix
transformRGBmatrix_opt = reshape(transformRGB_opt, 3, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 7: TRANSFORM CONTRAST IMAGE W OPTIMAL %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert LMS to RGB
triRGBCalFormat = Disp.M_cones2rgb * triLMSCalFormat;

% Convert RGB to contrast RGB
triRGBContrastCalFormat = (triRGBCalFormat - Disp.grayRGB) ./ Disp.grayRGB;

% Transform RGB contrast image
triRGBContrastCalFormat_T = transformRGBmatrix_opt * triRGBContrastCalFormat;
diRGBContrastCalFormat_T  = M_triRGBc2diRGBc * triRGBContrastCalFormat_T;

% Add back in gray before outputting the image
% diRGBCalFormat_T = (diRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;
% min(min(diRGBCalFormat_T(:)))
% max(max(diRGBCalFormat_T(:)))

trirgbLinCalFormat_T = (triRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;
if (max(trirgbLinCalFormat_T(:))>1)
    trirgbLinCalFormat_T(trirgbLinCalFormat_T>1)=1;
end

if (min(trirgbLinCalFormat_T(:))<0)
    trirgbLinCalFormat_T(trirgbLinCalFormat_T<1)=0;
end

triRGBCalFormat_T = rgbLin2RGB(trirgbLinCalFormat_T,Disp);

% Get LMS values to output
triLMSCalFormatOpt = Disp.M_rgb2cones * trirgbLinCalFormat_T;

end
