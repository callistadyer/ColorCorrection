function  [LMSDaltonizedCalFormat, rgbLinDaltonizedCalFormat, LMSDaltonizedRenderedCalFormat,rgbLinDaltonizedRenderedCalFormat, transformRGBmatrix_opt,...
    info, infoNormalized, distortion, distortionNormalized]...
    = compute(obj, LMSCalFormat, imgParams, dichromatType, ...
    useLambdaOrTargetInfo, lambdaOrTargetInfo, varargin)
% Generic compute method for the @daltonize class.
%
% Syntax:
%   [LMSDaltonizedCalFormat, LMSDaltonizedRenderedCalFormat] = compute(obj, ...
%       LMSCalFormat, imgParams, dichromatType, ...
%       useLambdaOrTargetInfo, lambdaOrTargetInfo, varargin);
%
% Description:
%    Compute method for the @daltonize class.  This is basically the parameter search code
%    written sufficiently generically that it can use the info, distortion, and render
%    functions specified when the @dalonize object was created.
%
% Inputs:
%    obj                    - the @daltonze object
%    LMSCalFormat           - the LMS excitations to daltonize in PTB CalFormat (3 by npixels matrix)
%    imgParams              - struct containing ancilliary information about the image.
%    dichromatType          - type of dichromat to daltonize for:
%                                       'Protaniopia'
%                                       'Deuteranopia'
%                                       'Tritanopia'
%    useLambdaOrTargetInfo
%    lambdaOrTargetInfo
%
% Optional key/value input arguments:
%    T_init -> defaults to identity eye(3,3)
%
% Outputs:
%   LMSDaltonizedCalFormat - the daltonized image LMS excitations in PTB CalFormat
%   LMSDaltonizedRenderedCalFormat - the daltonized image rendered for the dichromat LMS excitations in PTB CalFormat
%
% See Also:
%     t_daltonize

% History:
%    2025-07-14  dhb, cmd  Wrote it

% Key value pair
parser = inputParser;
parser.addParameter('T_init', eye(3), @(x) isnumeric(x) && isequal(size(x), [3 3]));
parser.addParameter('options', struct(), @isstruct);  % optional
parser.addParameter('infoNormalizer', [], @isnumeric);
parser.addParameter('distortionNormalizer', [], @isnumeric);
parser.parse(varargin{:});

T_init  = parser.Results.T_init;
options = parser.Results.options;

% Load Disp structure. This contains background LMS excitations
Disp = obj.Disp;

% Convert input image to contrast from passed excitations.  Use background LMS
% excitations in the Disp structure.
LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS)./Disp.grayLMS;

% Get normalizers for this image for the info and and distortion functions
%
% Note that you can set the normalizer to 1 and call the function to get an
% unnormalized value, which is what you need to do to get the normalizing value.
normalizerValueToGetRawValue = 1;

% Also note that for the normalizers, you first need to create
% the second, distorted image to compare to the original. The idea
% behind this is to render the image for all three types of dichromats.
% Then, the new image is the L-cone value obtained from protonopia
% simulation, M-cone value obtained from a deuteranopia simulation, and
% S-cone value obtained from a tritanopia simulation.
%
% Here we need to produce those three simulations to build the second
% image to compare to the original.
[calFormatLMS_prot,~,~] = DichromRenderLinear(LMSCalFormat, 'Protanopia', Disp);
[calFormatLMS_deut,~,~] = DichromRenderLinear(LMSCalFormat, 'Deuteranopia', Disp);
[calFormatLMS_trit,~,~] = DichromRenderLinear(LMSCalFormat, 'Tritanopia', Disp);

% Build new LMS to compare to old LMS
LMSCalFormat_new = [calFormatLMS_prot(1,:); calFormatLMS_deut(2,:); calFormatLMS_trit(3,:)];

% Get contrast LMS for info and distortion functions
LMSContrastCalFormat_new = (LMSCalFormat_new - Disp.grayLMS)./Disp.grayLMS;



if isempty(parser.Results.infoNormalizer)
    % Compute normalizers as usual
    infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
    distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams,          normalizerValueToGetRawValue, obj.distortionParams);

else
    % Use pre-passed normalizers
    infoNormalizer = parser.Results.infoNormalizer;
    distortionNormalizer = parser.Results.distortionNormalizer;
end

% Calculate normalizers
% infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
% distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, imgParams,          normalizerValueToGetRawValue, obj.distortionParams);

% Optimization function
[LMSDaltonizedCalFormat, rgbLinDaltonizedCalFormat, transformRGBmatrix_opt, info, infoNormalized, distortion, distortionNormalized] = ...
    colorCorrectionOptimize( ...
    useLambdaOrTargetInfo, lambdaOrTargetInfo, ... % Choose lambda or target info, then set the value
    LMSCalFormat,imgParams, ...                    % LMS image to be transformed
    dichromatType, ...                             % Dichromat type
    obj.infoFcn, obj.distortionFcn, ...            % info and distortion functions
    infoNormalizer, distortionNormalizer,...       % info and distortion normalizers
    Disp, ...                                      % Display  parameters
    'T_init', T_init);                             % Starting point for transformation matrix search


[LMSDaltonizedRenderedCalFormat,rgbLinDaltonizedRenderedCalFormat,~] = DichromRenderLinear(LMSDaltonizedCalFormat, dichromatType, Disp);

end
