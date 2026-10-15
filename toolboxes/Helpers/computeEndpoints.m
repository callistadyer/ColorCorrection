function [infoNormalized_0, infoNormalized_1, distortionNormalized0, distortionNormalized1] = computeEndpoints(saveSubdir, LMSCalFormat, imgParams, dichromatType, obj, infoNormalizer, distortionNormalizer, Disp)
% computeEndpoints  Compute and save sweep endpoint metrics (lambda = 0 and lambda = 1).
%
% Syntax:
%   [infoNormalized_0, infoNormalized_1, distortionNormalized0, distortionNormalized1] = ...
%       computeEndpoints(saveSubdir, LMSCalFormat, imgParams, dichromatType, ...
%                        obj, infoNormalizer, distortionNormalizer, Disp)
%
% Description:
%   This function computes the endpoint values used to define sweep targets
%   in computeSweep. Specifically, it evaluates the normalized information
%   and distortion metrics at:
%
%       - lambda = 0 (minimize distortion)
%       - lambda = 1 (maximize information)
%
%   Because these endpoint optimizations are expensive and invariant across
%   sweep steps, the results are cached to disk and reused on subsequent runs.
%   If cached values are found, they are loaded automatically and no
%   optimization is performed.
%
%   Endpoints are saved to:
%       <saveSubdir>/endpoints.mat
%
% Inputs:
%   saveSubdir            - Directory where sweep outputs are stored; used
%                           to load/save cached endpoint values
%   LMSCalFormat          - 3 x N LMS values of the input image (CalFormat)
%   imgParams             - Struct with image-related parameters (e.g., m, n)
%   dichromatType         - Type of dichromacy:
%                               'Protanopia'
%                               'Deuteranopia'
%                               'Tritanopia'
%   obj                   - Daltonizer object containing:
%                               .infoFcn
%                               .distortionFcn
%                               .infoParams
%                               .distortionParams
%   infoNormalizer        - Scalar used to normalize the information metric
%   distortionNormalizer  - Scalar used to normalize the distortion metric
%   Disp                  - Display struct used by the rendering and metric
%                           functions
%
% Outputs:
%   infoNormalized_0        - Normalized information value at lambda = 0
%   infoNormalized_1        - Normalized information value at lambda = 1
%   distortionNormalized0  - Normalized distortion value at lambda = 0
%   distortionNormalized1  - Normalized distortion value at lambda = 1



endpointFile = fullfile(saveSubdir, 'endpoints.mat');

% Load if available
if exist(endpointFile, 'file')
    S = load(endpointFile, 'infoNormalized_0', 'infoNormalized_1', 'distortionNormalized0', 'distortionNormalized1');

    if all(isfield(S, {'infoNormalized_0', 'infoNormalized_1','distortionNormalized0', 'distortionNormalized1'}))

        infoNormalized_0       = S.infoNormalized_0;
        infoNormalized_1       = S.infoNormalized_1;
        distortionNormalized0  = S.distortionNormalized0;
        distortionNormalized1  = S.distortionNormalized1;
        return;
    end
end

% Compute endpoints

% Note: when you do the lambda=0 endpoint, you are saying "find me the
% transformation with least distortion, don't even care about info." This
% should return the identity. Usually it does, but occassionally it returns
% something just barely different. This is annoying so we are just going to
% force it to be the identity here by adding this to the call:
% 'T_init', eye(3), 'skipFmincon', true
% See inside colorCorrectionOptimize that this just forced the transformation to be identity
[~,~,~,~, infoNormalized_0, ~, distortionNormalized0] = colorCorrectionOptimize("lambda", 0, [], [], LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, obj.infoParams, obj.distortionParams, infoNormalizer, distortionNormalizer, Disp,'T_init', eye(3), 'skipFmincon', true);

[~,~,~,~, infoNormalized_1, ~, distortionNormalized1] = colorCorrectionOptimize("lambda", 1, [], [], LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, obj.infoParams, obj.distortionParams, infoNormalizer, distortionNormalizer, Disp);

% Save
if ~exist(saveSubdir, 'dir'); mkdir(saveSubdir); end

save(endpointFile, ...
    'infoNormalized_0', 'infoNormalized_1', ...
    'distortionNormalized0', 'distortionNormalized1', ...
    '-v7.3');

end
