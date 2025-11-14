function [distortion, distortionNormalized] = computeDistortion_squared(LMS_old, LMS_new, imgParams, normalizingValue, paramsStruct)
% computeDistortion_squared  Computes sum of squared error between LMS images
%
% Syntax:
%   [distortion, distortionNormalized] = computeDistortion_squared(LMS_old, LMS_new, imgParams, normalizingValue, paramsStruct)
%
% Inputs:
%   LMS_old:            3 x N matrix of original LMS values
%   LMS_new:            3 x N matrix of transformed LMS values
%   imgParams:          image parameters
%   normalizingValue:   Value for normalizing the output
%   paramsStruct: 
%
% Outputs:
%   distortion:  Scalar – distortion in image via how much total LMS has changed
%                         (large value = high distortion)
%
% Example:
%   distortion = computeDistortion_squared(LMS_old, LMS_new);

diff = LMS_new(:) - LMS_old(:);       % subtract LMS vectors
distortion = sum(diff.^2);            % Sum of squared differences

distortionNormalized = distortion/normalizingValue;

% distortionNormalized = distortionNormalized/100;
end
