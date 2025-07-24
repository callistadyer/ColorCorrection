function [distortion, distortionNormalized] = computeDistortion_squared(LMS_old, LMS_new, imgParams, paramsStruct)
% computeDistortion_squared  Computes sum of squared error between LMS images
%
% Syntax:
%   [distortion, distortionNormalized] = computeDistortion_squared(LMS_old, LMS_new, imgParams, paramsStruct)
%
% Inputs:
%   LMS_old:     3 x N matrix of original LMS values
%   LMS_new:     3 x N matrix of transformed LMS values
%   imgParams:   image parameters
%                      imgParams.distortionNorm  --> value used to normalize distortion function
%   paramsStruct: ???
%
% Outputs:
%   distortion:  Scalar – distortion in image via how much total LMS has changed
%                         (large value = high distortion)
%
% Example:
%   distortion = computeDistortion_squared(LMS_old, LMS_new);

diff = LMS_new(:) - LMS_old(:);       % subtract LMS vectors
distortion = sum(diff.^2);            % Sum of squared differences

distortionNormalized = distortion/imgParams.distortionNorm;

end
