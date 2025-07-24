function [distortion,distortionNormalized] = computeDistortion_luminance(LMS_old, LMS_new, imgParams, paramsStruct)
% computeDistortion_luminance  Computes chromaticity difference to preserve hue
%
% Syntax:
%   [distortion,distortionNormalized] = computeDistortion_luminance(LMS_old, LMS_new, imgParams, paramsStruct)
%
% Inputs:
%   LMS_old:     3 x N matrix of original LMS values
%   LMS_new:     3 x N matrix of transformed LMS values
%   imgParams:   image parameters
%                      imgParams.distortionNorm  --> value used to normalize distortion function
%   paramsStruct: ???
%
% Outputs:
%   distortion:  Scalar – distortion in image via how much chromaticity (L/(L+M+S), M/(L+M+S)) has shifted
%                         (large value = more hue shift)
%
% Example:
%  

% Extract LMS from original image
L_old = LMS_old(1, :);     
M_old = LMS_old(2, :);     
S_old = LMS_old(3, :);     

% Extract LMS from transformed image
L_new = LMS_new(1, :);      
M_new = LMS_new(2, :);     
S_new = LMS_new(3, :);     

sum_old = L_old + M_old + S_old;     % Total luminance per pixel (original)
sum_new = L_new + M_new + S_new;     % Total luminance per pixel (transformed)

sum_old(sum_old == 0) = eps;         % Avoid division by zero
sum_new(sum_new == 0) = eps;

% Compute chromaticity coordinates
l_old = L_old ./ sum_old;            % l = L / (L+M+S)
m_old = M_old ./ sum_old;            % m = M / (L+M+S)

l_new = L_new ./ sum_new;            % l (transformed)
m_new = M_new ./ sum_new;            % m (transformed)

% Squared distance in chromaticity space
chroma_diff = (l_new - l_old).^2 + (m_new - m_old).^2;

distortion = sum(chroma_diff);       % Total chromaticity distortion across pixels

distortionNormalized = distortion/imgParams.distortionNorm;

end
