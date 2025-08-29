function [info, infoNormalized] = computeInfo_Wade(LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, dichromatType, normalizingValue, Disp, paramsStruct)
% computeInfo_Wade  
%   Measures how much original L–M contrast is redistributed into the (L+M) and S planes.
%
% Syntax:
%   [info, infoNormalized] = computeInfo_Wade(LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, dichromatType, normalizingValue, Disp, paramsStruct)
%
% Inputs:
%   LMSContrastCalFormat_old:   3 x N original LMS contrast
%   LMSContrastCalFormat_new:   3 x N transformed LMS contrast
%   imgParams:                  struct with image metadata
%   dichromatType:              'Protanopia' | 'Deuteranopia' | 'Tritanopia'
%   normalizingValue:           scalar normalization constant
%   Disp:                       display calibration struct (unused here)
%   paramsStruct:               additional parameters (unused here)
%
% Outputs:
%   info:                       scalar value: how much L–M contrast leaked into L+M and S planes
%   infoNormalized:            info divided by normalizingValue

% L, M, and S cone contrasts
L_old = LMSContrastCalFormat_old(1, :);
M_old = LMSContrastCalFormat_old(2, :);
S_old = LMSContrastCalFormat_old(3, :);

L_new = LMSContrastCalFormat_new(1, :);
M_new = LMSContrastCalFormat_new(2, :);
S_new = LMSContrastCalFormat_new(3, :);

% L-M plane
LM_orig = L_old - M_old;

% Delta in (L+M) and S planes
delta_LplusM = (L_new + M_new) - (L_old + M_old);  % (L+M)new - (L+M)old
delta_S      = S_new - S_old;                      % Snew - Sold

% Get best fitting coefficients that explain how much of the original 
% L–M contrast (LM_orig) appears in the delta(L+M) and delta(S) directions.
% Basically a linear regression with no intercept:
% delta ≈ coefficient × LM_orig.
a = LM_orig' \ delta_LplusM'; % how much L–M shows up in L+M plane
b = LM_orig' \ delta_S';      % how much L–M shows up in S plane

% Total info is sum of both projections
info = a^2 + b^2;

% Normalize
infoNormalized = info / normalizingValue;

end
