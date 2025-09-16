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
% delta_LplusM = (L_new + M_new) - (L_old + M_old);  % (L+M)new - (L+M)old
delta_L = (L_new) - (L_old);  % (L+M)new - (L+M)old

delta_S      = S_new - S_old;                      % Snew - Sold

% Get best fitting coefficients that explain how much of the original 
% L–M contrast (LM_orig) appears in the delta(L+M) and delta(S) directions.
% Basically a linear regression with no intercept:
% delta ≈ coefficient × LM_orig.
% a = LM_orig' \ delta_LplusM'; % how much L–M shows up in L+M plane
% b = LM_orig' \ delta_S';      % how much L–M shows up in S plane

x = LM_orig(:);                               % Nx1
Y = [delta_L(:)  delta_S(:)];            % Nx2


projection_deltaL = dot(x,delta_L(:))./(norm(x));
projection_deltaS = dot(x,delta_S(:))./(norm(x));

info = projection_deltaL.^2 + projection_deltaS.^2;

% ab = x \ Y;                                   % 1×2 slopes [a b]
% a  = ab(1);  % how much L–M shows up in L+M plane
% b  = ab(2);  % how much L–M shows up in S plane


% Total info is sum of both projections
% info = a^2 + b^2;

% Prediction error on LM_orig
% predictionError = sum( (delta_LplusM - a * LM_orig).^2 + (delta_S      - b * LM_orig).^2 );
% Vector length of LM_orig
% LMnorm = norm(LM_orig);
% info = 1 - ((predictionError)./(LMnorm).^2);

% Add penalty for negative b. We want positive b because that means if L-M
% is big, then there is more L, and if there is more L, we probably want
% that to be converted to bluer as opposed to greener. 
% if b<0
%     info = info - b^2;
% end

% Normalize
infoNormalized = info / normalizingValue;

end
