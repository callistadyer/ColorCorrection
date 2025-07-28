function [info, infoNormalized] = computeInfo_LMdifference(LMSContrastCalFormat_old, LMSContrastCalFormat_new, dichromatType, normalizingValue, Disp, imgParams, paramsStruct)
% computeInfo_LMdifference  Contrast loss weighted by L–M opponent contrast from the original image
%
% Syntax:
%   [info, infoNormalized] = computeInfo_LMdifference(LMSContrastCalFormat_old, LMSContrastCalFormat_new, dichromatType, normalizingValue, Disp, imgParams, paramsStruct)
%
% Inputs:
%   LMSContrastCalFormat_old:   3 x N matrix of original LMS contrast
%   LMSContrastCalFormat_new:   3 x N matrix of transformed LMS contrast 
%   dichromatType:              type of dichromat for this simulation   
%                                    "Protanopia"
%                                    "Deuteranopia"
%                                    "Tritanopia"
%   normalizingValue:           value for normalizing info
%   Disp:                       display parameters
%   imgParams:                  image parameters
%   paramsStruct:
%
% Outputs:
%   info:                       scalar information value based on L–M contrast loss
%
% Example:
%{

%}

% Get cone indices
switch dichromatType
    case 'Protanopia'
        availableConesIdx = [2 3];  % M and S
    case 'Deuteranopia'
        availableConesIdx = [1 3];  % L and S
    case 'Tritanopia'
        availableConesIdx = [1 2];  % L and M
    otherwise
        error('Unknown DichromatType: %s. Use ''Protanopia'', ''Deuteranopia'', or ''Tritanopia''.', DichromatType);
end

% Compute weights from L–M contrast in the original image
L = LMSContrastCalFormat_old(1, :);   % L-cone contrast (always from channel 1)
M = LMSContrastCalFormat_old(2, :);   % M-cone contrast (always from channel 2)
LM_diff = abs(L - M);                 % Opponent contrast
alpha = 2;                            % Emphasis exponent
weight = LM_diff .^ alpha;          

% Compute info
availableCones_old = LMSContrastCalFormat_old(availableConesIdx, :);
availableCones_new = LMSContrastCalFormat_new(availableConesIdx, :);
delta = availableCones_old - availableCones_new;

info = norm(delta .* weight).^2;

% Normalize
infoNormalized = info/normalizingValue;

end
