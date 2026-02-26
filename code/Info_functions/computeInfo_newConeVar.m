function [info,infoNormalized] = computeInfo_newConeVar(LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, dichromatType, normalizingValue, Disp, paramsStruct)
% computeInfo_newConeVar  Maximize variance on the available cones of the transformed LMS image
%
% Syntax:
%   [info,infoNormalized] = computeInfo_newConeVar(LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, dichromatType, normalizingValue, Disp, paramsStruct)
%
% Inputs:
%   LMSContrastCalFormat_old:   3 x N matrix of original LMS contrast
%   LMSContrastCalFormat_new:   3 x N matrix of transformed LMS contrast 
%   imgParams:                  image parameters
%   dichromatType:              type of dichromat for this simulation   
%                                    "Protanopia"
%                                    "Deuteranopia"
%                                    "Tritanopia"
%   normalizingValue:           value used to normalize info function
%   Disp:                       display parameters
%   paramsStruct:
%
% Outputs:
%   info:                       Scalar total variance across both available cones
%
% Example:
%{

%}

% Determine available cones
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

% Compute total variance on the available cones 
availableCones_new = LMSContrastCalFormat_new(availableConesIdx, :);
info = var(availableCones_new(1, :)) + var(availableCones_new(2, :));

infoNormalized = info/normalizingValue;


end
