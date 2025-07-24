function [info,infoNormalized] = computeInfo_newConeVar(LMSContrastCalFormat_old, LMSContrastCalFormat_new, dichromatType, Disp, imgParams, paramsStruct)
% computeInfo_newConeVar  Maximize variance on the available cones of the transformed LMS image
%
% Syntax:
%   info = computeInfo_newConeVar(LMSContrastCalFormat_old, LMSContrastCalFormat_new, DichromatType, normalizingValue, Disp, imgParams,paramsStruct)
%
% Inputs:
%   LMSContrastCalFormat_old:   3 x N matrix of original LMS contrast
%   LMSContrastCalFormat_new:   3 x N matrix of transformed LMS contrast 
%   DichromatType:              type of dichromat for this simulation   
%                                    "Protanopia"
%                                    "Deuteranopia"
%                                    "Tritanopia"
%   NormalizingValue:           value used to normalize info function
%   Disp:                       (unused here) display parameters
%   imgParams:                  (unused here) image parameters
%
% Outputs:
%   info:                       Scalar total variance across both available cones
%
% Example:
%{

%}

% Determine available cones
switch DichromatType
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

infoNormalized = info/imgParams.infoNorm;


end
