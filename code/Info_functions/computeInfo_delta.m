function [info,infoNormalized] = computeInfo_delta(LMSContrastCalFormat_old, LMSContrastCalFormat_new, dichromatType, Disp, imgParams, paramsStruct)
% computeInfo_delta  Weighted contrast loss using missing cone contrast as weight
%
% Syntax:
%   [info,infoNormalized] = computeInfo_delta(LMSContrastCalFormat_old, LMSContrastCalFormat_new, dichromatType, Disp, imgParams, paramsStruct)
%
% Inputs:
%   LMSContrastCalFormat_old:   3 x N matrix of original LMS contrast
%   LMSContrastCalFormat_new:   3 x N matrix of transformed LMS contrast 
%   dichromatType:              type of dichromat for this simulation   
%                                    "Protanopia"
%                                    "Deuteranopia"
%                                    "Tritanopia"
%   Disp:                       display parameters
%   imgParams:                  image parameters
%                                  imgParams.infoNorm --> value used to normalize info function
%
% Outputs:
%   info:                       scalar value – weighted contrast loss using missing cone contrast
%
% Example:
%{

%}

% Identify available and missing cones
switch DichromatType
    case 'Protanopia'
        missingConeIdx     = 1;
        availableConesIdx  = [2 3];  % M, S
    case 'Deuteranopia'
        missingConeIdx     = 2;
        availableConesIdx  = [1 3];  % L, S
    case 'Tritanopia'
        missingConeIdx     = 3;
        availableConesIdx  = [1 2];  % L, M
    otherwise
        error('Unknown DichromatType: %s. Use ''Protanopia'', ''Deuteranopia'', or ''Tritanopia''.', DichromatType);
end

% Extract contrast values
available_old = LMSContrastCalFormat_old(availableConesIdx, :);
available_new = LMSContrastCalFormat_new(availableConesIdx, :);
missing_old   = LMSContrastCalFormat_old(missingConeIdx, :);  % 1 x N

% Compute info
delta = available_old - available_new;        % 2 x N
info = norm(delta .* missing_old).^2;         % scalar

infoNormalized = info/imgParams.infoNorm;

end
