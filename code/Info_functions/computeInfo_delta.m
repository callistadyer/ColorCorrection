function info = computeInfo_delta(LMSContrastCalFormat_old, LMSContrastCalFormat_new, DichromatType, NormalizingValue, Disp, imgParams, varargin)
% computeInfo_delta  Weighted contrast loss using missing cone contrast as weight
%
% Syntax:
%   info = computeInfo_delta(LMSContrastCalFormat_old, LMSContrastCalFormat_new, DichromatType, NormalizingValue, Disp, imgParams)
%
% Inputs:
%   LMSContrastCalFormat_old:   3 x N matrix of original LMS contrast
%   LMSContrastCalFormat_new:   3 x N matrix of transformed LMS contrast 
%   DichromatType:              type of dichromat for this simulation   
%                                    "Protanopia"
%                                    "Deuteranopia"
%                                    "Tritanopia"
%   NormalizingValue:           (unused here) value used to normalize info function
%   Disp:                       (unused here) display parameters
%   imgParams:                  (unused here) image parameters
%
% Outputs:
%   info:                       scalar value – weighted contrast loss using missing cone contrast
%
% Optional key/value pairs:
%   'AlgParams' (struct):       algorithm-specific parameters (optional)
%
% Example:
%{
   info = computeInfo_delta(LMScontrast_old([1 3],:), LMScontrast_new([1 3],:), LMScontrast_old(2,:));
%}

% -------------------------------
% Parse optional key/value pairs
% -------------------------------
parser = inputParser;
parser.KeepUnmatched = true;
addParameter(parser, 'AlgParams', struct());  % Default: empty struct
parse(parser, varargin{:});

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

end
