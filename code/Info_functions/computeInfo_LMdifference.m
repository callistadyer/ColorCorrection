function info = computeInfo_LMdifference(LMSContrastCalFormat_old, LMSContrastCalFormat_new, DichromatType, NormalizingValue, Disp, imgParams, varargin)
% computeInfo_LMdifference  Contrast loss weighted by L–M opponent contrast from the original image
%
% Syntax:
%   info = computeInfo_LMdifference(LMSContrastCalFormat_old, LMSContrastCalFormat_new, DichromatType, NormalizingValue, Disp, imgParams)
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
%   info:                       scalar information value based on L–M contrast loss
%
% Optional key/value pairs:
%   'AlgParams' (struct):       algorithm-specific parameters (optional)
%
% Example:
%{
   info = computeInfo_LMdifference(LMScontrast_old([1 2],:), LMScontrast_new([1 2],:));
%}

% -------------------------------
% Parse optional key/value pairs
% -------------------------------
parser = inputParser;
parser.KeepUnmatched = true;
addParameter(parser, 'AlgParams', struct());  % Default: empty struct
parse(parser, varargin{:});

% Get cone indices
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

end
