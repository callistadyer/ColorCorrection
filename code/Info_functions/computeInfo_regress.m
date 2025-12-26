function [info,infoNormalized] = computeInfo_regress(LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, dichromatType, normalizingValue, Disp, paramsStruct)
%  Syntax:
%     [info,infoNormalized] = computeInfo_regress(LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, dichromatType, normalizingValue, Disp, paramsStruct)
%
% Description:
%
% Inputs:
%   LMSContrastCalFormat_old:     3 x N matrix of original LMS contrast
%   LMSContrastCalFormat_new:     3 x N matrix of transformed LMS contrast 
%   imgParams:                    image parameters
%   dichromatType:                type of dichromat for this simulation   
%                                           "Protanopia"
%                                           "Deuteranopia"
%                                           "Tritanopia"
%   normalizingValue:             Value for normalizing the output
%   Disp:                         display parameters
%   paramsStruct:
%
% Outputs:
%   info:                         a scalar reporting information in image – you want the info
%                                 to increase most in the available cone planes, especially where the
%                                 missing cone could not be well predicted from them.

if nargin < 7 || isempty(paramsStruct) || ...
   ~isfield(paramsStruct,'predictingWhat') || isempty(paramsStruct.predictingWhat) || ...
   ~isfield(paramsStruct,'predictingFromWhat') || isempty(paramsStruct.predictingFromWhat)

    switch dichromatType
        case 'Protanopia'    % Missing L; predict L from [M S]
            paramsStruct.predictingWhat     = 'L';
            paramsStruct.predictingFromWhat = 'M and S';
        case 'Deuteranopia'  % Missing M; predict M from [L S]
            paramsStruct.predictingWhat     = 'M';
            paramsStruct.predictingFromWhat = 'L and S';
        case 'Tritanopia'    % Missing S; predict S from [L M]
            paramsStruct.predictingWhat     = 'S';
            paramsStruct.predictingFromWhat = 'L and M';
        otherwise
            error('Unknown DichromatType: %s. Use ''Protanopia'', ''Deuteranopia'', or ''Tritanopia''.', dichromatType);
    end
end

L_old = LMSContrastCalFormat_old(1,:);  M_old = LMSContrastCalFormat_old(2,:);  S_old = LMSContrastCalFormat_old(3,:);
L_new = LMSContrastCalFormat_new(1,:);  M_new = LMSContrastCalFormat_new(2,:);  S_new = LMSContrastCalFormat_new(3,:);

% Build target y 
y = buildTarget(paramsStruct.predictingWhat, L_old, M_old, S_old);

% Build predictors X 
X = buildPredictors(paramsStruct.predictingFromWhat, ...
                    L_old, M_old, S_old, L_new, M_new, S_new);

smallestInterestingContrast = 1e-5;
X(abs(X) < smallestInterestingContrast) = 0;

N = size(X,1);
if size(y,1) ~= N && size(y,2) == N
    y = y.';  % transpose 1xN -> Nx1 or 3xN -> Nx3
end
if size(y,1) ~= N
    error('Row mismatch: size(X)=%s, size(y)=%s', mat2str(size(X)), mat2str(size(y)));
end

if all(X == 0, 'all')                 % If all predictors are zero, can't do regression
    residual = y;                     % ... so use the full missing cone contrast as residual
else
    beta = X \ y;                     % Linear least-squares regression to predict missing cone
    y_hat = X * beta;                 % Predicted missing cone contrast
    residual = y - y_hat;             % Difference between actual and predicted (unexplained stuff)
end

% info = -sum(residual.^2, 'all'); % when residual is 0, then info is 0
                                 % when residual is high, then info is negative
                                 % so "good" info will be close to 0 and "bad" info will
                                 % be very negative. Still maximizing info in this way

[distortion, ~] = computeDistortion_DE2000(LMSContrastCalFormat_old, y_hat', imgParams, 1, Disp, paramsStruct);
info = -distortion;

infoNormalized = info/normalizingValue;

end


function y = buildTarget(spec, L, M, S)
spec = strtrim(spec);
switch spec
    case 'L',     y = L;
    case 'M',     y = M;
    case 'S',     y = S;
    case 'L-M',   y = L - M;
    case 'L,M,S', y = [L(:), M(:), S(:)];
    otherwise
        error('Unsupported predictingWhat: %s', spec);
end
end

function X = buildPredictors(spec, L_old, M_old, S_old, L_new, M_new, S_new)
% Return N x P matrix (rows = pixels).
spec = strtrim(spec);
switch spec
    case 'L and S'
        X = [L_new(:), S_new(:)];
    case 'L and M'
        X = [L_new(:), M_new(:)];
    case 'M and S'
        X = [M_new(:), S_new(:)];
    case 'L,M,S'
        X = [L_new(:), M_new(:), S_new(:)];
    case 'deltaL and deltaS'
        X = [(L_new(:)-L_old(:)), (S_new(:)-S_old(:))];
    case 'deltaL+M and deltaS'
        X = [(L_new(:)+M_new(:) - (L_old(:)+M_old(:))), (S_new(:)-S_old(:))];
    otherwise
        error('Unsupported predictingFromWhat: %s', spec);
end
end

