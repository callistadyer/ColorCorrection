function info = computeInfo_regress(LMSContrastCalFormat_old, LMSContrastCalFormat_new, dichromatType, normalizingValue, Disp, imgParams, paramsStruct)
%  Syntax:
%     info  = computeInfo_regress(LMSContrastCalFormat_old, LMSContrastCalFormat_new, dichromatType, normalizingValue, Disp, imgParams, paramsStruct)
%
% Description:
%
% Inputs:
%   LMSContrastCalFormat_old:     3 x N matrix of original LMS contrast
%   LMSContrastCalFormat_new:   3 x N matrix of transformed LMS contrast 
%   dichromatType:                          type of dichromat for this simulation   
%                                                        "Protanopia"
%                                                        "Deuteranopia"
%                                                        "Tritanopia"
%   normalizingValue:                     value used to normalize info function
%   Disp:                                         display parameters
%   imgParams:                              image parameters
%
% Outputs:
%   info:                       a scalar reporting information in image – you want the info
%                                 to increase most in the available cone planes, especially where the
%                                 missing cone could not be well predicted from them.

switch dichromatType
    case 'Protanopia'
        % Missing L cone (1), available: M (2), S (3)
        missingConesIdx     = 1;
        availableConesIdx   = [2 3];
    case 'Deuteranopia'
        % Missing M cone (2), available: L (1), S (3)
        missingConesIdx     = 2;
        availableConesIdx   = [1 3];
    case 'Tritanopia'
        % Missing S cone (3), available: L (1), M (2)
        missingConesIdx     = 3;
        availableConesIdx   = [1 2];
    otherwise
        error('Unknown DichromatType: %s. Use ''Protanopia'', ''Deuteranopia'', or ''Tritanopia''.', dichromatType);
end


availableConesContrast_old = LMSContrastCalFormat_old(availableConesIdx);
missingConeContrast_old    = LMSContrastCalFormat_old(missingConesIdx);


X = availableConesContrast_old';      % N x 2 (predictors)
y = missingConeContrast_old';         % N x 1 (target)

smallestInterestingContrast = 1e-5;      % Threshold to ignore very small contrast values
X(X < smallestInterestingContrast) = 0; 

if all(X == 0, 'all')                 % If all predictors are zero, can't do regression
    residual = y;                     % ... so use the full missing cone contrast as residual
else
    beta = X \ y;                     % Linear least-squares regression to predict missing cone
    y_hat = X * beta;                 % Predicted missing cone contrast
    residual = y - y_hat;             % Difference between actual and predicted (unexplained stuff)
end

weight = residual'.^2;                % More weight where missing cone is poorly predicted

delta = availableConesContrast_old - availableConesContrast_new;  % Contrast difference in visible cones
info = norm(delta .* weight).^2;     


end

