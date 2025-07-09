function info = computeInfo_regress(availableConesContrast_old, availableConesContrast_new, missingConeContrast_old)
% computeInfo_regress  Contrast loss weighted by regression residuals from original contrast
%
% Syntax:
%   info = computeInfo_regress(availableConesContrast_old, availableConesContrast_new, missingConeContrast_old)
%
% Inputs:
%   availableConesContrast_old:   2 x N matrix of original LMS contrast (available cones)
%   availableConesContrast_new:   2 x N matrix of transformed LMS contrast (available cones)
%   missingConeContrast_old:      1 x N vector of original LMS contrast (missing cone)
%
% Outputs:
%   info:                         information in image – you want the info
%                                 to increase most in the available cone planes, especially where the
%                                 missing cone could not be well predicted from them.
%
% Example:
%{
   info = computeInfo_regress(LMScontrast_old([1 3],:), LMScontrast_new([1 3],:), LMScontrast_old(2,:));
%}

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

