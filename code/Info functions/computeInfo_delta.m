function info = computeInfo_delta(availableConesContrast_old, availableConesContrast_new, missingConeContrast_old)
% computeInfo_delta  Weighted contrast loss using missing cone contrast as weight
%
% Syntax:
%   info = computeContrastLoss_delta(availableConesContrast_old, availableConesContrast_new, missingConeContrast_old)
%
% Inputs:
%   availableConesContrast_old:   2 x N matrix of original LMS contrast (available cones)
%   availableConesContrast_new:   2 x N matrix of transformed LMS contrast (available cones)
%   missingConeContrast_old:      1 x N vector of original LMS contrast (missing cone)
%
% Outputs:
%   info:                         information in image - you want the info
%                                 to increase most in the available cone planes, mostly where there was a
%                                 lot of missing cone contrast
%
% Example:
%{
   info = computeInfo_delta(LMScontrast_old([1 3],:), LMScontrast_new([1 3],:), LMScontrast_old(2,:));
%}

delta = availableConesContrast_old - availableConesContrast_new;  % 2 x N
info = norm(delta .* missingConeContrast_old).^2;

end
