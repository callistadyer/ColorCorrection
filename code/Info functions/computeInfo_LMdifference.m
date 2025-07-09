function info = computeInfo_LMdifference(availableConesContrast_old, availableConesContrast_new,Lcontrast_old,Mcontrast_old)
% computeInfo_LMdifference  Contrast loss weighted by L-M difference in original contrast
%
% Syntax:
%   info = computeInfo_LMdifference(availableConesContrast_old, availableConesContrast_new)
%
% Inputs:
%   availableConesContrast_old:   2 x N matrix [L; M] of contrast from original LMS
%   availableConesContrast_new:   2 x N matrix [L; M] of contrast from transformed LMS
%
% Outputs:
%   info:                         information in image – you want the info
%                                 to increase most in the available cone planes, especially where
%                                 the L–M difference was large in the original image.
%
% Example:
%{
   info = computeInfo_LMdifference(LMScontrast_old([1 2],:), LMScontrast_new([1 2],:));
%}

L = Lcontrast_old;         % Extract L-cone contrast from the original image
M = Mcontrast_old;         % Extract M-cone contrast from the original image

LM_diff = abs(L - M);                         % Compute L–M opponent contrast (absolute difference)
alpha = 2;                                    % Exponent to emphasize stronger opponent differences
weight = LM_diff .^ alpha;                    

delta = availableConesContrast_old - availableConesContrast_new;  % Contrast difference in visible cones
info = norm(delta .* weight).^2;            


end
