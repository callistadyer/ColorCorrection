function info = computeInfo_newConeVar(availableCones_new)
% computeInfo_newConeVar  Maximize variance on the available cones of the transformed LMS image
%
% Syntax:
%   variance = computeInfo_newConeVar(availableCones_new)
%
% Inputs:
%   availableCones_new:   Transformed LMS values (available cones only) 2 x N
%
% Outputs:
%   info:                 Scalar total variance across both available cones
%
% Example:
%{
   info = computeInfo_newConeVar(LMS_new([1 3], :));
%}

info = var(availableCones_new(1, :)) + var(availableCones_new(2, :));

end
