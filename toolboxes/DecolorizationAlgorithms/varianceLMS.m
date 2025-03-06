function variance = varianceLMS(varianceType,renderType,LMS_old,LMS_new)
% Computes variance metric between original and transformed LMS image
%
% Syntax:
%   
%
% Inputs:
%   varianceType:      type of variance metric
%                           "newConeVar" -> maximize variance on the
%                            available cones of the transformed LMS image
%   renderType:        type of dichromacy
%   LMS_old:           original LMS values
%   LMS_new:           transformed LMS values
%
% Outputs:
%   triLMSCalFormatOpt: Transformed LMS values
%
% Constraints:
%   - RGB values must be between 0 and 1
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        index = [1 3];
    case 'Protanopia'   % l cone deficiency
        index = [2 3];
    case 'Tritanopia'   % s cone deficiency
        index = [1 2];
end

% Variance term
switch (varianceType)
    case 'newConeVar'
        variance = (var(LMS_new(index(1), :)) + var(LMS_new(index(2), :)));
    otherwise
        error('Unknown variance type specified');
end


end
