function variance = varianceLMS(varianceType,renderType,LMS_old,LMS_new,LMSc_old,LMSc_new)
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
%   LMSc_new           transformed LMS values, contrast
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
        missingIdx = 2;
    case 'Protanopia'   % l cone deficiency
        index = [2 3];
        missingIdx = 1;
    case 'Tritanopia'   % s cone deficiency
        index = [1 2];
        missingIdx = 3;
end

% Variance term
switch (varianceType)
    case 'newConeVar'
        variance = (var(LMS_new(index(1), :)) + var(LMS_new(index(2), :)));
    case 'contrast'
        delta1 = LMSc_old(index(1),:) - LMSc_new(index(1),:);
        delta2 = LMSc_old(index(2),:) - LMSc_new(index(2),:);
        missingCone = LMSc_old(missingIdx,:); % should this be old or new M??
        variance = norm([delta1 .* missingCone; delta2 .* missingCone])^2;
    otherwise
        error('Unknown variance type specified');
end


end
