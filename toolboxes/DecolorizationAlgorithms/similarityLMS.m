function similarity = similarityLMS(similarityType,LMS_old,LMS_new)
% Computes similarity metric between original and transformed LMS image
%
% Syntax:
%   
%
% Inputs:
%   similarityType:    type of similarity metric
%                           "angle"    -> cosine similarity
%                           "distance" -> normalized euclidean distance
%                           "squared"  -> sum of squared error
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

switch (similarityType)
    case 'angle'
        similarity = (LMS_new(:)'*LMS_old(:))/(norm(LMS_new(:))*norm(LMS_old(:)));
    case 'distance'
        similarity = -norm(LMS_new(:)-LMS_old(:))/norm(LMS_old(:));
    case "squared"
        similarity = -sum((LMS_new(:)-LMS_old(:)).^2);
    otherwise
        error('Unknown similarity type specified');
end


end
