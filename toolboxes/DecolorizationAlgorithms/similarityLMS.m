function similarityLMS()
% Computes similarity metric between original and transformed LMS image
%
% Syntax:
%   
%
% Inputs:
%   triLMSCalFormat:    original LMS values to be transformed
%   bPLOT:              1 -> Plot the result, 0 -> Don't plot
%   lambda_var:         Weight for maximizing variance (0 <= lambda_var <= 1)
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
        similarity_term_raw = (newLMSContrastCalFormatTran(:)'*LMSCalFormatTran(:))/(norm(newLMSContrastCalFormatTran(:))*norm(LMSCalFormatTran(:)));
    case 'distance'
        similarity_term_raw = norm(newLMSContrastCalFormat(:)-LMSCalFormatTran(:))/norm(LMSCalFormatTran(:));
    otherwise
        error('Unknown similarity type specified');
end


end
