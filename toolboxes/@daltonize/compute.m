function  [LMSDaltonizedCalFormat, LMSDaltonizedRenderedCalFormat] = compute(obj, ...
                LMSCalFormat, dichromatType, imgParams, options)
% Generic compute method for the @daltonize class.
%
% Syntax:
%   [LMSDaltonizedCalFormat, LMSDaltonizedRenderedCalFormat] = compute(obj, ...
%               LMSCalFormat, dichromatType, imgParams, options)
%
% Description:
%    Compute method for the @daltonize class.  This is basically the parameter search code
%    written sufficiently generically that it can use the info, distortion, and render
%    functions specified when the @dalonize object was created.
%
% Inputs:
%    obj                            - the @daltonze object  
%    LMSCalFormat         - the LMS excitations to daltonize in PTB CalFormat (3 by npixels matrix)
%    dichromatType          - type of dichromat to daltonize for: 'Protaniopia',
%                                       'Deuteranopia', 'Tritanopia'
%    imgParams               - struct containing ancilliary information about the image.
%
% Optional key/value input arguments:
%    None.
%
% Outputs:
%   LMSDaltonizedCalFormat - the daltonized image LMS excitations in PTB CalFormat
%   LMSDaltonizedRenderedCalFormat - the daltonized image rendered for the dichromat LMS excitations in PTB CalFormat
%
% See Also:
%     t_daltonize

% History:
%    2025-07-14  dhb, cmd  Wrote it
    
    % Convert input image to contrast from passed excitations.  Use background LMS
    % excitations in the Disp structure.

    % Get normalizers for this image for the info and and distortion functions
    %
    % Note that you can set the normalizer to 1 and call the function to get an
    % unnormalized value, which is what you need to do to get the normalzing value.

    % Use the render function to get the linear transformation needed to render an LMS
    % image for this type of dichromat, for viewing by a trichromat.  We need this to set
    % up the constraint matrix.

    % At this point, you can run the daltonizing search.  But we still need to either pass
    % in lambda, or a value of info to lock in and then minimize distortion with respect
    % to, or a value of distoriont to lock in and then maximize info with respec to.  Need
    % to think about how you'll do that.

    % Since we want this to work generically, probably you pull out of the LMS image the
    % two available cone classes and the one that is not available to the dicrhomat.

             
end
        