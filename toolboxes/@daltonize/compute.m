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
    
    % Load Disp structure. This contains background LMS excitations
    Disp = loadDisplay();
    % Or...
    Disp = obj.Disp;

    % Convert input image to contrast from passed excitations.  Use background LMS
    % excitations in the Disp structure.
    LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS)./Disp.grayLMS;
        
    % Get normalizers for this image for the info and and distortion functions
    %
    % Note that you can set the normalizer to 1 and call the function to get an
    % unnormalized value, which is what you need to do to get the normalizing value.
    infoNormalizer       = obj.infoFcn(args{:}, obj.infoParams);

    % Also note that for distortion normalizer, you first need to create
    % the second, distorted image to compare to the original. The idea
    % behind this is to render the image for all three types of dichromats.
    % Then, the new image is the L-cone value obtained from protonopia
    % simulation, M-cone value obtained from a deuteranopia simulation, and
    % S-cone value obtained from a tritanopia simulation.
    %
    % Here we need to produce those three simulations to build the second
    % image to compare to the original:

    distortionNormalizer = obj.distortionFcn(args{:}, obj.distortionParams);

    % Use the render function to get the linear transformation needed to render an LMS
    % image for this type of dichromat, for viewing by a trichromat.  We need this to set
    % up the constraint matrix.
    [calFormatDiLMS,calFormatDirgbLin,M_triToDi] = DichromSimulateLinear(calFormatLMS, cbType, Disp);

    % At this point, you can run the daltonizing search.  But we still need to either pass
    % in lambda, or a value of info to lock in and then minimize distortion with respect
    % to, or a value of distoriont to lock in and then maximize info with respec to.  Need
    % to think about how you'll do that.

    % Since we want this to work generically, probably you pull out of the LMS image the
    % two available cone classes and the one that is not available to the dicrhomat.

             
end
        