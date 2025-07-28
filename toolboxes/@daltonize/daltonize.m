classdef daltonize < handle
    % Define a neuralResponseEngine class
    %
    % Syntax:
    %   theDaltonizer =
    %      daltonize(infoFunctionHandle, infoFunctionParamsStruct, ...
    %                     distortionFunctionHandle, distortionFunctionParamsStruct, ...
    %                     renderFunctionHandle, renderFunctionParamsStruct)
    %
    % Description:
    %    The daltonize object is initialized with three key functions and their parameters
    %    (see Inputs section below).  These each have a standard API that allows different
    %    variants of each to be used.  The object then knows how to call these functions to
    %    daltonize any image, through its compute method.
    %
    % Inputs:
    %    infoFunctionHandle        - Function handle to the function that defines
    %                                          how much information about the missing cone
    %                                          class is carried by the cone classes that are
    %                                          there.
    %    infoFunctionParamsStruct  - Parameter structure that the info function understands.
    %    distortionFunctionHandle  - Function handle to the function that defines
    %                                          how much the daltonized image differs from the
    %                                          original
    %    distortionFunctionParamsStruct - Parameter structure that the distortion function understands.
    %    renderFunctionHandle - Function handle to the function that defines
    %                                          how how to linearly render a trichromatic image
    %                                          to approximate how a dichromat might see it.
    %    renderFunctionParamsStruct - Parameter structure that the render function understands.
    %    Disp                            - Standard display structure that defines what RGB means, the cone
    %                                         fundamentals, etc.
    %
    % Outputs:
    %    The created daltonize object.
    %
    % Optional key/value pairs:
    %    None
    %
    % See Also:
    %    t_daltonize

    % History:
    %    2025-07-14  dhb, cmd  Wrote it

    %% Public properties
    properties

    end

    %% Private properties
    properties (SetAccess=private)
        % Info function and parameters.
        infoFcn
        infoParams

        % Distortion function and parameters
        distortionFcn
        distortionParams

        % Render function and parameters
        renderFcn
        renderParams

        % The display structure
        Disp
    end

    % Public methods
    methods
        % Constructor
        function obj = daltonize(infoFunctionHandle, infoFunctionParamsStruct, ...
                distortionFunctionHandle, distortionFunctionParamsStruct, ...
                renderFunctionHandle, renderFunctionParamsStruct, ...
                Disp, ...
                options)

            % You could do argument checking and/or set optional key values pairs here.
            % arguments
            %
            % end

            % Set the key properties from the passed arguments
            obj.infoFcn = infoFunctionHandle;
            obj.infoParams = infoFunctionParamsStruct;
            obj.distortionFcn = distortionFunctionHandle;
            obj.distortionParams = distortionFunctionParamsStruct;
            obj.renderFcn = renderFunctionHandle;
            obj.renderParams = renderFunctionParamsStruct;
            obj.Disp = Disp;
        end

        % Compute method
        [LMSDaltonizedCalFormat, LMSDaltonizedRenderedCalFormat] = compute(obj, ...
            LMSCalFormat, imgParams, dichromatType, ...
            useLambdaOrTargetInfo, lambdaOrTargetInfo, varargin)
        % Compute info sweep method 
        [triLMSCalFormatOpt, trirgbLinCalFormatOpt,diLMSCalFormatOpt,dirgbLinCalFormatOpt, info, infoNormalized, transformRGBmatrixOpt, targetInfoVals] = computeInfoSweep(obj,...
            LMSCalFormat, imgParams, dichromatType, nSteps)



    end

end