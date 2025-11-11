function metricFolder = buildMetricFolderName(infoFcn, infoParams, distortionFcn)
% buildMetricFolderName  Construct a folder name for a (info, distortion) pair
%
% Syntax:
%   metricFolder = buildMetricFolderName(infoFcn, infoParams, distortionFcn)
%
% Description:
%   Returns a folder name that identifies the info metric and distortion metric. 
%   For non-regression info metrics, the folder is "<infoFcnName>__<distortionFcnName>". 
%   For the regression metric 'computeInfo_regress', the parameters are embedded so that different
%   regression setups save to different folders:
%   "<infoFcnName>__<predictingWhat>_from_<predictingFromWhat>__<distortionFcnName>"
%
% Inputs:
%   infoFcn        Function handle used for the information metric.
%                  (e.g., @computeInfo_Wade or @computeInfo_regress)
%   infoParams     Struct of parameters for the info metric. 
%                  For 'computeInfo_regress', this MUST include:
%                      .predictingWhat        (e.g., 'L,M,S')
%                      .predictingFromWhat    (e.g., 'L and S')
%   distortionFcn  Function handle used for the distortion metric.
%                  (e.g., @computeDistortion_squared)
%
% Outputs:
%   metricFolder   folder name
%                  (e.g., 'computeInfo_Wade__computeDistortion_squared'
%                  or 'computeInfo_regress__L_M_S_from_L_and_S__computeDistortion_squared').
%
% Optional key/value pairs:
%   None
%
% Examples:
%{
% Non-regression:
  metricFolder = buildMetricFolderName(@computeInfo_Wade, struct(), @computeDistortion_squared)

% Regression (params embedded in name):
  ip = struct('predictingWhat','L,M,S','predictingFromWhat','L and S');
  metricFolder = buildMetricFolderName(@computeInfo_regress, ip, @computeDistortion_squared)
%}
%
% History:
%   2025-11-07  cmd  Standalone utility with folder-name convention used by sweeps/grids.

    % Convert the function handle for the info function into a string name
    infoFcnName       = func2str(infoFcn);

    % Convert the function handle for the distortion function into a string name
    distortionFcnName = func2str(distortionFcn);

    % Regression case is a little annoying...
    if strcmp(infoFcnName, 'computeInfo_regress')

        % Regression mode needs to know what is being predicted from what
        if ~isfield(infoParams,'predictingWhat') || ~isfield(infoParams,'predictingFromWhat')
            % If either field is missing, throw an error
            error('infoParams must include predictingWhat and predictingFromWhat for computeInfo_regress.');
        end

        % Build a string describing the regression relationship
        % Example: "L-from-M"
        paramsStrRaw = sprintf('%s-from-%s', infoParams.predictingWhat, infoParams.predictingFromWhat);

        % Add some underscores, make it pretty
        paramsString   = regexprep(paramsStrRaw, '[^A-Za-z0-9]+', '_');  
        paramsString   = regexprep(paramsString, '^_+|_+$', '');           

        % Construct the folder name in the form:
        %   "<infoFcnName>__<paramsString>__<distortionFcnName>"
        % Example: "computeInfo_regress__L_from_M__computeDistortion_squared"
        metricFolder = sprintf('%s__%s__%s', infoFcnName, paramsString, distortionFcnName);

    else
        % For non-regression info functions, just join the function names with "__"
        % Example: "computeInfo_Wade__computeDistortion_squared"
        metricFolder = sprintf('%s__%s', infoFcnName, distortionFcnName);
    end
end
