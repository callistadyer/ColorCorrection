function [triLMScalFormatCorrected] = colorCorrect(method,triLMSCalFormat,renderType,lambda_var,Disp,bScale)

% We just need this function in order to call blockproc correctly. This
% simply calls the appropriate color correction algorithm.

% INPUTS:
%   method          - String specifying the color correction method:
%                     'linTransform' - Linear transformation optimization
%                     'easyPCA'      - Simple PCA-based correction
%                     'hardPCA'      - Hard PCA correction using two principal components
%   triLMSCalFormat - Input LMS cal format 
%   renderType      - Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   lambda_var      - Wavelength-related variable (used in 'linTransform')
%   Disp            - Structure containing display parameters
%   bScale          - Scaling factor (irrelevant usually)
%
% OUTPUT:
%   triLMScalFormatCorrected - Color-corrected LMS calibration format.


switch (method)
    case 'linTransform'
        % decolorOptimize does mean subtraction, then maximizes variance fmincon
        % expects x y z dimensions in rows and measurements in columns ie. [3 x 1000]
        [triLMScalFormatCorrected] = colorCorrectionOptimize(triLMSCalFormat,renderType,lambda_var,Disp,bScale);
    case 'easyPCA'
        triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,renderType,Disp,bScale);
    case 'hardPCA'
        numPCs = 2;
        triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat,numPCs,Disp);
end


end