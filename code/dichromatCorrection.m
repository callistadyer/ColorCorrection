
function [triRGBImgFormatCorrected,diRGBImgFormatCorrected,s_raw_P,v_raw_P,s_bal_P,v_bal_P,T] = dichromatCorrection(varargin)

% Transform trichromatic image so that dichromat can see more color
% contrast. Also want to try and preserve some naturalness. This is
% accomplished in colorCorrectionOptimize where we incorporate similarity
% to original in the loss function
%
% Syntax:
%   [triRGBImgFormatCorrected, diRGBImgFormatCorrected, s_raw_P, v_raw_P, s_bal_P, v_bal_P, T] = ...
%       dichromatCorrection('Name', Value, ...)
%
% Description:
%   Performs color correction on a trichromatic image for a specified
%   dichromat type using a linear transform or PCA method. Optimizes a
%   tradeoff between increasing color variance and preserving similarity.
%
% Key/Value Inputs (all optional):
%   'lambdaOrVar'     - 'lambda' or 'var'                                     (default: 'lambda')
%   'var'             - Scalar variance index for 'var' mode                  (default: [])
%   'lambda'          - Weight on variance term [0,1]                         (default: 0.5)
%   'img'             - Image identifier string (e.g., 'gray' 'flower1.mat')  (default: 'gray')
%   'renderType'      - Dichromacy type: 
%                              'Deuteranopia'
%                              'Protanopia'
%                              'Tritanopia'                                   (default: 'Deuteranopia')
%   'varianceType'    - Metric for variance/contrast: 
%                              'LMdifferenceContrast'
%                              'regress'
%                                                                             (default: 'LMdifferenceContrast')
%   'similarityType'  - Similarity metric: 
%                              'squared'
%                              'luminance'                                    (default: 'squared')
%   'correctionMethod'- 'linTransform', 'easyPCA', or 'hardPCA'               (default: 'linTransform')
%   'setType'         - Integer specifying image set type                     (default: 1)
%   'modType'         - Modulation type: 'M', 'L', 'S', or 'rand'             (default: 'M')
%   'constraintWL'    - Wavelength (nm) defining confusion line plane         (default: 585)
%   'T_prev'          - Initial 3x3 transformation matrix                     (default: eye(3))
%   'V0'              - Variance at lambda=0 (only for 'var' mode)            (default: [])
%   'V1'              - Variance at lambda=1 (only for 'var' mode)            (default: [])
%
% Outputs:
%   triRGBImgFormatCorrected - Corrected trichromatic RGB image
%   diRGBImgFormatCorrected  - Simulated dichromat RGB image
%   s_raw_P                  - Raw similarity for current lambda/var
%   v_raw_P                  - Raw variance for current lambda/var
%   s_bal_P                  - Balanced similarity = (1-lambda)*s_raw_P
%   v_bal_P                  - Balanced variance  = lambda*v_raw_P
%   T                        - Optimized 3x3 transformation matrix for this run
%
% Examples:
%{
  [triRGB, diRGB, s_raw, v_raw, s_bal, v_bal, T] = dichromatCorrection( ...
      'lambdaOrVar', 'lambda', ...
      'lambda', 0.3, ...
      'img', 'gray', ...
      'renderType', 'Deuteranopia', ...
      'varianceType', 'LMdifferenceContrast', ...
      'similarityType', 'squared');
%}
%{
%% Dichromat Correction Lambda Sweep - Key/Value Version
% TEST: lambda=0 → should give no transformation (original image)
[triRGB_0, diRGB_0, s_raw_0, v_raw_0, s_bal_0, v_bal_0, T_0] = ...
    dichromatCorrection( ...
        'lambdaOrVar', 'lambda', ...
        'lambda', 0, ...
        'img', 'gray', ...
        'renderType', 'Deuteranopia', ...
        'varianceType', 'LMdifferenceContrast', ...
        'similarityType', 'squared', ...
        'correctionMethod', 'linTransform', ...
        'setType', 1, ...
        'modType', 'M', ...
        'constraintWL', 585, ...
        'T_prev', eye(3));

% TEST: lambda=1 → should give maximum transformation (max change in gamut)
[triRGB_1, diRGB_1, s_raw_1, v_raw_1, s_bal_1, v_bal_1, T_1] = ...
    dichromatCorrection( ...
        'lambdaOrVar', 'lambda', ...
        'lambda', 1, ...
        'img', 'gray', ...
        'renderType', 'Deuteranopia', ...
        'varianceType', 'LMdifferenceContrast', ...
        'similarityType', 'squared', ...
        'correctionMethod', 'linTransform', ...
        'setType', 1, ...
        'modType', 'M', ...
        'constraintWL', 585, ...
        'T_prev', eye(3));

% Looping through multiple lambdas
nSteps = 10;
lambda_vals = linspace(0, 1, nSteps);
T_prev{1} = eye(3);  % Start with identity

for i = 1:nSteps
    lambda_val = lambda_vals(i);

    [triRGB{i}, diRGB{i}, s_raw_P(i), v_raw_P(i), s_bal_P(i), v_bal_P(i), T_prev{i+1}] = ...
        dichromatCorrection( ...
            'lambdaOrVar', 'lambda', ...
            'lambda', lambda_val, ...
            'img', 'gray', ...
            'renderType', 'Deuteranopia', ...
            'varianceType', 'LMdifferenceContrast', ...
            'similarityType', 'squared', ...
            'correctionMethod', 'linTransform', ...
            'setType', 1, ...
            'modType', 'M', ...
            'constraintWL', 585, ...
            'T_prev', T_prev{i});
end

%% Plot how similarity and variance change across lambdas
figure();
subplot(2,2,1);
plot(lambda_vals, s_raw_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('\lambda', 'FontSize', 20); ylabel('Similarity', 'FontSize', 20); title('Raw Similarity', 'FontSize', 25);

subplot(2,2,2);
plot(lambda_vals, s_bal_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('\lambda', 'FontSize', 20); ylabel('Similarity', 'FontSize', 20); title('Balanced Similarity', 'FontSize', 25);

subplot(2,2,3);
plot(lambda_vals, v_raw_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('\lambda', 'FontSize', 20); ylabel('Variance', 'FontSize', 20); title('Raw Variance', 'FontSize', 25);

subplot(2,2,4);
plot(lambda_vals, v_bal_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('\lambda', 'FontSize', 20); ylabel('Variance', 'FontSize', 20); title('Balanced Variance', 'FontSize', 25);

figure();
subplot(1,2,1);
plot(lambda_vals, s_raw_P + v_raw_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('\lambda', 'FontSize', 20); ylabel('Similarity + Variance', 'FontSize', 20); title('Raw Sum', 'FontSize', 25);

subplot(1,2,2);
plot(lambda_vals, s_bal_P + v_bal_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('\lambda', 'FontSize', 20); ylabel('Similarity + Variance', 'FontSize', 20); title('Balanced Sum', 'FontSize', 25);

%}
%{
%% Dichromat Correction VAR Sweep - Key/Value Version
% Compute V0: run dichromatCorrection at lambda=0
[~, ~, ~, v_raw_0, ~, ~, ~] = ...
    dichromatCorrection( ...
        'lambdaOrVar', 'lambda', ...
        'lambda', 0, ...
        'img', 'gray', ...
        'renderType', 'Deuteranopia', ...
        'varianceType', 'LMdifferenceContrast', ...
        'similarityType', 'squared', ...
        'correctionMethod', 'linTransform', ...
        'setType', 1, ...
        'modType', 'M', ...
        'constraintWL', 585, ...
        'T_prev', eye(3));
V0 = v_raw_0;

% Compute V1: run dichromatCorrection at lambda=1
[~, ~, ~, v_raw_1, ~, ~, ~] = ...
    dichromatCorrection( ...
        'lambdaOrVar', 'lambda', ...
        'lambda', 1, ...
        'img', 'gray', ...
        'renderType', 'Deuteranopia', ...
        'varianceType', 'LMdifferenceContrast', ...
        'similarityType', 'squared', ...
        'correctionMethod', 'linTransform', ...
        'setType', 1, ...
        'modType', 'M', ...
        'constraintWL', 585, ...
        'T_prev', eye(3));
V1 = v_raw_1;

% Loop through VAR indices
nSteps = 10;
T_prev{1} = eye(3);  % Start with identity

for i = 1:nSteps
    [triRGB_var{i}, diRGB_var{i}, s_raw_P_var(i), v_raw_P_var(i), ...
     s_bal_P_var(i), v_bal_P_var(i), T_prev{i+1}] = ...
        dichromatCorrection( ...
            'lambdaOrVar', 'var', ...
            'var', i, ...
            'img', 'gray', ...
            'renderType', 'Deuteranopia', ...
            'varianceType', 'LMdifferenceContrast', ...
            'similarityType', 'squared', ...
            'correctionMethod', 'linTransform', ...
            'setType', 1, ...
            'modType', 'M', ...
            'constraintWL', 585, ...
            'T_prev', T_prev{i}, ...
            'V0', V0, ...
            'V1', V1);
end

%% Plot how similarity and variance change across VAR
var_idx = 1:nSteps;
figure();
subplot(2,2,1);
plot(var_idx, s_raw_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('Similarity', 'FontSize', 20); title('Raw Similarity', 'FontSize', 25);

subplot(2,2,2);
plot(var_idx, s_bal_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('Similarity', 'FontSize', 20); title('Balanced Similarity', 'FontSize', 25);

subplot(2,2,3);
plot(var_idx, v_raw_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('Variance', 'FontSize', 20); title('Raw Variance', 'FontSize', 25);

subplot(2,2,4);
plot(var_idx, v_bal_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('Variance', 'FontSize', 20); title('Balanced Variance', 'FontSize', 25);

figure();
subplot(1,2,1);
plot(var_idx, s_raw_P_var + v_raw_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('Similarity + Variance', 'FontSize', 20); title('Raw Sum', 'FontSize', 25);

subplot(1,2,2);
plot(var_idx, s_bal_P_var + v_bal_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('Similarity + Variance', 'FontSize', 20); title('Balanced Sum', 'FontSize', 25);

%}

% See also: colorCorrectionOptimize, DichromSimulateLinear
%
% Examples are included within the code
%
% History
%   09/05/2024  cmd  Initial go.
%   06/27/2025  cmd  using key/value pairs
%


% Create input parser
p = inputParser;
p.KeepUnmatched = false; % Optional: error if unknown key passed

% Add key/value parameters with defaults and validators
%%%%%%%%%%%%%% PARAMETER           DEFAULT
p.addParameter('lambdaOrVar',      'lambda',       @(x) ischar(x) || isstring(x));
p.addParameter('var',              [],             @(x) isempty(x) || isnumeric(x));
p.addParameter('lambda',           0.5,            @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('img',              'gray',         @(x) ischar(x) || isstring(x));
p.addParameter('renderType',       'Deuteranopia', @(x) ischar(x) || isstring(x));
p.addParameter('varianceType',     'LMdifferenceContrast', @(x) ischar(x) || isstring(x));
p.addParameter('similarityType',   'squared',      @(x) ischar(x) || isstring(x));
p.addParameter('correctionMethod', 'linTransform', @(x) ischar(x) || isstring(x));
p.addParameter('setType',          1,              @(x) isnumeric(x) && isscalar(x));
p.addParameter('modType',          'M',            @(x) ischar(x) || isstring(x));
p.addParameter('constraintWL',     585,            @(x) isnumeric(x) && isscalar(x));
p.addParameter('T_prev',           eye(3),         @(x) isnumeric(x) && isequal(size(x),[3 3]));
p.addParameter('V0',               [],             @(x) isempty(x) || isnumeric(x));
p.addParameter('V1',               [],             @(x) isempty(x) || isnumeric(x));

% Parse the key/value pairs from varargin
p.parse(varargin{:});

% Assign parsed results to local variables
lambdaOrVar     = p.Results.lambdaOrVar;
var             = p.Results.var;
lambda          = p.Results.lambda;
img             = p.Results.img;
renderType      = p.Results.renderType;
varianceType    = p.Results.varianceType;
similarityType  = p.Results.similarityType;
correctionMethod= p.Results.correctionMethod;
setType         = p.Results.setType;
modType         = p.Results.modType;
constraintWL    = p.Results.constraintWL;
T_prev          = p.Results.T_prev;
V0              = p.Results.V0;
V1              = p.Results.V1;



% Define output directory.
% The key preference gets set up by the TbTb local hook.
projectName = 'ColorCorrection';
myFullPath = mfilename('fullpath');
[myPath,myName] = fileparts(myFullPath);

outputDir = getpref(projectName,'outputDir');
outputSubdir = fullfile(outputDir,myName);
if (~exist(outputSubdir,"dir"))
    mkdir(outputSubdir);
end

% Load display 
Disp = loadDisplay(img);

% WHERE TO PUT THIS?
% setType = 1;

% Load LMS values for this image
[triLMSCalFormat,triRGBCalFormat,Disp] = loadLMSvalues(img,renderType,setType,Disp);
[diLMSCalFormat,M_triToDi]             = DichromSimulateLinear(triLMSCalFormat, Disp.grayLMS,  constraintWL, renderType, Disp);

% Color Correction Algorithm
switch (correctionMethod)
    case 'linTransform'
        % decolorOptimize does mean subtraction, then maximizes variance fmincon
        % expects x y z dimensions in rows and measurements in columns ie. [3 x 1000]
        disp('Entering optimization function');
        [triLMScalFormatCorrected,s_raw, v_raw, s_bal, v_bal, T] = colorCorrectionOptimize(lambdaOrVar,var,lambda,triLMSCalFormat,renderType,varianceType,similarityType,constraintWL,T_prev,Disp,V0,V1);
            % triLMScalFormatCorrected_plate = triLMScalFormatCorrected;
            s_raw_P = s_raw;
            v_raw_P = v_raw;
            s_bal_P = s_bal;
            v_bal_P = v_bal;
    case 'easyPCA'
        triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,renderType,Disp);
        % triLMScalFormatCorrected_plate = colorCorrectionEasyPCA(triLMSCalFormat_plate,renderType,Disp);
    case 'hardPCA'
        numPCs = 2;
        triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat,numPCs,Disp);
        % triLMScalFormatCorrected_plate = colorCorrectionHardPCA(triLMSCalFormat_plate,numPCs,Disp);
end


% Imaging the transformation 
disp('callista!!!!! Need to gamma correct!!!!');

%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%
% Create RGB image from LMS
% Dichromat simulation of original image
diRGBCalFormatOrig = Disp.M_cones2rgb * diLMSCalFormat;
[diRGBCalFormatOrig]        = LMS2RGBCalFormat(diLMSCalFormat, Disp);

% Trichromat simulation of original image
triRGBcalFormatOrig = Disp.M_cones2rgb * triLMSCalFormat;
[triRGBcalFormatOrig]       = LMS2RGBCalFormat(triLMSCalFormat, Disp);

%%%%%%%%%%%%%%% CORRECTED %%%%%%%%%%%%%%%
% Corrected trichromat image
% triRGBcalFormatCorrected = Disp.M_cones2rgb * triLMScalFormatCorrected;
[triRGBcalFormatCorrected]        = LMS2RGBCalFormat(triLMScalFormatCorrected, Disp);

[diLMSCalFormatCorrected,~]        = DichromSimulateLinear(triLMScalFormatCorrected, Disp.grayLMS,  constraintWL, renderType, Disp);
% diRGBCalFormatCorrected            = Disp.M_cones2rgb * diLMSCalFormatCorrected;
diRGBCalFormatCorrected            = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp);


% Transform to RGB image format for viewing:
% original trichromat
triRGBImgFormatOrig              = CalFormatToImage(triRGBcalFormatOrig,Disp.m,Disp.n); % no modulation
% corrected trichromat
triRGBImgFormatCorrected         = CalFormatToImage(triRGBcalFormatCorrected,Disp.m,Disp.n); % no modulation
% original dichromat
diRGBImgFormatOrig               = CalFormatToImage(diRGBCalFormatOrig,Disp.m,Disp.n); % no modulation
% corrected dichromat
diRGBImgFormatCorrected          = CalFormatToImage(diRGBCalFormatCorrected,Disp.m,Disp.n); % no modulation


figure('Position',[161   302   562   552]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
imshow(triRGBImgFormatOrig);
title('trichromat original image');

nexttile
imshow(triRGBImgFormatCorrected);
title('trichromat corrected');

nexttile
imshow(diRGBImgFormatOrig);
title('dichromat original image');

nexttile
imshow(diRGBImgFormatCorrected);
title('dichromat corrected');


if strcmp(lambdaOrVar,'var')
sgtitle(['var = ' num2str(var) ', variance: ' varianceType])
elseif strcmp(lambdaOrVar,'lambda')
sgtitle(['lambdavar = ' num2str(lambda) ', variance: ' varianceType])
end


figure(); tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); nexttile
imshow(triRGBImgFormatCorrected);
title(['trichromat corrected, var = ' num2str(var)]);
nexttile
imshow(diRGBImgFormatCorrected);
title('dichromat corrected');


% Save the data
% save(outputSubdir);

end

