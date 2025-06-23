
function [triRGBImgFormatCorrected,diRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T] = dichromatCorrection(lambdaOrVar,var,lambda,img,renderType,varianceType,similarityType,plateType,correctionMethod,nSquares,modType,constraintWL,T_prev,V0,V1)
% Transform trichromatic image so that dichromat can see more color
% contrast. Also want to try and preserve some naturalness. This is
% accomplished in colorCorrectionOptimize where we incorporate similarity
% to original in the loss function
%
% Syntax:
%   [triRGBImgFormatCorrected,diRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T] = dichromatCorrection(lambdaOrVar,var,lambda,img,renderType,varianceType,similarityType,plateType,correctionMethod,nSquares,modType,constraintWL,T_prev,V0,V1)
%
% Description:
%
% Inputs:
% 
% T_prev,V0,V1)
%   lambdaOrVar:  - String. Use lambda range or specific variance (computed from lamdbas)  
%                   Optimize using lambda value (between 0 and 1) or var
%                   which samples linspace between the variances of
%                   lambda = 0 and lambda = 1
%                       'lambda'
%                       'var'
%   var:          - Double. Variance value. Leave empty if you are using
%                   lambda. Otherwise, we get the var value by
%                   interpolating between the variances at lambda=0 and
%                   lambda=1. Why are we doing this? We don't get nice
%                   evenly spaced image transformations when we use
%                   linspace(lambda=0, lambda=1, 10). Instead, we search
%                   for 10 images where we linspace between the variances
%                   at those endpoints:
%                   linspace(var(lambda=0),var(lambda=1),10) where var is
%                   taken from v_raw_P (or maybe v_bal_P?? not sure)
%   lambda:       - Weight on the variance term in the optimization [0 1]
%   img:          - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   varianceType: - String. Method for quantifying the variance or contrast
%                   enhancement term in optimization
%                       'LMdifferenceContrast' ** this one is pretty good 
%                       'regress' 
%                       'delta'                
%                       'detail'
%                       'newConeVar' ** original
%
%   similarityType: - type of similarity metric
%                        "squared"   -> sum of squared error
%                        "luminance" -> keep chromaticity similar only
%   plateType:    - Double. Only relevant for ishihara plates. 
%                        1 -> gray with missing cone mod
%                        2 -> background random inside with missing cone mod
%                        3 -> LS background, M inside
%                        4 -> like 2 but constrained between .3 and .7 colors so more room for modulation
%   correctionMethod:   - Color correction method:
%                       'linTransform'
%                       'easyPCA'
%                       'hardPCA'
%   nSquares:     - number of squares in isochromatic plate
%
%   modType       - type of isochromatic plate modulation
%                       'rand'
%                       'M'
%                       'L'
%                       'S'
%   constraintWL  - Wavelength that forms plane with gray that the
%                   dichromat image gets projected onto in
%                   DichromSimulateLinear.m 
%                        585 for deuteronopes
%   T_prev        - Initial transformation matrix for the RGB image 
%                   Start this at T_prev = eye(3,3). Usually this is most
%                   useful when you are looping over lambdas (see below)
%                   because you want to use the previous T solution to
%                   initialize the next optimization. This avoids some
%                   wonky failures in the fmincon routine. 
%   V0:           - variance at lambda = 0. Used for lambdaOrVar = 'var'
%   V1:           - variance at lambda = 1. Used for lambdaOrVar = 'var'
%
% Outputs:
%   triRGBImgFormatCorrected:  - Transformed RGB image after PCA and scaling. Also replaced missing cone as done in other code
%   s_raw_P                    - raw similarity values for current lambda (_P indicates that it is for the modulated image)  
%   v_raw_P                    - raw variance values for current lambda
%   s_bal_P                    - balanced similarity values for current lambda = (1-lambda) * s_raw_P   
%   v_bal_P                    - balanced variance values for current lambda = (lambda) * v_raw_P 
%   T                          - starting point for the transformation matrix. First go should always 
%                                be identity eye(3,3) but when doing loops through lambda or var, we 
%                                use the optimal T that was found from the previous run through the loop
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

% History
%   09/05/2024  cmd  Initial go.
%
% Examples:
%{

% EXAMPLE INPUTS:

lambdaOrVar = 'lambda';      % Use lambda, not var
var = [];                    % Keep empty because above is lambda
lambda = 0.5;                % Lambda value of 0.5
img = 'gray';                % Gray image type (square in the middle)
renderType = 'Deuteranopia'; % Color blind simulation -- M-cone deficient
varianceType = 'LMdifferenceContrast';      % How are we defining "color contrast" in optimization?  
similarityType = 'squared';                 % How are we defining "similarity" in optimization?  
plateType = []; % only for ishihara plates  % Type of ishihara plate 
correctionMethod = 'linTransform';          % Type of color correction. Here we pretty much stick to linTransform 
nSquares = 1;                % For 'gray' image type ONLY --> how many squares?
modType = 'M';               % For 'gray' image type ONLY --> what is that square color?  
constraintWL = 585;          % For simulating dichromacy. This number is good for deuteranopia
T_prev = eye(3,3);           % Previous transformation matrix. The purpose is to give a better starting point for searchin in optimization  
V0 = []; % These are only relavent when youre using var instead of lambda - what two variance values are you interpolating between?                    
V1 = []; % These are only relavent when youre using var instead of lambda - what two variance values are you interpolating between? 

% Call using the above inputs:
[triRGBImgFormatCorrected,diRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T] = dichromatCorrection(lambdaOrVar,var,lambda_var,img,renderType,varianceType,similarityType,plateType,correctionMethod,nSquares,modType,constraintWL,T_prev,V0,V1)


%%%%%%%%%%%%%%%% LOOPING THROUGH LAMBDA %%%%%%%%%%%%%%%%
% TEST: this one has a lambda of 0 so should be no transformation = original image
[triRGB_0, diRGB_0, s_raw_0, v_raw_0, s_bal_0, v_bal_0, T_0] = ...
    dichromatCorrection('lambda',[],0, ...
    'gray','Deuteranopia', ...
    'LMdifferenceContrast','squared',[], ...
    'linTransform',1,'M',585, ...
    eye(3),[],[]);

% TEST: this one has a lambda of 1 so should be max transformation = black and gray or white and gray image (max change in gamut) 
[triRGB_1, diRGB_1, s_raw_1, v_raw_1, s_bal_1, v_bal_1, T_1] = ...
    dichromatCorrection('lambda',[],1, ...
    'gray','Deuteranopia', ...
    'LMdifferenceContrast','squared',[], ...
    'linTransform',1,'M',585, ...
    eye(3),[],[]);

nSteps = 10;
lambda_vals = linspace(0,1,nSteps);
T_prev{1} = eye(3);  % Start with identity

for i = 1:nSteps
    lambda_val = lambda_vals(i);

    [triRGB{i}, diRGB{i}, s_raw_P(i), v_raw_P(i), s_bal_P(i), v_bal_P(i), T_prev{i+1}] = ...
        dichromatCorrection('lambda',[],lambda_val, ...
        'gray','Deuteranopia', ...
        'LMdifferenceContrast','squared',[], ...
        'linTransform',1,'M',585, ...
        T_prev{i},[],[]);
end

% Plot how raw and balanced variance and similarity changes across lambda
lambda = lambda_vals;  % Just alias for readability
figure();
subplot(2,2,1);
plot(lambda, s_raw_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('lambda', 'FontSize', 20); ylabel('similarity', 'FontSize', 20); title('raw', 'FontSize', 25);
subplot(2,2,2);
plot(lambda, s_bal_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('lambda', 'FontSize', 20); ylabel('similarity', 'FontSize', 20); title('balanced', 'FontSize', 25);
subplot(2,2,3);
plot(lambda, v_raw_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('lambda', 'FontSize', 20); ylabel('variance', 'FontSize', 20); title('raw', 'FontSize', 25);
subplot(2,2,4);
plot(lambda, v_bal_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('lambda', 'FontSize', 20); ylabel('variance', 'FontSize', 20); title('balanced', 'FontSize', 25);

figure();
subplot(1,2,1);
plot(lambda, s_raw_P + v_raw_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('lambda', 'FontSize', 20); ylabel('similarity + variance', 'FontSize', 20); title('raw', 'FontSize', 25);
subplot(1,2,2);
plot(lambda, s_bal_P + v_bal_P, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('lambda', 'FontSize', 20); ylabel('similarity + variance', 'FontSize', 20); title('balanced', 'FontSize', 25);



%%%%%%%%%%%%%%%% LOOPING THROUGH VAR %%%%%%%%%%%%%%%% 
you need to have already run the lambda of 0 and 1,
then used linspace to sample between the variances at 0 and 1. Then this
"var" setting essentially finds you solutions with the variances equally
spaced between the variance at lambda = 0 and lambda = 1

% Lambda = 0 → minimal transformation
[~, ~, ~, v_raw_0, ~, ~, ~] = ...
    dichromatCorrection('lambda',[],0, ...
    'gray','Deuteranopia', ...
    'LMdifferenceContrast','squared',[], ...
    'linTransform',1,'M',585, ...
    eye(3),[],[]);
V0 = v_raw_0;

% Lambda = 1 → maximal transformation
[~, ~, ~, v_raw_1, ~, ~, ~] = ...
    dichromatCorrection('lambda',[],1, ...
    'gray','Deuteranopia', ...
    'LMdifferenceContrast','squared',[], ...
    'linTransform',1,'M',585, ...
    eye(3),[],[]);
V1 = v_raw_1;

nSteps = 10;
T_prev{1} = eye(3);  % Start with identity

for i = 1:nSteps
    [triRGB_var{i}, diRGB_var{i}, s_raw_P_var(i), v_raw_P_var(i), ...
     s_bal_P_var(i), v_bal_P_var(i), T_prev{i+1}] = ...
        dichromatCorrection('var',i,[], ...
        'gray','Deuteranopia', ...
        'LMdifferenceContrast','squared',[], ...
        'linTransform',1,'M',585, ...
        T_prev{i},V0,V1);
end

% Plot how raw and balanced variance and similarity changes across var
var_idx = 1:nSteps;
figure();
subplot(2,2,1);
plot(var_idx, s_raw_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('similarity', 'FontSize', 20); title('raw', 'FontSize', 25);
subplot(2,2,2);
plot(var_idx, s_bal_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('similarity', 'FontSize', 20); title('balanced', 'FontSize', 25);
subplot(2,2,3);
plot(var_idx, v_raw_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('variance', 'FontSize', 20); title('raw', 'FontSize', 25);
subplot(2,2,4);
plot(var_idx, v_bal_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('variance', 'FontSize', 20); title('balanced', 'FontSize', 25);

figure();
subplot(1,2,1);
plot(var_idx, s_raw_P_var + v_raw_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('similarity + variance', 'FontSize', 20); title('raw', 'FontSize', 25);
subplot(1,2,2);
plot(var_idx, s_bal_P_var + v_bal_P_var, '-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'white');
xlabel('var index', 'FontSize', 20); ylabel('similarity + variance', 'FontSize', 20); title('balanced', 'FontSize', 25);


%}



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

% Load LMS values for this image
[triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues(img,renderType,modType,nSquares,constraintWL,plateType,Disp);

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

