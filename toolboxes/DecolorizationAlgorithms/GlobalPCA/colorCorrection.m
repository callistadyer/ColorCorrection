function [triLMScalFormatCorrected, triLMSmeans]  = colorCorrection(triLMSCalFormat,renderType,Disp,bScale)

% Correcting an image so a dichromat can see color contrasts that she could
% not see otherwise. Correction is happening via a PCA 
%
% Syntax:
%   [triLMScalFormatCorrected, triLMSmeans]  = colorCorrection(triLMSCalFormat,renderType,Disp,bScale)
%
% Description:
%
% Inputs:
%   img:          - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.            
%   originalLMS:  - original LMS values of a trichromat seeing img
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   cone_mean_orig- original mean of the missing cone value
%   bScale:       - Boolean. Scale the image values into display range (1
%                   or 0).  A good idea except for 'gray'.
%
% Outputs:
%   correctedLMS  - new LMS values after PCA manipulation (for trichromat)
%   K_opt         - optimal K scalars 
%   D_mnew        - LMS values before scaling
%   T_mean        - mean of LMS values 
%
% Optional key/value pairs:
%   None
%

% Perform PCA

% Perform the easy or hard way? (easy = matlab pca(), hard = fmincon)
method = 'linTransform';
% PCA (the easy way):
if strcmp(method,'easyPCA')
    coeff = pca(triLMSCalFormat',"Centered",true); % Centered true subtracts mean
    PC2D(:,1) = coeff(:, 1); % First principal component
    PC2D(:,2) = coeff(:, 2); % Second principal component

    % LMS means
    triLMSmeans = mean(triLMSCalFormat,2);

    % Mean subtracted LMS values
    % T_ms = triLMSCalFormat - triLMSmeans; % already subtracted off in "centered true" in PCA

    % Map mean subtracted LMS valud onto two principle components (linear regression)
    triLMSCalFormatOpt = PC2D\triLMSCalFormat;

    % Map mean of LMS onto two principle components (linear regression)
    % D_m = PC2D \ T_mean;
    triLMSCalFormatOpt(1,:) = triLMSCalFormatOpt(1,:);
    triLMSCalFormatOpt(2,:) = triLMSCalFormatOpt(1,:);
    triLMSCalFormatOpt(3,:) = triLMSCalFormatOpt(2,:);
elseif strcmp(method,'hardPCA')
    % PCA (the hard way):
    numPCs = 2;
    lambda_var = 0.8;
    [triLMSCalFormatOpt] = decolorOptimize(triLMSCalFormat,method,numPCs,renderType,0,lambda_var,Disp,bScale);
    triLMSCalFormatOpt(1,:) = triLMSCalFormatOpt(:,1);
    triLMSCalFormatOpt(2,:) = triLMSCalFormatOpt(:,1);
    triLMSCalFormatOpt(3,:) = triLMSCalFormatOpt(:,2);
    triLMSmeans = mean(triLMSCalFormat,2);
elseif strcmp(method,'linTransform')
    % decolorOptimize does mean subtraction, then maximizes variance fmincon 
    % expects x y z dimensions in rows and measurements in columns ie. [3 x 1000]  
    lambda_var = 0.5;
    triLMSmeans = mean(triLMSCalFormat,2);
    [triLMSCalFormatOpt] = decolorOptimize(triLMSCalFormat,method,0,renderType,0,lambda_var,Disp,bScale);

end

%% Find the scaling matrix that maps D_ms onto approximate cone values

if strcmp(method,'easyPCA')

    % Initial guess for K
    initialKvec = [0.1, 0.1, 0.1];

    % Define the constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = []; % Lower bounds
    ub = []; % Upper bounds

    % Options
    options = optimset('fmincon');
    options = optimset(options,'Diagnostics','off','Display','iter','LargeScale','off','Algorithm','active-set');

    % Call fmincon
    [K_optvec, fval] = fmincon(@(kVec) T_EstObjectiveFunction(kVec, triLMSCalFormatOpt, triLMSmeans, renderType, Disp, bScale), initialKvec, A, b, Aeq, beq, lb, ub, [], options);

    K_opt = diag(K_optvec);
    % Display results
    disp('Optimal K:');
    disp(K_optvec);
    disp('Objective Function Value:');
    disp(fval);

    disp(['Note: colorCorrectionPCA scaling in lines 100-120 do NOT work when you dont add mean back in'])
    % T_opt = K_opt * projected_data;
    % Add back in mean when you do easy pca
    T_opt = K_opt * triLMSCalFormatOpt + triLMSmeans;

    triLMScalFormatCorrected = T_opt;
else

triLMScalFormatCorrected = triLMSCalFormatOpt;

% OK let's try scaling here to stay in gamut... ideally we would do this in
% optimization but we cannot because we need ALL cone values to tell if it 
% is in gamut, but we are optimizing one at a time
% [correctedLMS_scaled, k] = scaleInGamut(correctedLMS,Disp,bScale);
% correctedLMS = correctedLMS_scaled;
% inGamutColorCorrectionPCA = checkGamut(correctedLMS,Disp,bScale);
% if inGamutColorCorrectionPCA == 0 
%     [correctedLMS_scaled, k] = scaleInGamut(correctedLMS,Disp,bScale);
%     correctedLMS = correctedLMS_scaled;
%     % error('PCA pushing out of gamut')
% end

end