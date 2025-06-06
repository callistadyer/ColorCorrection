function triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,renderType,Disp,bScale)

% Syntax:
%   triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,Disp,bScale)
%
% Description:
%
% Inputs:
%   triLMSCalFormat:     - LMS values in cal format
%   Disp:                - Display parameters
%   bScale:              - Boolean. Scale the image values into display range (1
%                          or 0).  A good idea except for 'gray'.
%
% Outputs:
%   triLMScalFormatCorrected: - Corrected LMS values
%
% Optional key/value pairs:
%   None
%

% Perform PCA
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

%% Find the scaling matrix that maps D_ms onto approximate cone values
scale = 1;
if scale == 1
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

end