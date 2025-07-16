function [triLMSCalFormatOpt,s_raw, v_raw, s_bal, v_bal, transformRGBmatrix_opt] = colorCorrectionOptimize(useLambdaOrTargetInfo,lambdaOrTargetInfo,triLMSCalFormat, renderType,infoType,distortionType, T_prev, Disp, imgParams, V0,V1)
% Optimizes a linear transformation to enhance color contrast for dichromats
%
% Syntax:
%   [triLMSCalFormatOpt, s_raw, v_raw, s_bal, v_bal, transformRGBmatrix_opt] = ...
%       colorCorrectionOptimize(lambdaOrVar, var, lambdaOrTargetInfo, triLMSCalFormat, ...
%       renderType, infoType, distortionType, constraintWL, T_prev, Disp, V0, V1)
%
% Inputs:
%   useLambdaOrTargetInfo: String. Optimization mode:
%                           'lambda'       – vary tradeoff weight λ (0 to 1)
%                           'targetInfo'   – vary based on target contrast info  
%
%   lambdaOrTargetInfo:   Either the lambda or the var value, depending on useLambdaOrTargetInfo 
%                                 If 'lambdaOrTargetInfo' = 'lambda': (0 ≤ λ ≤ 1). 
%                                    Tradeoff weight balancing contrast maximization
%                                    and similarity to the original image.
%                                 If 'lambdaOrTargetInfo' = 'targetInfo', this specifies
%                                    a contrast info target
%
%   triLMSCalFormat:    Nx3 matrix of original LMS-calibrated image data to be transformed.
%
%   renderType:         String. Type of color vision deficiency:
%                           'Deuteranopia' (M-cone missing)
%                           'Protanopia'   (L-cone missing)
%                           'Tritanopia'   (S-cone missing)
%
%   infoType:       String. Method used to define perceptual contrast variance.
%
%   distortionType:     String. Type of similarity metric to preserve naturalness.
%                           "squared"   -> sum of squared error
%                           "luminance" -> keep chromaticity similar only
%
%   constraintWL:       Scalar. Wavelength (in nm) used to define the confusion plane
%                       for dichromat projection (e.g., 585 for Deuteranopia).
%
%   T_prev:             3×3 matrix. Initial RGB transformation matrix. Usually start with eye(3).
%                       Useful for warm-starting optimization across a λ sweep.
%
%   Disp:               Struct with display parameters. Must include:
%                           .M_cones2rgb  – Matrix to convert LMS to RGB
%                           .T_cones      – Cone fundamentals
%                           .P_monitor    – Monitor spectral power distribution
%                           .wls          – Wavelength sampling
%
%   V0, V1:             Optional variables for internal optimization bookkeeping (used in
%                       variance/similarity functions).
%
% Outputs:
%   triLMSCalFormatOpt: Transformed LMS-calibrated image optimized for colorblind viewing.
%   s_raw:              Similarity of original image to itself (baseline).
%   v_raw:              Variance of original image (baseline).
%   s_bal:              Similarity of optimized image to original.
%   v_bal:              Variance of optimized image.
%   transformRGBmatrix_opt:
%                       3×3 transformation matrix that maps original RGB values to corrected RGB.
%
% Constraints:
%   - RGB values must be between 0 and 1
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

% History
%   09/14/2024  cmd  Initial go.
%
% Examples:
%{
Disp = loadDisplay()
% Load LMS values for this image
imgParams = buildSetParameters('gray',1,128,128);
[triLMSCalFormat,trirgbLinCalFormat,diLMSCalFormat,dirgbLinCalFormat,pathName] = loadLMSvalues('gray','Deuteranopia',Disp,imgParams);
T_prev = eye(3);
[triLMSOpt, s0, v0, s1, v1, T_opt] = colorCorrectionOptimize('lambda', 0.5, triLMSCalFormat,'Deuteranopia', 'LMdifferenceContrast', 'squared', T_prev, Disp,imgParams,[],[]);
%}

% Sample data to see how code works
if isempty(triLMSCalFormat)
    % Pretend data to visualize what the function does
    N = 100;                         % Number of points
    x = randn(1, N);                 % Random data for the first dimension
    y = 2 * x + randn(1, N) * 0.2;   % Strong correlation with the first dimension
    z = 0.5 * x + 0.3 * randn(1, N); % Weaker correlation with the first dimension
    triLMSCalFormat = [x; y; z];
end

%% Find optimal transformation matrix
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        cbType = 2;
    case 'Protanopia'   % l cone deficiency
        cbType = 1;
    case 'Tritanopia'   % s cone deficiency
        cbType = 3;
end

rng(1);

% Number of pixels
nPix   = size(triLMSCalFormat,2);

% Start building the linear constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 1: GET TRICHROMAT RGB VALUES FROM LMS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triRGBCalFormat = Disp.M_cones2rgb * triLMSCalFormat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 2: GET TRI CONTRAST RGB FROM RGB      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast RGB, trichromat
triRGBContrastCalFormat = (triRGBCalFormat - Disp.grayRGB)./Disp.grayRGB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 3: MAKE BLOCKED DIAGONAL MATRIX OF CONTRAST RGB %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast version: the contrast image A
% [(nPix x 3) x 9] =    [nPix x 3]                 [nPix x 3]              [nPix x 3]
% A_upper = blkdiag(triRGBContrastCalFormat', triRGBContrastCalFormat', triRGBContrastCalFormat');      % Upper constraint blocks

% ALTERNATIVE:
% Attempt at making the [3 x 9] matrix for each pixel, then append every
% pixel on the bottom so that the total matrix is size [(3*nPix) x 9]
constraintA = [];
for i = 1:size(triRGBContrastCalFormat,2)
    constraintA = [constraintA; eye(3)*triRGBContrastCalFormat(1,i), eye(3)*triRGBContrastCalFormat(2,i), eye(3)*triRGBContrastCalFormat(3,i)];
end
A_upper = constraintA;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 4: MAKE LOWER MATRIX JUST NEGATIVE OF UPPER %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lower is always negative of upper
A_lower = -(A_upper);              % for -I * X <= 0
% Create A from upper and lower
A = double([A_upper; A_lower]);

% b defines the gamut. In regular RGB, it must stay between 0 and 1
% b = [ones(nPix * 3, 1); zeros(nPix * 3, 1)];
% % Update b to deal with contrast image
% % (contrastA * T) * grayRGB + grayRGB <= b
% % (contrastA * T)                     <= (b - grayRGB)/grayRGB
% b = (b - grayRGB(1))./grayRGB(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 5: RIGHT HAND SIDE OF CONSTRAINT, POST CONTRAST %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = [ones(nPix * 3, 1); ones(nPix * 3, 1)]; % only for the case where gray is background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 6: BUILD CONSTRAINT FOR DICHROMAT     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix that transforms contrast RGB trichromat input into contrast RGB
% dichromat output --> this also needs to be in gamut
% Contrast RGB, dichromat
[calFormatDiLMS,calFormatDirgbLin,M_triToDi] = DichromSimulateLinear(triLMSCalFormat, renderType, Disp);

%%% Gamut check %%%
% dirgb = inv(Disp.M_rgb2cones) * calFormatDiLMS;
% max(dirgb(:))
% min(dirgb(:))

% Matrix that converts LMS contrast into rgb contrast
M_conesC2rgbC = diag(1./Disp.grayRGB) * inv(Disp.M_rgb2cones) * diag(Disp.grayLMS);
M_all = M_conesC2rgbC * M_triToDi * inv(M_conesC2rgbC); % left apply this to tri rgb contrast image

% NOTE: TRY TO DO THIS OUTISDE OF THE OPT. PROCEDURE TO SPEED IT UP ... JUT
% NEED THE CONSTRAINT WL
% Make big diagonal matrix of Ms to multiply at each pixel
contRGB = M_all * triRGBContrastCalFormat;
RGB = (contRGB.*Disp.grayRGB)+ Disp.grayRGB;
max(RGB(:))
min(RGB(:))

bigM = kron(speye(nPix), M_all);

% Dichromat constraint is just trichromat constraint transformed
% Left multuply the big M transformation matrix
A_di_upper = bigM * constraintA;
% Lower is just negative of upper
A_di_lower = -A_di_upper;
% Combine lower and upper
A_di = double([A_di_upper; A_di_lower]);
% b (right hand side of constraint AT<b) is constructed the same as orig
b_di = [ones(nPix * 3, 1); ones(nPix * 3, 1)];

% Combine trichromat and dichromat constraints
A_total = [A; A_di];
b_total = [b; b_di];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial guess at transformation matrix - start with identity
T_I    = eye(3, 3);
T_prev = T_prev;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Just reached optimization')
% Optimization - start with identity transformation matrix
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'StepTolerance', 1e-10, 'Display', 'iter','MaxIterations',200);
% fmincon
[transformRGB_opt_TI, fval] = fmincon(@(transformRGB) loss_function(useLambdaOrTargetInfo,lambdaOrTargetInfo,transformRGB, triLMSCalFormat, renderType,infoType,distortionType, Disp,imgParams), ...
    T_I(:), A_total, b_total, [], [], [], [], [], options);
% Test loss function with final transformation matrix values
[fValOpt_TI, s_raw, v_raw, s_bal, v_bal] = loss_function(useLambdaOrTargetInfo,lambdaOrTargetInfo,transformRGB_opt_TI, triLMSCalFormat,renderType,infoType,distortionType, Disp,imgParams);

% If previous transformation is the identity, then skip this step
if (~isequal(T_prev,T_I)) && (~isequal(T_prev,eye(3,3)))
    % Optimization - start with previous transformation matrix
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'StepTolerance', 1e-10, 'Display', 'iter','MaxIterations',200);
    % fmincon
    [transformRGB_opt_Tprev, fval] = fmincon(@(transformRGB) loss_function(useLambdaOrTargetInfo,lambdaOrTargetInfo,transformRGB, triLMSCalFormat, renderType,infoType,distortionType, Disp,imgParams), ...
        T_prev(:), A_total, b_total, [], [], [], [], [], options);
    % Test loss function with final transformation matrix values
    [fValOpt_Tprev, s_raw, v_raw, s_bal, v_bal] = loss_function(useLambdaOrTargetInfo,lambdaOrTargetInfo, transformRGB_opt_Tprev, triLMSCalFormat,renderType,infoType,distortionType, Disp,imgParams);
    % Choose the transformation that gets a lower loss
    if fValOpt_TI <= fValOpt_Tprev
        % transformRGB_opt = transformRGB_opt_TI;
        transformRGB_opt = transformRGB_opt_Tprev;
    elseif fValOpt_Tprev  < fValOpt_TI
        transformRGB_opt = transformRGB_opt_Tprev;
    end
else
    transformRGB_opt = transformRGB_opt_TI;
end
disp('Just finished optimization')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See if constraint worked
bCheck = A_total*transformRGB_opt(:);
if (any(bCheck > b_total))
    fprintf('Failed to satisfy constraint\n');
end  
 
% Reshape optimal solution into matrix
transformRGBmatrix_opt = reshape(transformRGB_opt, 3, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 7: TRANSFORM CONTRAST IMAGE W OPTIMAL %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform contrast image
triRGBContrastCalFormat_T = transformRGBmatrix_opt * triRGBContrastCalFormat;
diRGBContrastCalFormat_T = M_all * triRGBContrastCalFormat_T;

% Add back in gray before outputting the image
diRGBCalFormat_T = (diRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;
min(min(diRGBCalFormat_T(:)))
max(max(diRGBCalFormat_T(:)))

triRGBCalFormat_T = (triRGBContrastCalFormat_T.*grayRGB) + grayRGB;
if (max(triRGBCalFormat_T(:))>1)
    triRGBCalFormat_T(triRGBCalFormat_T>1)=1;
end

% Get LMS values to output
triLMSCalFormatOpt = Disp.M_rgb2cones * triRGBCalFormat_T;

%% Functions

% OBJECTIVE FUNCTION
    function [loss, s_raw, v_raw, s_bal, v_bal] = loss_function(useLambdaOrTargetInfo,lambdaOrTargetInfo,t_vec, LMSCalFormat, renderType,infoType, distortionType, Disp,imgParams)
        T = reshape(t_vec, 3, 3);       % RESHAPE x_vec INTO 3x3 MATRIX

        % I - LMS
        % O - linear RGB
        % M - LMS2RGB
        % T - RGB TRANSFORMATION
        % [nPizels x 3]     = [3 x nPixels] x [3 x 3] x [3 x 3]

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transformation matrix tips %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % M_rgb2cones = Disp.T_cones*Disp.P_monitor;
        % M_cones2rgb = inv(M_rgb2cones)
        % M_rgb2cones -->     [3 x nPixels] x [nPixels x 3]  APPLY ON LEFT OF RGB VALUES WHEN IN CAL FORMAT
        % M_cones2rgb --> inv([3 x nPixels] x [nPixels x 3]) APPLY ON LEFT OF LMS VALUES WHEN IN CAL FORMAT
        %
        % Must apply M on RIGHT with a TRANSPOSE when LMS is in cal format TRANSPOSE
        % LMSCalFormatTran * M_cones2rgb' --> converts LMS to RGB values
        % T scales the RGB values
        % altogether, (LMSCalFormatTran * M_cones2rgb' * T) returns scaled RGB values in calFormatTransposed format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Convert into RGB where gray is removed
        RGBCalFormat = Disp.M_cones2rgb * LMSCalFormat;

        % Create original contrast image
        RGBContrastCalFormat     = (RGBCalFormat - Disp.grayRGB)./Disp.grayRGB;
        
        % Convert original LMS to contrast image
        LMSContrastCalFormat = (LMSCalFormat-Disp.grayLMS) ./ Disp.grayLMS;

        %%%%%%%% Transformation on gray subtracted image %%%%%%%%
        newRGBContrastCalFormat = T * RGBContrastCalFormat;

        % Add gray back in
        newRGBCalFormat = (newRGBContrastCalFormat.*Disp.grayRGB) + Disp.grayRGB;

        % Convert to LMS excitations 
        newLMSCalFormat = inv(Disp.M_cones2rgb)*newRGBCalFormat;

        % Convert to LMS contrast 
        newLMSContrastCalFormat =  (newLMSCalFormat-Disp.grayLMS) ./ Disp.grayLMS;

        % Get into cal format
        newLMSCalFormatTran = newLMSCalFormat';
        LMSCalFormatTran    = LMSCalFormat';

        % Which cones are missing? available?
        switch (renderType)
            case 'Deuteranopia' % m cone deficiency
                index = [1 3];
                missingIdx = 2;
            case 'Protanopia'   % l cone deficiency
                index = [2 3];
                missingIdx = 1;
            case 'Tritanopia'   % s cone deficiency
                index = [1 2];
                missingIdx = 3;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%% Info term %%%%%%%%%%%%%%%%%%%%%%%%%%
        switch (infoType)
            case 'LMdifferenceContrast'  % use contrast

                LMSold = LMSContrastCalFormat;
                LMSnew = newLMSContrastCalFormat;
                availableConesContrast_old = LMSold(index,:);
                availableConesContrast_new = LMSnew(index,:);
                Lcontrast_old = LMSold(1,:);
                Mcontrast_old = LMSold(2,:);
                % Compute information available
                info = computeInfo_LMdifference(availableConesContrast_old, availableConesContrast_new,Lcontrast_old,Mcontrast_old);

            case 'regress'  % use LMS contrast

                LMSold = LMSContrastCalFormat;
                LMSnew = newLMSContrastCalFormat;
                availableConesContrast_old = LMSold(index,:);
                availableConesContrast_new = LMSnew(index,:);
                missingConeContrast_old    = LMSold(missingIdx,:);
                % Compute information available
                info = computeInfo_regress(availableConesContrast_old, availableConesContrast_new, missingConeContrast_old);

            case 'delta'   % use LMS contrast

                LMSold = LMSContrastCalFormat;
                LMSnew = newLMSContrastCalFormat;
                availableConesContrast_old = LMSold(index,:);
                availableConesContrast_new = LMSnew(index,:);
                missingConeContrast_old    = LMSold(missingIdx,:);
                % Compute information available
                info = computeInfo_delta(availableConesContrast_old, availableConesContrast_new, missingConeContrast_old);

            case 'newConeVar' % use LMS excitations

                LMSold = LMSCalFormat;
                LMSnew = newLMSCalFormat;
                availableCones_new = LMSnew(index,:);
                % Compute information available
                info = computeInfo_newConeVar(availableCones_new);

        end
               
        % Weight by lambda (use this when trying to find range of variances)
        if strcmp(useLambdaOrTargetInfo,'lambda')
            infoWeighted = lambdaOrTargetInfo*info;
        elseif strcmp(useLambdaOrTargetInfo,'targetInfo')
            infoWeighted = info;
        end

        % Normalization to get info and distortion on the same scale
        infoNormalized = infoWeighted./imgParams.infoNorm;

       
        %%%%%%%%%%%%%%%%%%%%%%%%%% Distortion term %%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISTORTION IN CONTRAST OR REGULAR LMS???
        switch (distortionType)
            case 'squared'

                LMSold = LMSCalFormat;
                LMSnew = newLMSCalFormat;
                % Compute distortion
                distortion = computeDistortion_squared(LMSold, LMSnew);

            case 'luminance'  

                LMSold = LMSCalFormat;
                LMSnew = newLMSCalFormat;
                % Compute distortion
                distortion = computeDistortion_luminance(LMSold, LMSnew);

        end

        % Weight by lambda
        if strcmp(useLambdaOrTargetInfo,'lambda')
            distortion_weighted = (1-lambdaOrTargetInfo)*distortion;
        elseif strcmp(useLambdaOrTargetInfo,'targetInfo')
            distortion_weighted = distortion;
        end

        % Normalization to get info and distortion on the same scale
        distortionNormalized = distortion_weighted./imgParams.distortionNorm;

        % Maybe next we normalize by something else:
        % Compute similarity and variance normalizing constants by taking
        % the original image, getting the daltonized version (L cone
        % deficient as well as S cone deficient) then computing the
        % similarity and variance scores between the original and those.
        % This could be a good reference point for scale.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale loss so that it is small enough to make fmincon happy but not
        % so small that it is unhappy. How to determine this?
        fminconFactor = 1e8;

        % To enforce a certain variance value, put it into the loss function
        if strcmp(useLambdaOrTargetInfo,'targetInfo')

            info_diff = infoNormalized - lambdaOrTargetInfo;
            % Scale the difference by some amount so that fmincon prioritizes it
            info_scalar = 1e20;

            loss = fminconFactor*((info_scalar*(info_diff.^2)));

        else  % Otherwise, just minimize this loss
            % want low loss
            % want +info_diff small
            % want info big, so -info small, min (-)
            % want distortion small
            % when distortion is "squared", it is big positive when bad, so
            % min (+); when distortion is luminance, is 0 when good, is big
            % positive when bad - so min (+)
            loss = fminconFactor*((-infoNormalized) + distortionNormalized);
        end

    end

end
