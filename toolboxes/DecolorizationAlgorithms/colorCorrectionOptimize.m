function [triLMSCalFormatOpt,s_raw, v_raw, s_bal, v_bal] = colorCorrectionOptimize(var, triLMSCalFormat, renderType, lambda_var, constraintWL, Disp)
% Optimizes linear transformation of the original cone values
%
% Syntax:
%   [triLMSCalFormatOpt] = colorCorrectionOptimize(triLMSCalFormat, renderType, lambda_var,Disp,bScale)
%
% Inputs:
%   triLMSCalFormat:    original LMS values to be transformed
%   bPLOT:              1 -> Plot the result, 0 -> Don't plot
%   lambda_var:         Weight for maximizing variance (0 <= lambda_var <= 1)
%
% Outputs:
%   triLMSCalFormatOpt: Transformed LMS values
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
wls = (400:10:700)';
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);
Disp.T_cones   = T_cones;
Disp.d         = d;
Disp.P_monitor = P_monitor;
Disp.wls       = wls;
[triLMScalFormatCorrected] = decolorOptimize([],'Deuteranopia',0.5,Disp,0);
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

% Converting from cones to rgb and vice versa
% NOTE: MUST APPLY THIS M MATRIX ON THE LEFT OF CAL FORMAT DATA
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);

% Create gray
grayRGB = [0.5 0.5 0.5]';
grayLMS = M_rgb2cones*grayRGB;

% Number of pixels
nPix   = size(triLMSCalFormat,2);

% Start building the linear constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 1: GET TRICHROMAT RGB VALUES FROM LMS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triRGBCalFormat = M_cones2rgb * triLMSCalFormat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 2: GET TRI CONTRAST RGB FROM RGB      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast RGB, trichromat
triRGBContrastCalFormat = (triRGBCalFormat - grayRGB)./grayRGB;

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
[calFormatDiLMS,M_triToDi] = DichromSimulateLinear(triLMSCalFormat, grayLMS,  constraintWL, renderType, Disp);
dirgb = inv(M_rgb2cones) * calFormatDiLMS;
max(dirgb(:))
min(dirgb(:))

% Matrix that converts LMS contrast into rgb contrast
% M_rgbC2conesC = diag(1./grayLMS) * M_rgb2cones * diag(grayRGB);
M_conesC2rgbC = diag(1./grayRGB) * inv(M_rgb2cones) * diag(grayLMS);

M_all = M_conesC2rgbC * M_triToDi * inv(M_conesC2rgbC); % left apply this to tri rgb contrast image

% NOTE: TRY TO DO THIS OUTISDE OF THE OPT. PROCEDURE TO SPEED IT UP ... JUT
% NEED THE CONSTRAINT WL
% Make big diagonal matrix of Ms to multiply at each pixel

contRGB = M_all * triRGBContrastCalFormat;
RGB = (contRGB.*grayRGB)+grayRGB;
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
T0 = eye(3, 3);
bCheck = A_total*T0(:);
if (any(bCheck > b_total))
    fprintf('Failed to satisfy constraint\n');
end

% Optimization
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'StepTolerance', 1e-10, 'Display', 'iter','MaxIterations',200);
% fmincon
[transformRGB_opt, fval] = fmincon(@(transformRGB) loss_function(var, transformRGB, triLMSCalFormat, M_cones2rgb, lambda_var, renderType, Disp), ...
    T0(:), A_total, b_total, [], [], [], [], [], options);
% Test loss function with final transformation matrix values
[fValOpt, s_raw, v_raw, s_bal, v_bal] = loss_function(var,transformRGB_opt, triLMSCalFormat, M_cones2rgb, lambda_var,renderType, Disp);

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
diRGBCalFormat_T = (diRGBContrastCalFormat_T.*grayRGB) + grayRGB;
min(min(diRGBCalFormat_T(:)))
max(max(diRGBCalFormat_T(:)))

triRGBCalFormat_T = (triRGBContrastCalFormat_T.*grayRGB) + grayRGB;

% Get LMS values to output
triLMSCalFormatOpt = M_rgb2cones * triRGBCalFormat_T;

%% Functions

% OBJECTIVE FUNCTION
    function [loss, s_raw, v_raw, s_bal, v_bal] = loss_function(var, t_vec, LMSCalFormat, M_cones2rgb, lambda, renderType, Disp)
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


        % Weber contrast image
        grayRGB = [0.5 0.5 0.5]';

        % Convert into RGB where gray is removed
        RGBCalFormat = M_cones2rgb * LMSCalFormat;

        % Create contrast image
        RGBContrastCalFormat     = (RGBCalFormat - grayRGB)./grayRGB;
        % Want to use non contrast version? Uncomment below:
        % RGBContrastCalFormat = RGBCalFormat;

        %%%%%%%% Transformation on gray subtracted image %%%%%%%%
        newRGBContrastCalFormat_noGray = T * RGBContrastCalFormat;

        % Add gray back in
        newRGBContrastCalFormat = (newRGBContrastCalFormat_noGray.*grayRGB) + grayRGB;

        % Convert into LMS
        newLMSContrastCalFormat = inv(M_cones2rgb) * newRGBContrastCalFormat;

        % Get into cal format
        newLMSContrastCalFormatTran = newLMSContrastCalFormat';
        LMSCalFormatTran = LMSCalFormat';
        %%%%%%%% Variance term %%%%%%%%
        var_term_raw = varianceLMS("newConeVar",renderType,[],newLMSContrastCalFormat);
        % Weight by lambda (use this when trying to find range of variances)
        var_term = lambda*var_term_raw;
        % Normalize via total variance in white noise image
        totalVariance = whiteNoiseVariance("newConeVar",renderType,Disp);
        % Normalize by whiteNoiseVariance
        var_term_balance = var_term/totalVariance;
        % Eventually, try to avoid lambda and choose variance wisely
        % var_term = var_term_raw;

        %%%%%%%% Similarity term %%%%%%%%
        similarity_term_raw = similarityLMS('squared',LMSCalFormatTran,newLMSContrastCalFormatTran);
        % similarity_term_raw = similarityLMS('squared',round(LMSCalFormatTran,5),round(newLMSContrastCalFormatTran,5));
        % Normalize by whiteNoiseSimilarity
        totalSimilarity = abs(whiteNoiseSimilarity('squared',Disp));
        % Weight by lambda
        similarity_term = (1-lambda)*similarity_term_raw;
        similarity_term_balance = similarity_term/totalSimilarity;
        % Eventually, try to avoid lambda and choose similarity wisely
        % similarity_term = similarity_term_raw; % Getting rid of lambda for now

        %%%%%%%% Loss %%%%%%%%
        % Scale loss so that it is small enough to make fmincon happy but not
        % so small that it is unhappy. How to determine this?
        % fminconFactor = 1e4;
        fminconFactor = 1e6;

        % To enforce a certain variance value, put it into the loss function
        varSpecificLoss = 0;
        if varSpecificLoss == 1
            % Scale the difference by some amount so that fmincon prioritizes it
            var_scalar = 1e6;
            loss = -fminconFactor*((var_scalar*(var_diff.^2)) + var_term_balance);
            % loss = -fminconFactor*((var_scalar*(var_diff.^2)) + similarity_term);
            % Otherwise, just minimize this loss
        else
            % loss = -fminconFactor*(var_term_balance + similarity_term)
            loss = -fminconFactor*(var_term_balance + similarity_term_balance);
        end
        s_raw = -similarity_term_raw;
        v_raw = var_term_raw;
        s_bal = -similarity_term_balance;
        v_bal = var_term_balance;
    end

end
