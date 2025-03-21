function [triLMSCalFormatOpt,s_raw, v_raw, s_bal, v_bal] = colorCorrectionOptimize(var, triLMSCalFormat, renderType, lambda_var,Disp)
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

% NOTE: I = LMS values [nPix x 3]
triLMSCalFormatTran = triLMSCalFormat'; % data = [3 x nPix]
% M = transformation matrix [3x3]'
% NOTE: MUST APPLY THIS M MATRIX ON THE LEFT OF CAL FORMAT DATA
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);
% Number of pixels
nPix   = size(triLMSCalFormat,2);

% Constraint matrix (A, includes lots of I iterations) and vector (b)
% triLMSCalFormatTran: trichromat LMS values in [nPix x 3] form
triRGBCalFormatTran = triLMSCalFormat' * M_cones2rgb';
% Perhaps round to get rid of small discrepancies during gray subtraction
triRGBCalFormatTran = round(triRGBCalFormatTran,4);
% Create contrast image
grayRGB = [0.5 0.5 0.5]';
% Contrast RGB
triRGBContrastCalFormatTran = (triRGBCalFormatTran - grayRGB')./grayRGB';

% Linear Constraint for staying in gamut
% Contrast version: the contrast image A
A_upper = blkdiag(triRGBContrastCalFormatTran, triRGBContrastCalFormatTran, triRGBContrastCalFormatTran);      % Upper constraint blocks
% % Non contrast version: regular RGB for A
% A_upper = blkdiag(triRGBCalFormatTran, triRGBCalFormatTran, triRGBCalFormatTran);      % Upper constraint blocks
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
b = [ones(nPix * 3, 1); ones(nPix * 3, 1)]; % only for the case where gray is background

% Initial guess at transformation matrix - start with identity
T0 = eye(3, 3);
% See if trying a different starting point would be helpful?
% T0 = T0 * 0.8;
% T0(T0==0) = .1;
% T0 = T0(:);

% OPTIMIZATION SETUP
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter','MaxIterations',200);
bCheck = A*T0(:);
if (any(bCheck > b))
    fprintf('Failed to satisfy constraint\n');
end

% T0 = [1.3220   -1.2685   -0.6358;
%       3.1295   -3.4495   -1.5253;
%      -7.6800   -2.6694    3.1589];
% T0 = T0*1.0e+03;

% Optimization with linear constraints A,b:
[transformRGB_opt, fval] = fmincon(@(transformRGB) loss_function(var, transformRGB, triLMSCalFormatTran, M_cones2rgb, lambda_var, renderType, Disp), ...
    T0, A, b, [], [], [], [], ...,
    @(transformRGB) nonlin_con(var, transformRGB, triLMSCalFormatTran, M_cones2rgb, cbType, Disp), options);

[fValOpt, s_raw, v_raw, s_bal, v_bal] = loss_function(var,transformRGB_opt, triLMSCalFormatTran, M_cones2rgb, lambda_var,renderType, Disp);
bCheck = A*transformRGB_opt(:);
if (any(bCheck > b))
    fprintf('Failed to satisfy constraint\n');
end

% RESHAPE THE OPTIMAL SOLUTION INTO MATRIX FORM
transformRGBmatrix_opt = reshape(transformRGB_opt, 3, 3);

% Contrast manipulation
grayRGB = [0.5 0.5 0.5]';

% Transform contrast image
% newRGBContrastCalFormatTranContrast_out = newRGBContrastCalFormatTranContrast_out * transformRGBmatrix_opt;
newRGBContrastCalFormatTranContrast_out = triRGBContrastCalFormatTran * transformRGBmatrix_opt;

% Add back in gray before outputting the image
triRGBCalFormatTranOpt = (newRGBContrastCalFormatTranContrast_out.*grayRGB') + grayRGB';
triRGBCalFormatOpt     = triRGBCalFormatTranOpt';

% Check for small perterburbations around 0 and 1 in gamut...
% transformations might have pushed it slightly out
% if (min(triRGBCalFormatOpt(:)) < 0) && (min(triRGBCalFormatOpt(:)) > -.01)
%     triRGBCalFormatOpt(triRGBCalFormatOpt<0) = 0;
% end

% Calculate optimal output from input * optimal transform
% triRGBCalFormatTranOpt = triLMSCalFormatTran * M_cones2rgb' * transformRGBmatrix_opt;

% Get LMS values to output
triLMSCalFormatOpt = M_rgb2cones * triRGBCalFormatOpt;

% Check if is in gamut
% inGamutAfterTransform = checkGamut(triLMSCalFormatOpt,Disp,0);
% if inGamutAfterTransform == 0
%     error(['ERROR: decolorOptimize, constraint is not working, transformation is pushing out of gamut'])
% end

%% Functions

% OBJECTIVE FUNCTION
    function [loss, s_raw, v_raw, s_bal, v_bal] = loss_function(var, t_vec, LMSCalFormatTran, M_cones2rgb, lambda, renderType, Disp)
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Must apply M on RIGHT with a TRANSPOSE when LMS is in cal format TRANSPOSE
        % LMSCalFormatTran * M_cones2rgb' --> converts LMS to RGB values
        % T scales the RGB values
        % altogether, (LMSCalFormatTran * M_cones2rgb' * T) returns scaled RGB values in calFormatTransposed format

        % Weber contrast image
        grayRGB = [0.5 0.5 0.5]';

        % Convert into RGB where gray is removed
        RGBCalFormatTran = LMSCalFormatTran * M_cones2rgb';

        % Create contrast image
        RGBContrastCalFormatTran = (RGBCalFormatTran - grayRGB')./grayRGB';
        % Want to use non contrast version? Uncomment below:
        % RGBContrastCalFormatTran = RGBCalFormatTran;

        % Transformation on gray subtracted image
        newRGBContrastCalFormatTran_noGray = RGBContrastCalFormatTran * T;

        % Add gray back in
        newRGBContrastCalFormatTran = (newRGBContrastCalFormatTran_noGray.*grayRGB') + grayRGB';
        % Convert into LMS
        newLMSContrastCalFormatTran = newRGBContrastCalFormatTran*inv(M_cones2rgb)';

        % % Old code, before contrast manipulation:
        % newRGBContrastCalFormatTran = LMSContrastCalFormatTran * M_cones2rgb' * T;
        % newRGBCalFormatTran = LMSCalFormatTran * M_cones2rgb' * T;
        % % Convert to LMS
        % newLMSContrastCalFormatTran = newRGBContrastCalFormatTran*inv(M_cones2rgb)';
        % newLMSCalFormatTran = newRGBCalFormatTran*inv(M_cones2rgb)';
        % Get into cal format
        % newLMSCalFormat = newLMSCalFormatTran';
        % newRGBCalFormat = newRGBCalFormatTran';
        % Check in gamut?
        % minRGB = min(newRGBCalFormatTran(:));
        % maxRGB = max(newRGBCalFormatTran(:));

        % Get into cal format
        newLMSContrastCalFormat = newLMSContrastCalFormatTran';

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
        totalSimilarity = whiteNoiseSimilarity('squared',Disp);
        % Weight by lambda
        similarity_term = (1-lambda)*similarity_term_raw;
        similarity_term_balance = similarity_term/totalSimilarity;
        % Eventually, try to avoid lambda and choose similarity wisely
        % similarity_term = similarity_term_raw; % Getting rid of lambda for now


        % Loss
        % Scale loss so that it is small enough to make fmincon happy but not
        % so small that it is unhappy.
        % How to determine this?
        fminconFactor = 1e6;
        % fminconFactor = 1e20;

        % Variance range: use this to sample in between extreme variances
        % v0 = 5.0686e-14;
        % v1 = 5.4905e-06;
        %
        % s0 = 1;
        % s1 = 0.9798;
        %
        % similarity_range  = linspace(s1, s0, 10);
        % variance_range    = linspace(v0, v1, 10);

        % Difference between the current variance (var_term_raw) and desired variance (variance_range)
        % var_diff = var_term_raw - variance_range(var);
        % sim_diff = similarity_term_raw - similarity_range(var);

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
        s_raw = similarity_term_raw;
        v_raw = var_term_raw;
        s_bal = similarity_term_balance;
        v_bal = var_term_balance;
    end


% This ensures that dichromat rendering does not go out of gamut. Calls
% DichromatSimulateBrettel and checks RGB values
    function [c, ceq] = nonlin_con(var, t_vec, LMSCalFormatTran, M_cones2rgb, cbType, Disp)

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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Must apply M on RIGHT with a TRANSPOSE when LMS is in cal format TRANSPOSE
        % LMSCalFormatTran * M_cones2rgb' --> converts LMS to RGB values
        % T scales the RGB values
        % altogether, (LMSCalFormatTran * M_cones2rgb' * T) returns scaled RGB values in calFormatTransposed format

        % Weber contrast image
        grayRGB = [0.5 0.5 0.5]';

        % Convert to RGB to subtract gray in RGB space
        RGBCalFormatTran = LMSCalFormatTran * M_cones2rgb';
        % Get contrast image
        RGBContrastCalFormatTran = (RGBCalFormatTran - grayRGB')./grayRGB';
        % Transformation
        newRGBContrastCalFormatTran = RGBContrastCalFormatTran * T;
        % Add back in the gray
        newRGBCalFormatTran = (newRGBContrastCalFormatTran.*grayRGB') + grayRGB';
        newRGBCalFormat = newRGBCalFormatTran';

        % Convert into LMS
        newLMSContrastCalFormatTran = newRGBContrastCalFormatTran*inv(M_cones2rgb)';

        % Once you have v0 and v1, now you can tell the function to search
        % only for the solution with these values
        % Nonlinear equality constraint: variance must match one in my chosen
        % range
        varSpecificNonlinear = 0;
        if varSpecificNonlinear == 1
            % For dichromat plate image
            % Obtain v0 and v1 by running the code with lambda still included.
            % Make lambda = 0 and 1 and collect values for var_term_raw in loss function
            v0 = 5.0686e-14;
            v1 = 5.4905e-06;

            s0 = 1;
            s1 = 0.9798;
            similarity_range    = linspace(s1, s0, 10);
            similarity_term_raw = similarityLMS('angle',LMSCalFormatTran,newLMSContrastCalFormatTran);

            variance_range = linspace(v0, v1, 10);
            var_term_raw   = varianceLMS("newConeVar",renderType,[],newLMSContrastCalFormatTran');

            % ceq = var_term_raw - variance_range(var);
            ceq = similarity_term_raw - similarity_range(var);

            tol = 1e-2;

            if abs(ceq) < tol
                ceq = 0; % Treat it as zero within tolerance
            end
        else
            ceq = [];
        end

        triLMSCalFormatTran_new = newRGBCalFormat'*inv(M_cones2rgb)';

        [diLMSCalFormat] = tri2dichromatLMSCalFormat(triLMSCalFormatTran_new',renderType,Disp,0);

        diRGBLinCalFormat = diLMSCalFormat' * M_cones2rgb';

        % Matrix to convert from rgb to xyz
        % Note: this matrix must be applied on the LEFT!!
        % M_rgb2xyz = Disp.T_xyz*Disp.P_monitor;
        % M_xyz2rgb = inv(M_rgb2xyz);
        %
        % % Linear RGB --> XYZ (Brettel takes in XYZ)
        % triXYZCalFormat = M_rgb2xyz * newRGBCalFormat;
        % % Cal Format --> Image Format
        % triXYZImgFormat = CalFormatToImage(triXYZCalFormat,Disp.m,Disp.n);
        %
        % % Convert with Brettel
        % [diXYZ] = DichromatSimulateBrettel(triXYZImgFormat, cbType, []);
        %
        % % Image Format --> CalFormat
        % diXYZCalFormat = ImageToCalFormat(diXYZ);
        % % XYZ --> Linear RGB
        % diRGBLinCalFormat = M_xyz2rgb * diXYZCalFormat;

        % Nonlin constraints must be <= 0
        % Get c(1) and c(2) from min and max
        % Nonlinear constraints for dichromat
        c(1) = -1*min(diRGBLinCalFormat(:));
        c(2) = max(diRGBLinCalFormat(:))-1;

    end

end
