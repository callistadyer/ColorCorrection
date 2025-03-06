function [triLMSCalFormatOpt] = colorCorrectionOptimize(triLMSCalFormat, renderType, lambda_var,Disp)
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
A_upper = blkdiag(triRGBCalFormatTran, triRGBCalFormatTran, triRGBCalFormatTran);      % Upper constraint blocks
A_lower = -A_upper;              % for -I * X <= 0
A = double([A_upper; A_lower]);
b = [ones(nPix * 3, 1); zeros(nPix * 3, 1)];

% Initial guess at transformation matrix - start with identity
T0 = eye(3, 3);
% See if trying a different starting point would be helpful?
% T0 = T0 * 0.8;
% T0(T0==0) = .1;
% T0 = T0(:);

% OPTIMIZATION SETUP
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter','MaxIterations',70);
[transformRGB_opt, fval] = fmincon(@(transformRGB) loss_function(transformRGB, triLMSCalFormatTran, M_cones2rgb, lambda_var, renderType, Disp), ...
    T0, A, b, [], [], [], [], ...,
    @(transformRGB) nonlin_con(transformRGB, triLMSCalFormatTran, M_cones2rgb, cbType, Disp), options);

fValOpt = loss_function(transformRGB_opt, triLMSCalFormatTran, M_cones2rgb, lambda_var,renderType, Disp);
bCheck = A*transformRGB_opt(:);
if (any(bCheck > b))
    fprintf('Failed to satisfy constraint\n');
end

% RESHAPE THE OPTIMAL SOLUTION INTO MATRIX FORM
transformRGBmatrix_opt = reshape(transformRGB_opt, 3, 3);

% Contrast manipulation
grayRGB = [0.5 0.5 0.5]';

% Transform into RGB space so you can subtract out gray in RGB (doesn't
% work exactly right when you subtract in LMS... some odd rounding errors)
newRGBCalFormatTran_out = triLMSCalFormatTran * M_cones2rgb';

% Get contrast image by subtracting out gray
newRGBContrastCalFormatTranContrast_out = (newRGBCalFormatTran_out - grayRGB')./grayRGB';

% Transform contrast image
newRGBContrastCalFormatTranContrast_out = newRGBContrastCalFormatTranContrast_out * transformRGBmatrix_opt;

% Add back in gray before outputting the image
triRGBCalFormatTranOpt = (newRGBContrastCalFormatTranContrast_out.*grayRGB') + grayRGB';
triRGBCalFormatOpt     = triRGBCalFormatTranOpt';

% Check for small perterburbations around 0 and 1 in gamut...
% transformations might have pushed it slightly out
if (min(triRGBCalFormatOpt(:)) < 0) && (min(triRGBCalFormatOpt(:)) > -.01)
    triRGBCalFormatOpt(triRGBCalFormatOpt<0) = 0;
end

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
    function loss = loss_function(t_vec, LMSCalFormatTran, M_cones2rgb, lambda, renderType, Disp)
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

        % NORMAL... NON CONTRAST RN
        % LMSContrastCalFormatTran = LMSCalFormatTran;

        % Convert into RGB where gray is removed
        RGBCalFormatTran = LMSCalFormatTran * M_cones2rgb';

        % Create contrast image
        RGBContrastCalFormatTran = (RGBCalFormatTran - grayRGB')./grayRGB';

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

        % Check in gamut?
        minRGB = min(newLMSContrastCalFormatTran(:));
        maxRGB = max(newLMSContrastCalFormatTran(:));

        % Penalize out of gamut
        if minRGB < 0
            lossGamutTerm = 100000000*(minRGB.^2);
        elseif (maxRGB > 1)
            lossGamutTerm = 100000000*((maxRGB-1).^2);
        else
            lossGamutTerm = 0;
        end

        % Variance term
        var_term_raw = varianceLMS("newConeVar",renderType,[],newLMSContrastCalFormat);
        var_term = lambda*var_term_raw;
        totalVariance = whiteNoiseVariance(Disp);
        % Normalize via total variance in white noise image
        var_term_balance = var_term/totalVariance;

        % Similarity term
        similarity_term_raw = similarityLMS('angle',LMSCalFormatTran,newLMSContrastCalFormatTran);
        similarity_term = (1-lambda)*similarity_term_raw;

        % Loss
        % Scale loss so that it is small enough to make fmincon happy but not
        % so small that it is unhappy.
        balanceFactor = 10e5;
        fminconFactor = 10^11/balanceFactor;
        loss = -fminconFactor*(var_term_balance + similarity_term) + lossGamutTerm;

    end


% This ensures that dichromat rendering does not go out of gamut. Calls
% DichromatSimulateBrettel and checks RGB values
    function [c, ceq] = nonlin_con(t_vec, LMSCalFormatTran, M_cones2rgb, cbType, Disp)

        ceq = [];

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

        % Matrix to convert from rgb to xyz
        % Note: this matrix must be applied on the LEFT!!
        M_rgb2xyz = Disp.T_xyz*Disp.P_monitor;
        M_xyz2rgb = inv(M_rgb2xyz);

        % Linear RGB --> XYZ (Brettel takes in XYZ)
        triXYZCalFormat = M_rgb2xyz * newRGBCalFormat;
        % Cal Format --> Image Format
        triXYZImgFormat = CalFormatToImage(triXYZCalFormat,Disp.m,Disp.n);

        % Convert with Brettel
        [diXYZ] = DichromatSimulateBrettel(triXYZImgFormat, cbType, []);

        % Image Format --> CalFormat
        diXYZCalFormat = ImageToCalFormat(diXYZ);
        % XYZ --> Linear RGB
        diRGBLinCalFormat = M_xyz2rgb * diXYZCalFormat;

        % Get c(1) and c(2) from min and max
        c(1) = -1*min(diRGBLinCalFormat(:));
        c(2) = max(diRGBLinCalFormat(:))-1;

    end

% Eventually try and implement something like this to capture better lambda
% value differences
    function [c, ceq] = nonlin_var_constraint(t_vec, LMSCalFormatTran, M_cones2rgb, renderType, v0, v1)
        % NONLIN_VAR_CONSTRAINT ensures that the variance of the transformed data
        % falls within a predefined range of values.
        %
        % Inputs:
        %   t_vec - Flattened 3x3 transformation matrix
        %   LMSCalFormatTran - Original LMS values in calibration format (transposed)
        %   M_cones2rgb - Transformation matrix from LMS to RGB
        %   Disp - Display structure containing transformation matrices
        %   v0, v1 - Lower and upper bounds for variance range
        %
        % Outputs:
        %   c - Inequality constraint (not used)
        %   ceq - Equality constraint enforcing variance within specified range

        % Reshape transformation vector into a 3x3 matrix
        T = reshape(t_vec, 3, 3);

        % Convert LMS to RGB using the transformation matrix
        newRGBCalFormatTran = LMSCalFormatTran * M_cones2rgb' * T;

        % Convert RGB back to LMS
        newLMSCalFormatTran = newRGBCalFormatTran * inv(M_cones2rgb)';

        newLMSCalFormat = newLMSCalFormatTran';
        % Compute variance for each LMS channel
        switch (renderType)
            case 'Deuteranopia' % m cone deficiency
                index = [1 3];
            case 'Protanopia'   % l cone deficiency
                index = [2 3];
            case 'Tritanopia'   % s cone deficiency
                index = [1 2];
        end

        % Variance term
        total_variance = (var(newLMSCalFormat(index(1), :)) + var(newLMSCalFormat(index(2), :)));

        % var_lms = var(newLMSCalFormatTran, 0, 1);
        % total_variance = sum(var_lms);

        % Range of variance vals based on present v0 v1
        variance_range = linspace(v0, v1, 10);

        % Find the closest variance value within the range
        [~, closest_index] = min(abs(variance_range - total_variance));
        closest_variance = variance_range(closest_index);

        % Nonlinear equality constraint: variance must match one in my chosen
        % range
        ceq = total_variance - closest_variance;

        c = [];
    end

end
