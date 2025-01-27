function [triLMSCalFormatOpt] = colorCorrectionOptimize(triLMSCalFormat, renderType, lambda_var,Disp,bScale)
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
    T0 = T0(:);

    % OPTIMIZATION SETUP
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter','MaxIterations',50);
    [transformRGB_opt, fval] = fmincon(@(transformRGB) loss_function(transformRGB, triLMSCalFormatTran, M_cones2rgb, lambda_var, renderType, Disp), ...
        T0, A, b, [], [], [], [], [], options);

    fValOpt = loss_function(transformRGB_opt, triLMSCalFormatTran, M_cones2rgb, lambda_var,renderType, Disp);
    bCheck = A*transformRGB_opt(:);
    if (any(bCheck > b))
        fprintf('Failed to satisfy constraint\n');
    end

    % RESHAPE THE OPTIMAL SOLUTION INTO MATRIX FORM
    transformRGBmatrix_opt = reshape(transformRGB_opt, 3, 3);

    % Calculate optimal output from input * optimal transform
    triRGBCalFormatTranOpt = triLMSCalFormatTran * M_cones2rgb' * transformRGBmatrix_opt;
    triRGBCalFormatOpt     = triRGBCalFormatTranOpt';

    % Check if given O is in gamut
    % Get LMS values
    triLMSCalFormatOpt = M_rgb2cones * triRGBCalFormatOpt;
    % Check if is in gamut
    inGamutAfterTransform = checkGamut(triLMSCalFormatOpt,Disp,bScale);
    if inGamutAfterTransform == 0
        error(['ERROR: decolorOptimize, constraint is not working, transformation is pushing out of gamut'])
    end


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

    newRGBCalFormatTran = LMSCalFormatTran * M_cones2rgb' * T;

    % Convert to LMS
    newLMSCalFormatTran = newRGBCalFormatTran*inv(M_cones2rgb)';

    % Get into cal format
    newLMSCalFormat = newLMSCalFormatTran';
    % newRGBCalFormat = newRGBCalFormatTran';

    % Check in gamut?
    minRGB = min(newRGBCalFormatTran(:));
    maxRGB = max(newRGBCalFormatTran(:));

    % Penalize out of gamut
    if minRGB < 0
        lossGamutTerm = 100000000*(minRGB.^2);
    elseif (maxRGB > 1)
        lossGamutTerm = 100000000*((maxRGB-1).^2);
    else
        lossGamutTerm = 0;
    end

    switch (renderType)
        case 'Deuteranopia' % m cone deficiency
            index = [1 3];
        case 'Protanopia'   % l cone deficiency
            index = [2 3];
        case 'Tritanopia'   % s cone deficiency
            index = [1 2];
    end

    % Variance term
    var_term_raw = (var(newLMSCalFormat(index(1), :)) + var(newLMSCalFormat(index(2), :)));
    var_term = lambda*var_term_raw;
    
    % Set a balance factor that brings the variance term to order 1
    balanceFactor = 10e4;
    var_term_balance = balanceFactor * var_term;

    % Similarity term
    similarityType = 'angle';
    switch (similarityType)
        case 'angle'
            similarity_term_raw = (newLMSCalFormatTran(:)'*LMSCalFormatTran(:))/(norm(newLMSCalFormatTran(:))*norm(LMSCalFormatTran(:)));
        case 'distance'
            similarity_term_raw = norm(newLMSCalFormatTran(:)-LMSCalFormatTran(:))/norm(LMSCalFormatTran(:));
        otherwise
            error('Unknown similarity type specified');
    end
    similarity_term = (1-lambda)*similarity_term_raw;
    
    % Loss
    %
    % Scale loss so that it is small enough to make fmincon happy but not
    % so small that it is unhappy.
    fminconFactor = 10^11/balanceFactor;
    loss = -fminconFactor*(var_term_balance + similarity_term) + lossGamutTerm;

    % Scale loss into a reasonable range
    % loss = 10e11*loss;
    
    end

end
