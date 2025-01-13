function [triLMSCalFormatOpt] = decolorOptimize(triLMSCalFormat,method, numPCs, renderType, bPLOT, lambda_var,Disp,bScale)
% Optimizes PCA projections to balance maximizing variance and similarity to original data.
%
% Syntax:
%   [PCs projected_data] = decolorOptimize(data, numPCs, bPLOT, lambda_var, lambda_dot,Disp,bScale)
%
% Inputs:
%   data:       Data to be projected (dimensions: variables x samples)
%   numPCs:     Number of principal components to project onto
%   bPLOT:      1 -> Plot the result, 0 -> Don't plot
%   lambda_var: Weight for maximizing variance (0 <= lambda_var <= 1)
%   lambda_dot: Weight for maximizing similarity (0 <= lambda_dot <= 1)
%
% Outputs:
%   PCs:            Principal components (column vectors)
%   projected_data: Data projected onto PCs
%
% Constraints:
%   - ||PC||^2 = 1 (unit norm)
%   - Dot product must be >= 0
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
[PCs projected_data] = decolorOptimize([],2,1, 0.5, 0.5)
[PCs projected_data] = decolorOptimize([],2,1)
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
if strcmp(method,"linTransform")

    % NOTE: I = LMS values [nPix x 3]
    triLMSCalFormatTran = triLMSCalFormat'; % data = [3 x nPix]
    % M = transformation matrix [3x3]'
    % NOTE: MUST APPLY THIS M MATRIX ON THE LEFT OF CAL FORMAT DATA
    M_rgb2cones = Disp.T_cones*Disp.P_monitor;
    M_cones2rgb = inv(M_rgb2cones);

    % Weight for variance term
    lambda = 0.9;             
    % Number of pixels
    nPix   = size(triLMSCalFormat,2);

    % Constraint matrix (A, includes lots of I iterations) and vector (b)
    % triLMSCalFormatTran: trichromat LMS values in [nPix x 3] form
    triRGBCalFormatTranOpt = triLMSCalFormat' * M_cones2rgb';
    A_upper = blkdiag(triRGBCalFormatTranOpt, triRGBCalFormatTranOpt, triRGBCalFormatTranOpt);      % Upper constraint blocks
    A_lower = -A_upper;              % for -I * X <= 0
    A = double([A_upper; A_lower]); 
    b = [ones(nPix * 3, 1); zeros(nPix * 3, 1)];

    % Initial guess at transformation matrix - start with identity
    T0 = eye(3, 3);
    T0 = T0(:);

    % OPTIMIZATION SETUP
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter','MaxIterations',30);
    [transformRGB_opt, fval] = fmincon(@(transformRGB) loss_function(transformRGB, triLMSCalFormatTran, M_cones2rgb, lambda, renderType, Disp), ...
        T0, A, b, [], [], [], [], [], options);

    fValOpt = loss_function(transformRGB_opt, triLMSCalFormatTran, M_cones2rgb', lambda,renderType, Disp)
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
end


%% PCA the hard way
if strcmp(method,"hardPCA")
    disp(['Note: "pca the hard way"... problem because variance is 0 at starting point'])

    % Subtract mean
    mean_data     = mean(triLMSCalFormat, 2);
    centered_data = round(triLMSCalFormat,4) - round(mean_data,4);

    % Normalize data for initial variance calculation
    % Why do this? Magnitude of variance might be very different than dot prod
    initial_variance = var(centered_data(:)); % Total variance for normalization
    if initial_variance == 0
        initial_variance = eps; % Small positive value to prevent division by zero
    end

    PCs = zeros(size(triLMSCalFormat, 1), numPCs); % Matrix to store principal components
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');

    theRemainingData = centered_data;
    projected_dataLMS = zeros(size(centered_data, 2), numPCs);

    % Loop to get PCs
    for pcIdx = 1:numPCs

        %%%%%%%%%%%%%% OLD WAY: ONLY MAXIMIZING VARIANCE %%%%%%%%%%%%%%%
        % Maximize variance
        variance = @(X) -var(theRemainingData' * X); % Maximize variance = minimize negative variance
        %
        % % Initial guess that satisfies constraint
        X0 = [1; 0; 0];
        %
        % % Find current principal component via minimizing the negative variance
        [PC, fval] = fmincon(variance, X0, [], [], [], [], [], [], @constraint_function, options);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Objective function: Combined weights of variance and dot product
        % objective = @(X) - ( ...
        %     lambda_var * var(theRemainingData' * X) / initial_variance + ...
        %     (1-lambda_var) * dot(theRemainingData' * X, theRemainingData' * X) ...
        %     );


        % Store PC
        PCs(:, pcIdx) = PC;

        % Define the nullspace of the
        theNullSpace = null((PCs(:,1:pcIdx))');

        % Project the data onto the null space.  This in essence 'gets rid'
        % of the part of the data that are explained by the PCA components that
        % we have so far, so that we can then find the direction that has
        % maximum variance when we project this part onto it.
        theRemainingData = (theNullSpace*(theNullSpace'*centered_data));

        projected_dataLMS(:,pcIdx)     = centered_data' * PC;                  % Projection onto current PC
    end
end


%% Plot? 

% Check with pca function: NOTE, is same but has sign flips
% [pcaAuto] = pca(data',"Centered",true);
% not the same as [pcaAuto] = pca(centered_data',"Centered",false); WHY?

% % Uncomment to plot original data
% Plot original data
% figure;
% scatter3(data(1,:),data(2,:),data(3,:),'filled')

% if bPLOT == 1
%     mean_data     = mean(data, 2);
%     centered_data = round(data,4) - round(mean_data,4);
%     % Plot the projected data and PCS
%     figure;
%     scatter3(centered_data(1, :), centered_data(2, :), centered_data(3, :), 'filled');
%     hold on;
% 
%     % Plot PCs
%     origin = mean_data; % Data mean = the origin for the PCs
%     for i = 1:numPCs
%         quiver3(origin(1), origin(2), origin(3), PCs(1, i), PCs(2, i), PCs(3, i),'LineWidth', 2, 'MaxHeadSize', 0.5);
%     end
% 
%     title('3D Data and Principal Components');
%     xlabel('X1');
%     ylabel('X2');
%     zlabel('X3');
%     axis equal;
%     legend({'Data', 'PC1', 'PC2'}, 'Location', 'Best');
%     view(3);
%     rotate3d on;
% end

%% Functions 

% OBJECTIVE FUNCTION
    function loss = loss_function(t_vec, LMSCalFormatTran, M_cones2rgb, lambda, renderType, Disp)
    T = reshape(t_vec, 3, 3);       % RESHAPE x_vec INTO 3x3 MATRIX
    % X = M_cones2rgb * T 
    % O = I * X
    % O = I * X'*inv(M');                     % CALCULATE O = X * I
   
    % I - LMS
    % O - linear RGB
    % M - LMS2RGB
    % T - RGB TRANSFORMATION
    % [nPizels x 3]     = [3 x nPixels] x [3 x 3] x [3 x 3]
    RGBCalFormatTran = LMSCalFormatTran * M_cones2rgb' * T;

    % How to implement constraints?


    % Check if given O is in gamut
    % Convert to LMS
    LMSCalFormatTran = RGBCalFormatTran*inv(M_cones2rgb');

    % Get into cal format
    LMSCalFormat = LMSCalFormatTran';
    RGBCalFormat = RGBCalFormatTran';

    minRGB = min(RGBCalFormatTran(:));
    maxRGB = max(RGBCalFormatTran(:));

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
    var_term = lambda * (var(LMSCalFormat(index(1), :)) + var(LMSCalFormat(index(2), :)));

    % Similarity term
    dot_term = (1 - lambda) *   (dot(LMSCalFormatTran(:, 1), LMSCalFormatTran(:, 1)) / norm(LMSCalFormatTran(:, 1))) + ...
                                (dot(LMSCalFormatTran(:, 2), LMSCalFormatTran(:, 2)) / norm(LMSCalFormatTran(:, 2))) + ...
                                (dot(LMSCalFormatTran(:, 3), LMSCalFormatTran(:, 3)) / norm(LMSCalFormatTran(:, 3)));
    % Loss
    loss = -(var_term + dot_term) + lossGamutTerm;
    
    end

end

% LINEAR INEQUALITY MATRIX CONSTRAINT FUNCTION 
% function [c, ceq] = constraint_function(x, A, b)
%     x_vec = reshape(x', 9, 1);
%     c = A * x_vec - b;    % INEQUALITY CONSTRAINT: A * x_vec <= b
%     ceq = [];             % NO EQUALITY CONSTRAINTS
% end

% % MAXIMIZING VARIANCE CONSTRAINT FUNCTION
function [c, ceq] = constraint_function(X)
    % c    : Inequality constraints (none)
    % ceq  : Equality constraints (||X||^2 - 1 = 0 to ensure unit norm)

    c = []; 
    ceq = sum(X.^2) - 1; % Equality constraint ||X||^2 = 1
end
