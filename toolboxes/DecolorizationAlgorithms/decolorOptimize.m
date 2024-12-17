function [PCs, projected_data] = decolorOptimize(data, numPCs, bPLOT, lambda_var, lambda_dot,Disp,bScale)
% Optimizes PCA projections to balance maximizing variance and similarity to original data.
%
% Syntax:
%   [PCs projected_data] = decolorOptimize(data, numPCs, bPLOT, lambda_var, lambda_dot)
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

% Check if lambda values sum to 1
if abs(lambda_var + lambda_dot - 1) > 1e-6
    error('lambda_var + lambda_dot must equal 1');
end

% Sample data to see how code works
if isempty(data)
    % Pretend data to visualize what the function does
    N = 100;                         % Number of points
    x = randn(1, N);                 % Random data for the first dimension
    y = 2 * x + randn(1, N) * 0.2;   % Strong correlation with the first dimension
    z = 0.5 * x + 0.3 * randn(1, N); % Weaker correlation with the first dimension
    data = [x; y; z];
end

% Subtract mean 
mean_data     = mean(data, 2);
centered_data = round(data,4) - round(mean_data,4); 

% Normalize data for initial variance calculation
% Why do this? Magnitude of variance might be very different than dot prod
initial_variance = var(centered_data(:)); % Total variance for normalization
if initial_variance == 0
    initial_variance = eps; % Small positive value to prevent division by zero
end


PCs = zeros(size(data, 1), numPCs); % Matrix to store principal components
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');

theRemainingData = centered_data;
projected_data = zeros(size(centered_data, 2), numPCs);

% Loop to get PCs
for pcIdx = 1:numPCs

    % Objective function: Combined weights of variance and dot product 
    % objective = @(X) - ( ...
    %     lambda_var * var(theRemainingData' * X) / initial_variance + ...
    %     lambda_dot * dot(theRemainingData' * X, theRemainingData' * X) ...
    %     );

    %%%%%%%%%%%%%%% OLD WAY: ONLY MAXIMIZING VARIANCE %%%%%%%%%%%%%%%
    % Maximize variance
    variance = @(X) -var(theRemainingData' * X); % Maximize variance = minimize negative variance
    % 
    % % Initial guess that satisfies constraint
    X0 = [1; 0; 0]; 
    % 
    % % Find current principal component via minimizing the negative variance
    [PC, fval] = fmincon(variance, X0, [], [], [], [], [], [], @constraint_function, options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % MAKING THE ASSUMPTION THAT FIRST CONE MEAN WILL BE ADDED TO FIRST PC
    % objective = @(X) compute_objective_with_penalty(theRemainingData, X, mean_data(numPCs), lambda_var, initial_variance,Disp,bScale);
    
    % Initial guess: Ensure dot product starts near 1????
    % X0 = [1; 0; 0]; 

    % Find current principal component with penalty enforcement
    % [PC, fval] = fmincon(objective, X0, [], [], [], [], [], [], @constraint_function, options);

    % Constraint to ensure ||PC|| = 1 and dot product >= 0
    % constraint = @(X) constraint_function_with_dot(theRemainingData, X);

    % Find current principal component
    % [PC, ~] = fmincon(objective, X0, [], [], [], [], [], [], constraint, options);

    % Store PC
    PCs(:, pcIdx) = PC;
    
    % Define the nullspace of the 
    theNullSpace = null((PCs(:,1:pcIdx))');
    
    % Project the data onto the null space.  This in essence 'gets rid'
    % of the part of the data that are explained by the PCA components that
    % we have so far, so that we can then find the direction that has
    % maximum variance when we project this part onto it.
    theRemainingData = (theNullSpace*(theNullSpace'*centered_data));

    projected_data(:,pcIdx)     = centered_data' * PC;                  % Projection onto current PC
end


%% Find the scaling matrix that maps D_ms onto approximate cone values

% Check with pca function: NOTE, is same but has sign flips
[pcaAuto] = pca(data',"Centered",true);
% not the same as [pcaAuto] = pca(centered_data',"Centered",false); WHY?

% disp('All Principal Components (column-wise):');
% disp(PCs);

% % Uncomment to plot original data
% Plot original data
% figure;
% scatter3(data(1,:),data(2,:),data(3,:),'filled')

if bPLOT == 1
    % Plot the projected data and PCS
    figure;
    scatter3(centered_data(1, :), centered_data(2, :), centered_data(3, :), 'filled');
    hold on;

    % Plot PCs
    origin = mean_data; % Data mean = the origin for the PCs
    for i = 1:numPCs
        quiver3(origin(1), origin(2), origin(3), PCs(1, i), PCs(2, i), PCs(3, i),'LineWidth', 2, 'MaxHeadSize', 0.5);
    end

    title('3D Data and Principal Components');
    xlabel('X1');
    ylabel('X2');
    zlabel('X3');
    axis equal;
    legend({'Data', 'PC1', 'PC2'}, 'Location', 'Best');
    view(3);
    rotate3d on;
end


% function [c, ceq] = constraint_function_with_dot(data, X)
%     % Project data onto X
%     projected_data = data' * X; % Size: m x 1
% 
%     % Compute the dot product between original and projected data
%     c = - dot(data' * X, projected_data); % Enforce dot product >= 0
%     ceq = sum(X.^2) - 1; % Unit norm constraint ||X||^2 = 1
% end

function [c, ceq] = constraint_function(X)
    % c    : Inequality constraints (none)
    % ceq  : Equality constraints (||X||^2 - 1 = 0 to ensure unit norm)

    c = []; 
    ceq = sum(X.^2) - 1; % Equality constraint ||X||^2 = 1
end
end