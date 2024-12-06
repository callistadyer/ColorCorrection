function [PCs, projected_data] = decolorOptimize(data,numPCs,bPLOT)
% Uses fmincon to find the axis of projection that maximizes variance 
%
% Syntax:
%   [PCs projected_data] = decolorOptimize(data,numPCs,bPLOT)
%
% Description:
%
% Inputs:
%   data:   data to be projected
%   numPCs: number of principal components to project onto
%   bPLOT:  1 -> plot
%           0 -> dont
%
% Outputs:
%   PCs:            principal components
%   projected_data: data projected onto PCs
%
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
[PCs projected_data] = decolorOptimize([],2,1)
%}

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
mean_data = mean(data, 2);
centered_data = round(data,4) - round(mean_data,4); 

PCs = zeros(size(data, 1), numPCs); % Matrix to store principal components
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');

theRemainingData = centered_data;
% Loop to get PCs
for pcIdx = 1:numPCs
    % Maximize variance
    variance = @(X) -var(theRemainingData' * X); % Maximize variance = minimize negative variance
    
    % Initial guess that satisfies constraint
    X0 = [1; 0; 0]; 

    % Find current principal component via minimizing the negative variance
    [PC, fval] = fmincon(variance, X0, [], [], [], [], [], [], @constraint_function, options);

    % Store PC
    PCs(:, pcIdx) = PC;
    
    theNullSpace = null((PCs(:,1:pcIdx))');
    
    % Project the data onto the null space.  This in essence 'gets rid'
    % of the part of the data that are explained by the PCA components that
    % we have so far, so that we can then find the direction that has
    % maximum variance when we project this part onto it.
    theRemainingData = (theNullSpace*(theNullSpace'*centered_data));

    projected_data(:,pcIdx)     = centered_data' * PC;                  % Projection onto current PC
end

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


end

function [c, ceq] = constraint_function(X)
    % c    : Inequality constraints (none)
    % ceq  : Equality constraints (||X||^2 - 1 = 0 to ensure unit norm)
    
    c = []; 
    ceq = sum(X.^2) - 1; % Equality constraint ||X||^2 = 1
end
