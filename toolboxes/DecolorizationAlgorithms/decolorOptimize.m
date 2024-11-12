function [] = decolorOptimize(data,numPCs)
% MATLAB Script For Iterative PCA Using Fmincon Without The PCA Function
% Input: A 3xN Matrix "data"
% Output: Principal Components Stored In Matrix "PCs", Visualized Over Data

if isempty(data)
    % Pretend data to visualize what the function does
    N = 100;                         % Number of points
    x = randn(1, N);                 % Random data for the first dimension
    y = 2 * x + randn(1, N) * 0.2;   % Strong correlation with the first dimension
    z = 0.5 * x + 0.3 * randn(1, N); % Weaker correlation with the first dimension
    data = [x; y; z];
end

% Subtract the mean (center the data)
mean_data = mean(data, 2);
centered_data = data - mean_data; 

% Initialize variables
PCs = zeros(size(data, 1), numPCs); % Matrix to store principal components.

% Options for fmincon
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');

% Loop to extract principal components (loop for each numPCs)
for pcIdx = 1:numPCs
    % Function to maximize variance
    variance_function = @(X) -var(centered_data' * X); % Maximize variance = minimize negative variance.
    
    % Nonlinear constraint function
    constraint_function = @(X) deal([], abs(sum(X.^2)) - 1); % Constraint ||X||^2 = 1.
    
    % Initial guess that satisfies constraint
    X0 = [1; 0; 0]; % Impulse at first position
    
    % Find current principal component via minimizing the negative variance
    [PC, ~] = fmincon(@(X) variance_function(X), X0, [], [], [], [], [], [], constraint_function, options);
    
    % Store PC
    PCs(:, pcIdx) = PC;
    
    % Remove current PC
    projected_data = centered_data' * PC; % Projection onto current PC
    centered_data = centered_data - PC * projected_data'; % Remove variance explained
end

disp('All Principal Components (column-wise):');
disp(PCs);

% Plot the data and PCS
figure;
scatter3(centered_data(1, :), centered_data(2, :), centered_data(3, :), 'filled');
hold on;

% Plot PCs
origin = mean_data; % Use the data mean as the origin for the PCs.
for i = 1:numPCs
    quiver3(origin(1), origin(2), origin(3), ...
        PCs(1, i), PCs(2, i), PCs(3, i), ...
        'LineWidth', 2, 'MaxHeadSize', 0.5); % Draw PC as an arrow.
end

title('3D Data and Principal Components');
xlabel('X1');
ylabel('X2');
zlabel('X3');
grid on;
axis equal;
legend({'Data', 'PC1', 'PC2'}, 'Location', 'Best');
view(3);
rotate3d on;

end
