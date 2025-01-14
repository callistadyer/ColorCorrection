function triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat,numPCs)
% Perform color correction using hard PCA on trichromatic LMS data
%
% Syntax:
%   triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat, numPCs)
%
% Description:
%   This function applies principal component analysis (PCA) to the input 
%   trichromatic LMS data using a "hard PCA" method. It iteratively computes 
%   principal components by maximizing variance while projecting the data 
%   onto the null space of previously computed components.
%
% Inputs:
%   triLMSCalFormat:   - Matrix (NxM). Input trichromatic LMS data with N dimensions 
%                        and M samples.
%   numPCs:            - Integer. Number of principal components to compute.
%
% Outputs:
%   triLMScalFormatCorrected: - Matrix (3xM). Corrected trichromatic LMS data 
%                               transformed by the computed PCA components.
%
% Optional key/value pairs:
%   None
%
% Notes:
%   - Uses a custom optimization routine (`fmincon`) to compute the PCA 
%     components.
%   - Ensures that the computed principal components are unit-norm.
%   - Removes components explained by previously computed principal components 
%     through null space projection.
%   - If the initial variance of the data is zero, a small positive value (`eps`) 
%     is used to avoid division by zero.
%
% Examples:
%{
triLMSCalFormat = rand(3, 100); % Example trichromatic LMS data
numPCs = 2; % Number of principal components
triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat, numPCs);
%}

% History:
%   01/14/2025  cmd creates func that was previously embedded in
%   colorCorrectionPCA
%

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

    triLMScalFormatCorrected(1,:) = projected_dataLMS(:,1);
    triLMScalFormatCorrected(2,:) = projected_dataLMS(:,1);
    triLMScalFormatCorrected(3,:) = projected_dataLMS(:,2);

function [c, ceq] = constraint_function(X)
    % c    : Inequality constraints (none)
    % ceq  : Equality constraints (||X||^2 - 1 = 0 to ensure unit norm)

    c = []; 
    ceq = sum(X.^2) - 1; % Equality constraint ||X||^2 = 1
end
end