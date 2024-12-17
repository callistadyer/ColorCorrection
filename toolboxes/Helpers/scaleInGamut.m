function [correctedLMS_scaled, k] = scaleInGamut(correctedLMS,Disp,bScale)

% Initialize scaling factor k
k = 1;                % Start with no scaling
step_size = 0.1;      % Increment step size for k
max_iterations = 100; % Prevent infinite loops

% Apply scaling until checkGamut does not return 0
iteration = 0;
while true
    % Scale correctedLMS by k
    scaledLMS = correctedLMS * k;
    
    % Check if gamut condition is satisfied
    if checkGamut(scaledLMS,Disp,bScale) ~= 0
        % Break the loop when checkGamut no longer returns 0
        fprintf('Scaling successful at k = %.4f\n', k);
        break;
    end
    
    % Update scaling factor k
    k = k - step_size;
    iteration = iteration + 1;

    % Prevent infinite loops
    if iteration > max_iterations
        error('Scaling failed: Maximum iterations reached. Consider increasing step size or checking data.');
    end
end

% Final scaled LMS
correctedLMS_scaled = scaledLMS;

% Optionally display the result
disp('Final scaled correctedLMS:');
disp(correctedLMS_scaled);
end