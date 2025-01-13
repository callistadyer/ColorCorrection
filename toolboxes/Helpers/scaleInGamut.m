function [correctedLMS_scaled, k] = scaleInGamut(correctedLMS, Disp, bScale)

% Initialize scaling factor k
k = .01;                % Start with no scaling
step_size = 0.0001;      % Increment step size for k
max_iterations = 500; % Prevent infinite loops
check_point = 500;    % Iteration to check for max value comparison

% Initialize variables
iteration = 0;
scaledLMS = correctedLMS; % Initialize scaled LMS
maxValues = zeros(1, max_iterations); % Track max values of scaled LMS

while true
    % Scale correctedLMS by k
    scaledLMS = correctedLMS * k;

    % Check if gamut condition is satisfied
    if checkGamut(scaledLMS, Disp, bScale) ~= 0
        % Break the loop when checkGamut no longer returns 0
        fprintf('Scaling successful at k = %.4f\n', k);
        break;
    end

    % Track the maximum value of scaledLMS
    maxValues(iteration + 1) = max(scaledLMS(:));

    % % Check after 100 iterations if max value at 100th iteration > initial max value
    % if iteration == check_point
    %     if maxValues(check_point) > maxValues(1)
    %         fprintf('Restarting with reversed scaling direction at iteration %d\n', iteration);
    %         % Restart scaling process
    %         k = 1;                % Reset scaling factor
    %         step_size = -step_size; % Reverse scaling direction
    %         iteration = 0;        % Reset iteration counter
    %         maxValues = zeros(1, max_iterations); % Reset max values
    %         continue;
    %     end
    % end

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




% function [correctedLMS_scaled, k] = scaleInGamut(correctedLMS,Disp,bScale)
% 
% % Initialize scaling factor k
% k = 1;                % Start with no scaling
% step_size = 0.1;      % Increment step size for k
% max_iterations = 100; % Prevent infinite loops
% 
% % Apply scaling until checkGamut does not return 0
% iteration = 0;
% while true
%     % Scale correctedLMS by k
%     scaledLMS = correctedLMS * k;
% 
%     % Check if gamut condition is satisfied
%     if checkGamut(scaledLMS,Disp,bScale) ~= 0
%         % Break the loop when checkGamut no longer returns 0
%         fprintf('Scaling successful at k = %.4f\n', k);
%         break;
%     end
% 
%     % Update scaling factor k
%     k = k - step_size;
%     iteration = iteration + 1;
% 
%     % Prevent infinite loops
%     if iteration > max_iterations
%         error('Scaling failed: Maximum iterations reached. Consider increasing step size or checking data.');
%     end
% end
% 
% % Final scaled LMS
% correctedLMS_scaled = scaledLMS;
% 
% % Optionally display the result
% disp('Final scaled correctedLMS:');
% disp(correctedLMS_scaled);
% end