function adjustedLMS = correctedLMSadjust(correctedLMS,originalLMS)

% Second adjustment after PCA to make new/corrected LMS values closer to
% the original LMS values of the input image (minimize the difference
% between the two inputs, correctedLMS and originalLMS)
%
% Syntax:
%   adjustedLMS = correctedLMSadjust(correctedLMS,originalLMS)
%
% Description:
%
% Inputs:
%   correctedLMS:  - LMS after PCA
%   originalLMS:   - LMS before PCA
%
% Outputs:
%   adjustedLMS    - LMS after the difference between the inputs is
%                    minimized
%
% Optional key/value pairs:
%   None
%


% Initial guess for K
initialKscale = [0.1, 0.1, 0.1];

% Define the constraints
A = [];
b = [];
Aeq = [];
beq = [];
lb = []; % Lower bounds
ub = []; % Upper bounds

% Options
 options = optimset('fmincon');
 options = optimset(options,'Diagnostics','off','Display','iter','LargeScale','off','Algorithm','active-set');

% Call fmincon
[K_optvec, fval] = fmincon(@(kScale) adjustLMSObjectiveFunction(kScale,correctedLMS,originalLMS), initialKscale, A, b, Aeq, beq, lb, ub, [], options);

% Scale corrected LMS values by optimal K scalar
adjustedLMS = K_optvec' .* correctedLMS;

function obj = adjustLMSObjectiveFunction(kScale,correctedLMS,originalLMS)
obj = sqrt( sum( ((kScale'.*correctedLMS) - originalLMS ).^2,"all") );
end


end




