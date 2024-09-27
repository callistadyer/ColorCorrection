function adjustedLMS = correctedLMSadjust(correctedLMS,originalLMS)



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

adjustedLMS = K_optvec' .* correctedLMS;

% K_opt = diag(K_optvec);

function obj = adjustLMSObjectiveFunction(kScale,correctedLMS,originalLMS)

obj = sqrt( sum( ((kScale'.*correctedLMS) - originalLMS ).^2,"all") );


end



























end




