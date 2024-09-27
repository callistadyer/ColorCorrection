function [correctedLMS, K_opt, D_mnew, T_mean]  = colorCorrectionPCA(img,originalLMS,renderType,cone_mean_orig,bScale)


% Perform PCA
coeff = pca(originalLMS');
PC2D(:,1) = coeff(:, 1); % First principal component
PC2D(:,2) = coeff(:, 2); % Second principal component

% Plot cloud of points visualization
% figure();
% x = originalLMS(1,:);
% y = originalLMS(2,:);
% z = originalLMS(3,:);
% 
% scatter3(x, y, z, 'filled'); hold on;
% xlabel('L')
% ylabel('M')
% zlabel('S')

% Map PCA onto available cones
% LMS image
T = originalLMS;

% LMS means
T_mean = mean(T,2);

% Mean subtracted LMS values
% T_ms = T - T_mean;
T_ms = T; % Looks better when we don't mean subtract here?

% Map mean subtracted LMS valud onto two principle components (linear regression) 
D_ms = PC2D\T_ms;

% Map mean of LMS onto two principle components (linear regression)
% D_m = PC2D \ T_mean;

%% Find the scaling matrix that maps D_ms onto approximate cone values

D_mnew(1,:) = D_ms(1,:);
D_mnew(2,:) = D_ms(1,:);
D_mnew(3,:) = D_ms(2,:);

% Get image
[hyperspectralImage wls d P_monitor] = loadImage(img);
% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Initial guess for K
initialKvec = [0.1, 0.1, 0.1];

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
[K_optvec, fval] = fmincon(@(kVec) T_EstObjectiveFunction(kVec, D_mnew, T_mean, d, T_cones, P_monitor, m, n, bScale), initialKvec, A, b, Aeq, beq, lb, ub, [], options);

K_opt = diag(K_optvec);
% Display results
disp('Optimal K:');
disp(K_optvec);
disp('Objective Function Value:');
disp(fval);

T_opt = K_opt * D_mnew;
% T_opt = K_opt * D_mnew + T_mean;

T_est_rgbImg = LMS2rgbLinimg(T_opt, d, T_cones, P_monitor, m, n, bScale);

% Get values for missing cone
correctedLMS = T_opt;

end