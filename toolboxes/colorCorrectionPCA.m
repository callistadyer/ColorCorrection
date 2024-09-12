function correctedLMS = colorCorrectionPCA(img,originalLMS,renderType,cone_mean_orig)


% Perform PCA
coeff = pca(originalLMS');
PC2D(:,1) = coeff(:, 1); % First principal component
PC2D(:,2) = coeff(:, 2); % Second principal component

% Plot cloud of points visualization
figure();
x = originalLMS(1,:);
y = originalLMS(2,:);
z = originalLMS(3,:);

scatter3(x, y, z, 'filled'); hold on;
xlabel('L')
ylabel('M')
zlabel('S')

% Map PCA onto available cones
% LMS image
T = originalLMS;

% LMS means
T_mean = mean(T,2);

% Mean subtracted LMS values
T_ms = T - T_mean;

% Map mean subtracted LMS valud onto two principle components (linear regression) 
D_ms = PC2D\T_ms;

% Map mean of LMS onto two principle components (linear regression)
D_m = PC2D \ T_mean;

% D_mnew(1,:) = D_ms(1,:);
% D_mnew(2,:) = D_ms(1,:);
% D_mnew(3,:) = D_ms(2,:);

%% Find the scaling matrix that maps D_ms onto approximate cone values

D_mnew(1,:) = D_ms(1,:);
D_mnew(2,:) = D_ms(1,:);
D_mnew(3,:) = D_ms(2,:);

% Get image
[hyperspectralImage wls d P_monitor] = loadImage(img);
% Get cone spectral sensitivities
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);
bScale = 0;
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Initial guess for K
initialK = [1, 1, 1];

% Define the constraints (if any). Here we assume no constraints.
A = [];
b = [];
Aeq = [];
beq = [];
lb = []; % Lower bounds
ub = []; % Upper bounds

% Call fmincon
[K_opt, fval] = fmincon(@(K) T_EstObjectiveFunction(K, D_mnew, T_mean, d, T_cones, P_monitor, m, n, bScale), initialK, A, b, Aeq, beq, lb, ub);

% Display results
disp('Optimal K:');
disp(K_opt);
disp('Objective Function Value:');
disp(fval);

T_opt = K_opt * D_mnew + T_mean;

T_est_rgbImg = LMS2rgbLinimg(T_opt, d, T_cones, P_monitor, m, n, bScale);


%% 

% CHECK: It is at least true that D has two different values
% figure();
% Dimg = CalFormatToImage(D_ms(1,:),256,256);
% imagesc(Dimg); 
% 
% figure();
% Dimg2 = CalFormatToImage(D_ms(2,:),256,256);
% imagesc(Dimg2); 

switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        coneGone = 2;
    case 'Protanopia'   % l cone deficiency
        coneGone = 1;
    case 'Tritanopia'   % s cone deficiency
        coneGone = 3;
end

% Choose M to minimize the difference between Din and Dout
% D_in = T(1:end ~= coneGone,:); % original LMS values, but only 2 rows for the available cones  
% 
% % Get M via linear regression
% M = [D_ms + D_m]' \ D_in';
% 
% % Get estimate of two cones available
% D_est = M' * (D_ms + D_m);

% How to decide which D is paired with which LMS? 
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        LMSimageCalFormat(1,:) = T_opt(1,:);
        LMSimageCalFormat(2,:) = ones(1,length(LMSimageCalFormat(1,:)));
        LMSimageCalFormat(3,:) = T_opt(2,:);
    case 'Protanopia'   % l cone deficiency
        LMSimageCalFormat(2,:) = T_opt(1,:);
        LMSimageCalFormat(1,:) = ones(1,length(LMSimageCalFormat(2,:)));
        LMSimageCalFormat(3,:) = T_opt(2,:);
    case 'Tritanopia'   % s cone deficiency
        LMSimageCalFormat(1,:) = T_opt(1,:);
        LMSimageCalFormat(2,:) = T_opt(2,:);
        LMSimageCalFormat(3,:) = ones(1,length(LMSimageCalFormat(1,:)));
end

% Get values for missing cone
[correctedLMS] = tri2dichromatLMS(LMSimageCalFormat,renderType,cone_mean_orig);

end