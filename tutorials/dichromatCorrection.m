
function [triRGBImgFormatCorrected,diRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T] = dichromatCorrection(lambdaOrVar,var,lambda_var,img,renderType,varianceType,similarityType,plateType,method,nSquares,modType,constraintWL,T_prev,V0,V1)
% Transform trichromatic image so that dichromat can see more color
% contrast. Also want to try and preserve some naturalness. This is
% accomplished in colorCorrectionOptimize where we incorporate similarity
% to original in the loss function
%
% Syntax:
%   [triRGBImgFormatCorrected,diRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T] = dichromatCorrection(lambdaOrVar,var,lambda_var,img,renderType,varianceType,similarityType,plateType,method,nSquares,modType,constraintWL,T_prev,V0,V1)
%
% Description:
%
% Inputs:
% 
% T_prev,V0,V1)
%   lambdaOrVar:  - String. Use lambda range or specific variance (computed from lamdbas)  
%                   Optimize using lambda value (between 0 and 1) or var
%                   which samples linspace between the variances of
%                   lambda = 0 and lambda = 1
%                       'lambda'
%                       'var'
%   var:          - Double. Variance value. Leave empty if you are using
%                   lambda. Otherwise, we get the var value by
%                   interpolating between the variances at lambda=0 and
%                   lambda=1. Why are we doing this? We don't get nice
%                   evenly spaced image transformations when we use
%                   linspace(lambda=0, lambda=1, 10). Instead, we search
%                   for 10 images where we linspace between the variances
%                   at those endpoints:
%                   linspace(var(lambda=0),var(lambda=1),10) where var is
%                   taken from v_raw_P (or maybe v_bal_P?? not sure)
%   lambda_var:   - Weight on the variance term in the optimization [0 1]
%   img:          - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   varianceType: - String. Method for quantifying the variance or contrast
%                   enhancement term in optimization
%                       'LMdifferenceContrast' ** this one is pretty good 
%                       'regress' 
%                       'delta'                
%                       'detail'
%                       'newConeVar' ** original
%
%   similarityType: - FILL IN
%   plateType:    - FILL IN
%   method:       - Color correction method:
%                       'linTransform'
%                       'easyPCA'
%                       'hardPCA'
%   nSquares:     - number of squares in isochromatic plate
%
%   modType       - type of isochromatic plate modulation
%                       'rand'
%                       'M'
%                       'L'
%                       'S'
%   constraintWL  - Wavelength that forms plane with gray that the
%                   dichromat image gets projected onto in
%                   DichromSimulateLinear.m 
%                        585 for deuteronopes
%   T_prev        - Initial transformation matrix for the RGB image 
%                   Start this at T_prev = eye(3,3). Usually this is most
%                   useful when you are looping over lambdas (see below)
%                   because you want to use the previous T solution to
%                   initialize the next optimization. This avoids some
%                   wonky failures in the fmincon routine. 
%   V0:           - FILL IN
%   V1:           - FILL IN
%
%{

% EXAMPLE INPUTS:

lambdaOrVar = 'lambda';
var = [];
lambda_var = 0.5;
img = 'gray';
renderType = 'Deuteranopia';
varianceType = 'LMdifferenceContrast';
similarityType = 'squared';
plateType = []; % only for ishihara plates
method = 'linTransform';
nSquares = 1;
modType = 'M';
constraintWL = 585;
T_prev = eye(3,3); % change this to be last run
V0 = [];
V1 = []; % These are only relavent when youre using var instead of lambda

                    for i = 1:10
                        T{1} = eye(3,3);
                        T_P{1} = eye(3,3);
                        [RGBImage_dichromat,s_raw_P(i), v_raw_P(i), s_bal_P(i), v_bal_P(i),T{i+1},T_P{i+1}] = dichromatCorrection('var',i,lambda,'ishihara','Deuteranopia','linTransform',1,'M',585,T{i},T_P{i});
                    end
                      % See how variance and similarity changes:
                      figure();
                      subplot(2,2,1); plot(lambda,s_raw_P,'-o',linewidth=2,markersize=10,markerfacecolor='white'); xlabel('lambda',fontsize=20); ylabel('similarity',fontsize=20); title('raw',fontsize=25);
                      subplot(2,2,2); plot(lambda,s_bal_P,'-o',linewidth=2,markersize=10,markerfacecolor='white'); xlabel('lambda',fontsize=20); ylabel('similarity',fontsize=20); title('balanced',fontsize=25);
                      subplot(2,2,3); plot(lambda,v_raw_P,'-o',linewidth=2,markersize=10,markerfacecolor='white'); xlabel('lambda',fontsize=20); ylabel('variance',fontsize=20); title('raw',fontsize=25);
                      subplot(2,2,4); plot(lambda,v_bal_P,'-o',linewidth=2,markersize=10,markerfacecolor='white'); xlabel('lambda',fontsize=20); ylabel('variance',fontsize=20); title('balanced',fontsize=25);
                      figure();
                      subplot(1,2,1);
                      plot(lambda,s_raw_P+v_raw_P,'-o',linewidth=2,markersize=10,markerfacecolor='white'); xlabel('lambda',fontsize=20); ylabel('similarity+variance',fontsize=20); title('raw',fontsize=25);
                      subplot(1,2,2);
                      plot(lambda,s_bal_P+v_bal_P,'-o',linewidth=2,markersize=10,markerfacecolor='white'); xlabel('lambda',fontsize=20); ylabel('similarity+variance',fontsize=20); title('balanced',fontsize=25);

%}
%
% Outputs:
%   triRGBImgFormatCorrected:  - Transformed RGB image after PCA and scaling. Also replaced missing cone as done in other code
%   s_raw_P                    - raw similarity values for current lambda (_P indicates that it is for the modulated image)  
%   v_raw_P                    - raw variance values for current lambda
%   s_bal_P                    - balanced similarity values for current lambda = (1-lambda) * s_raw_P   
%   v_bal_P                    - balanced variance values for current lambda = (lambda) * v_raw_P 
%   T
%   T_P
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

% History
%   09/05/2024  cmd  Initial go.
%
% Examples:
%{

%%%%%% MANIPULATING LAMBDA %%%%%%
% TEST: this one has a lambda of 0 so should be no transformation = original image
[RGBImage_dichromat,s_raw_P, v_raw_P, s_bal_P, v_bal_P] = dichromatCorrection('lambda',[],0,'gray','Deuteranopia','linTransform',1,'M',585,eye(3,3),eye(3,3));
% TEST: this one has a lambda of 1 so should be max transformation = black and gray or white and gray image (max change in gamut) 
[RGBImage_dichromat,s_raw_P, v_raw_P, s_bal_P, v_bal_P] = dichromatCorrection('lambda',[],1,'gray','Deuteranopia','linTransform',1,'M',585,eye(3,3),eye(3,3));
% Loop over lambdas:
for i = 1:30
    lambda = linspace(0,1,30)
    T{1} = eye(3,3);
    T_P{1} = eye(3,3);
    [RGBImage_dichromat,s_raw_P(i), v_raw_P(i), s_bal_P(i), v_bal_P(i),T{i+1},T_P{i+1}] = dichromatCorrection('lambda',[],lambda(i),'gray','Deuteranopia','linTransform',1,'M',585,T{i},T_P{i});
                                                                                          dichromatCorrection('lambda',[],lambda(i),'gray','Deuteranopia','LMdifferenceContrast','squared',plateType,'linTransform',nSquares,modType,585,T_prev,V0,V1)
end

%%%%%% MANIPULATING VAR %%%%%% 
you need to have already run the lambda of 0 and 1,
then used linspace to sample between the variances at 0 and 1. Then this
"var" setting essentially finds you solutions with the variances equally
spaced between the variance at lambda = 0 and lambda = 1

[RGBImage_dichromat,s_raw_P, v_raw_P, s_bal_P, v_bal_P] = dichromatCorrection('var',10,[],'ishihara','Deuteranopia','linTransform',1,'M',585,eye(3,3),eye(3,3));
% Loop over vars:
for i = 1:10
    T{1} = eye(3,3);
    T_P{1} = eye(3,3);
    [RGBImage_dichromat,s_raw_P(i), v_raw_P(i), s_bal_P(i), v_bal_P(i),T{i+1},T_P{i+1}] = dichromatCorrection('var',i,[],'ishihara','Deuteranopia','linTransform',1,'M',585,T{i},T_P{i});
end

%}



% Load display 
Disp = loadDisplay(img);

% Load LMS values for this image
[triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues(img,renderType,modType,nSquares,constraintWL,plateType,Disp);

% Color Correction Algorithm
switch (method)
    case 'linTransform'
        % decolorOptimize does mean subtraction, then maximizes variance fmincon
        % expects x y z dimensions in rows and measurements in columns ie. [3 x 1000]
        disp('Entering optimization function');
        [triLMScalFormatCorrected,s_raw, v_raw, s_bal, v_bal, T] = colorCorrectionOptimize(lambdaOrVar,var,lambda_var,triLMSCalFormat,renderType,varianceType,similarityType,constraintWL,T_prev,Disp,V0,V1);
            % triLMScalFormatCorrected_plate = triLMScalFormatCorrected;
            s_raw_P = s_raw;
            v_raw_P = v_raw;
            s_bal_P = s_bal;
            v_bal_P = v_bal;
    case 'easyPCA'
        triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,renderType,Disp);
        % triLMScalFormatCorrected_plate = colorCorrectionEasyPCA(triLMSCalFormat_plate,renderType,Disp);
    case 'hardPCA'
        numPCs = 2;
        triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat,numPCs,Disp);
        % triLMScalFormatCorrected_plate = colorCorrectionHardPCA(triLMSCalFormat_plate,numPCs,Disp);
end


% Imaging the transformation 
disp('callista!!!!! Need to gamma correct!!!!');

%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%
% Create RGB image from LMS
% Dichromat simulation of original image
diRGBCalFormatOrig = Disp.M_cones2rgb * diLMSCalFormat;
[diRGBCalFormatOrig]        = LMS2RGBCalFormat(diLMSCalFormat, Disp);

% Trichromat simulation of original image
triRGBcalFormatOrig = Disp.M_cones2rgb * triLMSCalFormat;
[triRGBcalFormatOrig]       = LMS2RGBCalFormat(triLMSCalFormat, Disp);

%%%%%%%%%%%%%%% CORRECTED %%%%%%%%%%%%%%%
% Corrected trichromat image
% triRGBcalFormatCorrected = Disp.M_cones2rgb * triLMScalFormatCorrected;
[triRGBcalFormatCorrected]        = LMS2RGBCalFormat(triLMScalFormatCorrected, Disp);

[diLMSCalFormatCorrected,~]        = DichromSimulateLinear(triLMScalFormatCorrected, Disp.grayLMS,  constraintWL, renderType, Disp);
% diRGBCalFormatCorrected            = Disp.M_cones2rgb * diLMSCalFormatCorrected;
diRGBCalFormatCorrected            = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp);


% Transform to RGB image format for viewing:
% original trichromat
triRGBImgFormatOrig              = CalFormatToImage(triRGBcalFormatOrig,Disp.m,Disp.n); % no modulation
% corrected trichromat
triRGBImgFormatCorrected         = CalFormatToImage(triRGBcalFormatCorrected,Disp.m,Disp.n); % no modulation
% original dichromat
diRGBImgFormatOrig               = CalFormatToImage(diRGBCalFormatOrig,Disp.m,Disp.n); % no modulation
% corrected dichromat
diRGBImgFormatCorrected          = CalFormatToImage(diRGBCalFormatCorrected,Disp.m,Disp.n); % no modulation


figure('Position',[161   302   562   552]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
imshow(triRGBImgFormatOrig);
title('trichromat original image');

nexttile
imshow(triRGBImgFormatCorrected);
title('trichromat corrected');

nexttile
imshow(diRGBImgFormatOrig);
title('dichromat original image');

nexttile
imshow(diRGBImgFormatCorrected);
title('dichromat corrected');


if strcmp(lambdaOrVar,'var')
sgtitle(['var = ' num2str(var) ', variance: ' varianceType])
elseif strcmp(lambdaOrVar,'lambda')
sgtitle(['lambdavar = ' num2str(lambda_var) ', variance: ' varianceType])
end


figure(); tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); nexttile
imshow(triRGBImgFormatCorrected);
title(['trichromat corrected, var = ' num2str(var)]);
nexttile
imshow(diRGBImgFormatCorrected);
title('dichromat corrected');


end

