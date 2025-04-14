
function [triRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T, T_P] = dichromatCorrection(lambdaOrVar,var,lambda_var,img,renderType,method,nSquares,modType,constraintWL,T_prev,T_prev_P)
% Transform trichromatic image so that dichromat can see more color
% contrast. Also want to try and preserve some naturalness. This is
% accomplished in colorCorrectionOptimize where we incorporate similarity
% to original in the loss function
%
% Syntax:
%   [triRGBImgFormatCorrected] = dichromatCorrection(img,renderType,method,nSquares,modType,lambda_var)
%
% Description:
%
% Inputs:
%   lambdaOrVar:  - String. Use lambda range or specific variance (computed from lamdbas)  
%                   Optimize using lambda value (between 0 and 1) or var
%                   which samples linspace between the variances of
%                   lambda = 0 and lambda = 1
%                       'lambda'
%                       'var'
%   img:          - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
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
%   lambda_var    - Weight on the variance term in the optimization [0 1]
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
%   T_prev_P      - Initial transformation matrix for the RGB image(isochromatic plate version)  
%                   Same thing as above, but you have to get another start
%                   point for the isochromatic plate version of the image.
%{
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

if strcmp(img,'ishihara')

    plateType = 3;
    [insideColors, outsideColors] = chooseIshiharaColors(renderType,plateType,Disp);

    % Generate plate now that you have the correct colors
    ishiharaRGB = generateIshiharaPlate('74', insideColors, outsideColors,Disp.m);
    ishiharaRGB = im2double(ishiharaRGB);

    % Plot modified RGB Image 
    figure();imagesc(ishiharaRGB)
    axis square;

    % Put modified image into LMS 
    triRGBCalFormat    = ImageToCalFormat(ishiharaRGB);
    triLMSCalFormat       = Disp.M_rgb2cones * triRGBCalFormat;
    triLMSCalFormat_plate = triLMSCalFormat;

    % Run the modulated image through the linear dichromat simulation
    [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, grayLMS,  constraintWL, renderType, Disp);
    [diLMSCalFormat_plate,M_triToDi] = DichromSimulateLinear(triLMSCalFormat_plate, grayLMS,  constraintWL, renderType, Disp);

    % Check
    % rgb1 = inv(Disp.M_rgb2cones) * diLMSCalFormat;
    % image = CalFormatToImage(rgb1,Disp.m,Disp.n);
    % figure();imagesc(image); axis square;
    
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)

    img_rgb = im2double(imread(img));           % Load and convert image to double
    if Disp.m > 128 || Disp.n > 128
        scaleFactor = 0.3;                      % Downsample
        img_rgb = imresize(img_rgb, scaleFactor);
        [rows, cols, ~] = size(img_rgb);
        Disp.m         = cols;
        Disp.n         = rows;
    end

    % Cal format RGB 
    triRGBCalFormat = ImageToCalFormat(img_rgb); 
    % Convert to LMS 
    triLMSCalFormat = Disp.M_rgb2cones * triRGBCalFormat;

    [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, Disp.grayLMS,  constraintWL, renderType, Disp);
    
    % check
    rgb1 = inv(Disp.M_rgb2cones) * diLMSCalFormat;
    imageDi = CalFormatToImage(rgb1, Disp.m, Disp.n);
    figure();imagesc(imageDi);

else
    % Get trichromatic (LMS) image
    [triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(img,renderType,0,nSquares,modType);
    clear triLMSCalFormat;
    clear diLMSCalFormat;
    triLMSCalFormat = triLMSCalFormat_plate; % do this when you just want to see the isochromatic plate square version (other is just gray)
    diLMSCalFormat  = diLMSCalFormat_plate;
end
%%%%%%%%%%%%%%%%%% PUT ALL OF THIS JUNK IN ANOTHER FUNCTION ^^^^^^ 


% This is also happening in colorCorrect, which is called by the block
% processing function
switch (method)
    case 'linTransform'
        % decolorOptimize does mean subtraction, then maximizes variance fmincon
        % expects x y z dimensions in rows and measurements in columns ie. [3 x 1000]
        disp('Entering optimization function');
        [triLMScalFormatCorrected,s_raw, v_raw, s_bal, v_bal, T] = colorCorrectionOptimize(lambdaOrVar,var,lambda_var,triLMSCalFormat,renderType,constraintWL,T_prev,Disp);
            % triLMScalFormatCorrected_plate = triLMScalFormatCorrected;
            s_raw_P = s_raw;
            v_raw_P = v_raw;
            s_bal_P = s_bal;
            v_bal_P = v_bal;
            T_P = T;
    case 'easyPCA'
        triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,renderType,Disp);
        % triLMScalFormatCorrected_plate = colorCorrectionEasyPCA(triLMSCalFormat_plate,renderType,Disp);
    case 'hardPCA'
        numPCs = 2;
        triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat,numPCs,Disp);
        % triLMScalFormatCorrected_plate = colorCorrectionHardPCA(triLMSCalFormat_plate,numPCs,Disp);
end
% Create RGB image from LMS
%%%%% ORIGINAL IMAGE AND PLATE %%%%%
% Dichromat simulation of original image
diRGBCalFormatOrig = M_cones2rgb * diLMSCalFormat;
% [diRGBCalFormatOrig]        = LMS2RGBCalFormat(diLMSCalFormat, Disp);

% Trichromat simulation of original image
triRGBcalFormatOrig = Disp.M_cones2rgb * triLMSCalFormat;
% [triRGBcalFormatOrig]       = LMS2RGBCalFormat(triLMSCalFormat, Disp);

%%%%% CORRECTED IMAGE AND PLATE %%%%%
% Corrected trichromat image
triRGBcalFormatCorrected = M_cones2rgb * triLMScalFormatCorrected;
% [triRGBcalFormatCorrected]        = LMS2RGBCalFormat(triLMScalFormatCorrected, Disp);


[diLMSCalFormatCorrected,~]        = DichromSimulateLinear(triLMScalFormatCorrected, Disp.grayLMS,  constraintWL, renderType, Disp);
diRGBCalFormatCorrected            = Disp.M_cones2rgb * diLMSCalFormatCorrected;
% diRGBCalFormatCorrected            = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp);

% diLMSCalFormatCorrected                  = tri2dichromatLMSCalFormat(triLMScalFormatCorrected,renderType,Disp);
% diRGBCalFormatCorrected                  = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp);

disp('callista!!!!! Need to gamma correct!!!!');

% Transform from cal format to image for viewing
% original trichromat
triRGBImgFormatOrig              = CalFormatToImage(triRGBcalFormatOrig,Disp.m,Disp.n); % no modulation

% corrected trichromat
triRGBImgFormatCorrected         = CalFormatToImage(triRGBcalFormatCorrected,Disp.m,Disp.n); % no modulation

% original dichromat
diRGBImgFormatOrig               = CalFormatToImage(diRGBCalFormatOrig,Disp.m,Disp.n); % no modulation

% corrected dichromat
diRGBImgFormatCorrected          = CalFormatToImage(diRGBCalFormatCorrected,Disp.m,Disp.n); % no modulation


figure('Position',[161   302   562   552]);

% Create a tiled layout
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add images to the tiles
% nexttile
% imshow(triRGBImgFormatOrig);
% title('trichromat');
%
% nexttile
% imshow(diRGBImgFormatOrig);
% title('dichromat');
%
% nexttile
% imshow(triRGBImgFormatCorrected);
% title('trichromat corrected');
%
% nexttile
% imshow(diRGBImgFormatCorrected);
% title('dichromat corrected');

% nexttile
% imshow(triRGBImgFormatOrig_plate);
% title('trichromat - plate');
%
% nexttile
% imshow(diRGBImgFormatOrig_plate);
% title('dichromat - plate');
nexttile
imshow(triRGBImgFormatOrig);
title('trichromat original image');

nexttile(diRGBImgFormatOrig);
title('dichromat original image');

nexttile
imshow(triRGBImgFormatCorrected);
title('trichromat corrected');

nexttile
imshow(diRGBImgFormatCorrected);
title('dichromat corrected');

% sgtitle(['lambdavar = ' num2str(lambda_var)])
sgtitle(['var = ' num2str(var)])

end

