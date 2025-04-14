
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

    [insideColors, outsideColors] = chooseIshiharaColors(renderType,plateType,Disp);
    % plateType: choose which kind you wanna make!
    %       1  -> gray background with only missing cone modulation in number
    %       2  -> some random colors with only missing cone modulation in number
    %       3  -> outside colors are just cone modulations of present
    %             cones, inside colors are just cone modulation of missing cones

    % Make initial inside (number) and outside (background) colors the
    % same... then add to the inside colors by only adding the M cone color
    % insideColors = [
    %     0.85    0.45    0.30;
    %     0.90    0.65    0.50;
    %     0.75    0.40    0.35;
    %     ];
    % LS_directions = [
    %     .2  0  0;
    %     .2  0 .2;
    %     0  0  .2;
    %     ]';

% % Normalize each direction to unit length in LMS contrast space
% LS_directions = LS_directions ./ vecnorm(LS_directions);
% 
% rgbColors = zeros(3, 3);
% 
% for i = 1:3
%     modulation_LMS = LS_directions(:,i);
% 
%     % convert to RGB
%     modulation_RGB = M_cones2rgb * modulation_LMS;
% 
%     % Maximize the amount of that RGB direction we can add to gray
%     scaleFactor = MaximizeGamutContrast(modulation_RGB, grayRGB);
%     scaleFactor = scaleFactor/3; % dont go all the way to the edge of gamut    
%     % Final color: background + scaled RGB direction
%     rgbModulated = grayRGB + scaleFactor * modulation_RGB;
% 
%     % Store
%     rgbColors(i, :) = rgbModulated';
% end
% 
% insideColors  = rgbColors;
%     insideColors = [
%     0.85    0.45    0.30;
%     0.90    0.65    0.50;
%     0.75    0.40    0.35;
%     ];
% outsideColors = insideColors;
% 
%     % % See if it's working... 
%     % imgRGB = generateIshiharaPlate('74', insideColors, outsideColors,imgSize);
%     % imgRGB = im2double(imgRGB);
%     % figure();imagesc(imgRGB)
%     % axis square;
% 
%     switch (modType)
%         case 'rand'
%             % % This code case allows the squares to be random colors
%             % modulationDirection_LMS = modulationDirection_LMS;
%         case 'M' % m cone deficiency
%             modulationDirection_LMS = [0 1 0]';
%         case 'L'   % l cone deficiency
%             modulationDirection_LMS = [1 0 0]';
%         case 'S'   % s cone deficiency
%             modulationDirection_LMS = [0 0 1]';
%     end
% 
%     % Direction of color (e.g., move in M cone direction to make a plate
%     % that deuteranopes cannot see
%     modulationDirection_rgb      = Disp.M_cones2rgb*modulationDirection_LMS;
% 
%     % NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
%     % is different from scaleFactor that goes into rgb->LMS conversions (and
%     % vice versa) which scales the image values to a sensible number
%     % Background is the value of each pixel: this determines cone contrast separately for each pixel
%     for i = 1:size(insideColors,1)
%         scaleFactor_rgb(i)      = MaximizeGamutContrast(modulationDirection_rgb,insideColors(i,:)'); % bg is in rgb cal format
%         % Stay away from the very edge
%         % toleranceFactor = 0.9;
%         % Scale modulation direction by scale factor to get modulation=
%         modulation_rgb(i,:)      = scaleFactor_rgb(i).*modulationDirection_rgb;
%     end
% 
%     % New colors. Outside colors stay the same. Inside colors simply add
%     % missing cone information to existing background colors.
%     % insideColorsMod  = insideColors + modulation_rgb.*[.9 .7 .5]';
%     insideColorsMod  = insideColors + modulation_rgb;

    % Generate plate now that you have the correct colors
    imgRGBmod = generateIshiharaPlate('74', insideColors, outsideColors,imgSize);
    imgRGBmod = im2double(imgRGBmod);

    % Plot modified RGB Image 
    figure();imagesc(imgRGBmod)
    axis square;

    % Put modified image into LMS 
    imgRGBmodCalFormat    = ImageToCalFormat(imgRGBmod);
    triLMSCalFormat       = M_rgb2cones * imgRGBmodCalFormat;
    triLMSCalFormat_plate = triLMSCalFormat;

    % Run the modulated image through the linear dichromat simulation
    [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, grayLMS,  constraintWL, renderType, Disp);
    [diLMSCalFormat_plate,M_triToDi] = DichromSimulateLinear(triLMSCalFormat_plate, grayLMS,  constraintWL, renderType, Disp);

    % check
    rgb1 = inv(Disp.M_rgb2cones) * diLMSCalFormat;
    image = CalFormatToImage(rgb1,Disp.m,Disp.n);
    figure();imagesc(image); axis square;
    
elseif strcmp(img,'pic')

    % img_rgb = im2double(imread('slide.png'));       % Replace with your filename
    % [rows, cols, ~] = size(img_rgb);                 % Get correct dimensions


    img_rgb = im2double(imread('slide.png'));    % Load and convert image to double
    scaleFactor = 0.3;                           % Downsample to 50% of original size
    img_rgb = imresize(img_rgb, scaleFactor);   % Resize the image
    [rows, cols, ~] = size(img_rgb);            % Get new dimensions

    wls = (400:10:700)';
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);
    load T_cones_ss2;
    T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);

    Disp.wls = wls;
    Disp.d = d;
    Disp.P_monitor = P_monitor;
    Disp.m = cols;
    Disp.n = rows;
    Disp.T_cones = T_cones;
    M_rgb2cones = Disp.T_cones*Disp.P_monitor;
    M_cones2rgb = inv(M_rgb2cones);
    grayRGB = [0.5 0.5 0.5]';
    grayLMS = M_rgb2cones*grayRGB;

    calFormatRGB = ImageToCalFormat(img_rgb);  
    triLMSCalFormat = M_rgb2cones * calFormatRGB;
    triLMSCalFormat_plate = triLMSCalFormat;

    [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, grayLMS,  constraintWL, renderType, Disp);
    
    % check
    rgb1 = inv(M_rgb2cones) * diLMSCalFormat;
    imageDi = CalFormatToImage(rgb1, Disp.m, Disp.n);

    figure();imagesc(imageDi);

    diLMSCalFormat_plate = diLMSCalFormat;
else
    % Get trichromatic (LMS) image
    [triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(img,renderType,0,nSquares,modType);
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
        if strcmp(img,'ishihara') || strcmp(img,'pic')
            triLMScalFormatCorrected_plate = triLMScalFormatCorrected;
            s_raw_P = s_raw;
            v_raw_P = v_raw;
            s_bal_P = s_bal;
            v_bal_P = v_bal;
            T_P = T;
        else
        [triLMScalFormatCorrected_plate,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T_P] = colorCorrectionOptimize(lambdaOrVar,var,lambda_var,triLMSCalFormat_plate,renderType,constraintWL,T_prev_P, Disp);
        end
    case 'easyPCA'
        triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,renderType,Disp);
        triLMScalFormatCorrected_plate = colorCorrectionEasyPCA(triLMSCalFormat_plate,renderType,Disp);
    case 'hardPCA'
        numPCs = 2;
        triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat,numPCs,Disp);
        triLMScalFormatCorrected_plate = colorCorrectionHardPCA(triLMSCalFormat_plate,numPCs,Disp);
end

%%%% Old colorCorrection function %%%%
% [triLMScalFormatCorrected       triLMSmeans]                   = colorCorrection(triLMSCalFormat,renderType,Disp);   % Original image
% [triLMScalFormatCorrected_plate triLMSmeans_plate]             = colorCorrection(triLMSCalFormat_plate,renderType,Disp); % Image with plate
% correctedLMS = K_opt_plate * D_mnew + T_mean_plate;

% Create RGB image from LMS
%%%%% ORIGINAL IMAGE AND PLATE %%%%%
% Dichromat simulation of original image
diRGBCalFormatOrig = M_cones2rgb * diLMSCalFormat;
diRGBCalFormatOrig_plate = M_cones2rgb * diLMSCalFormat_plate;

% [diRGBCalFormatOrig]        = LMS2RGBCalFormat(diLMSCalFormat, Disp);
% [diRGBCalFormatOrig_plate]  = LMS2RGBCalFormat(diLMSCalFormat_plate, Disp);         % isochromatic plate

% Trichromat simulation of original image
triRGBcalFormatOrig = M_cones2rgb * triLMSCalFormat;
triRGBcalFormatOrig_plate = M_cones2rgb * triLMSCalFormat_plate;
% 
% [triRGBcalFormatOrig]       = LMS2RGBCalFormat(triLMSCalFormat, Disp);
% [triRGBcalFormatOrig_plate] = LMS2RGBCalFormat(triLMSCalFormat_plate, Disp);  % isochromatic plate

%%%%% CORRECTED IMAGE AND PLATE %%%%%
% Corrected trichromat image
triRGBcalFormatCorrected = M_cones2rgb * triLMScalFormatCorrected;
triRGBcalFormatCorrected_plate = M_cones2rgb * triLMScalFormatCorrected_plate;
% 
% [triRGBcalFormatCorrected]        = LMS2RGBCalFormat(triLMScalFormatCorrected, Disp);
% [triRGBcalFormatCorrected_plate]  = LMS2RGBCalFormat(triLMScalFormatCorrected_plate, Disp);       % isochromatic plate

% Corrected dichromat image
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);
grayRGB = [0.5 0.5 0.5]';
grayLMS = M_rgb2cones*grayRGB;

[diLMSCalFormatCorrected,~]        = DichromSimulateLinear(triLMScalFormatCorrected, grayLMS,  constraintWL, renderType, Disp);
diRGBCalFormatCorrected            = M_cones2rgb * diLMSCalFormatCorrected;
% diRGBCalFormatCorrected            = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp);

[diLMSCalFormatCorrected_plate,~]  = DichromSimulateLinear(triLMScalFormatCorrected_plate, grayLMS,  constraintWL, renderType, Disp);
diRGBCalFormatCorrected_plate            = M_cones2rgb * diLMSCalFormatCorrected_plate;
% diRGBCalFormatCorrected_plate      = LMS2RGBCalFormat(diLMSCalFormatCorrected_plate, Disp); % isochromatic plate

% diLMSCalFormatCorrected                  = tri2dichromatLMSCalFormat(triLMScalFormatCorrected,renderType,Disp);
% diRGBCalFormatCorrected                  = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp);
% diLMSCalFormatCorrected_plate            = tri2dichromatLMSCalFormat(triLMScalFormatCorrected_plate,renderType,Disp);      % isochromatic plate
% diRGBCalFormatCorrected_plate            = LMS2RGBCalFormat(diLMSCalFormatCorrected_plate, Disp); % isochromatic plate


disp('callista!!!!! Need to gamma correct!!!!');

% Transform from cal format to image for viewing
% original trichromat
triRGBImgFormatOrig              = CalFormatToImage(triRGBcalFormatOrig,Disp.m,Disp.n); % no modulation
triRGBImgFormatOrig_plate        = CalFormatToImage(triRGBcalFormatOrig_plate,Disp.m,Disp.n); % modulation

% corrected trichromat
triRGBImgFormatCorrected         = CalFormatToImage(triRGBcalFormatCorrected,Disp.m,Disp.n); % no modulation
triRGBImgFormatCorrected_plate   = CalFormatToImage(triRGBcalFormatCorrected_plate,Disp.m,Disp.n); % modulation

% original dichromat
diRGBImgFormatOrig               = CalFormatToImage(diRGBCalFormatOrig,Disp.m,Disp.n); % no modulation
diRGBImgFormatOrig_plate         = CalFormatToImage(diRGBCalFormatOrig_plate,Disp.m,Disp.n); % modulation

% corrected dichromat
diRGBImgFormatCorrected          = CalFormatToImage(diRGBCalFormatCorrected,Disp.m,Disp.n); % no modulation
diRGBImgFormatCorrected_plate    = CalFormatToImage(diRGBCalFormatCorrected_plate,Disp.m,Disp.n); % plate


figure('Position',[161   302   562   552]);

% Create a tiled layout
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

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
imshow(triRGBImgFormatCorrected_plate);
title('trichromat corrected - plate');

nexttile
imshow(diRGBImgFormatCorrected_plate);
title('dichromat corrected - plate');

% sgtitle(['lambdavar = ' num2str(lambda_var)])
sgtitle(['var = ' num2str(var)])

end

