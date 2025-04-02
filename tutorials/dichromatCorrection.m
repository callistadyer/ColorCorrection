
function [triRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T, T_P] = dichromatCorrection(img,renderType,bScale,method,nSquares,modType,lambda_var,constraintWL,T_prev,T_prev_P)
% Transform trichromatic image so that dichromat can see more color
% contrast. Also want to try and preserve some naturalness. This is
% accomplished in colorCorrectionOptimize where we incorporate similarity
% to original in the loss function
%
% Syntax:
%   [triRGBImgFormatCorrected] = dichromatCorrection(img,renderType,bScale,method,nSquares,modType,lambda_var)
%
% Description:
%
% Inputs:
%   img:          - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   bScale:       - Boolean. Scale the image values into display range (1
%                   or 0).  A good idea except for 'gray'.
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
%{
                      for i = 1:30
                          lambda = linspace(0,1,30)
                          T{1} = eye(3,3);
                          T_P{1} = eye(3,3);
                          [RGBImage_dichromat,s_raw_P(i), v_raw_P(i), s_bal_P(i), v_bal_P(i),T{i+1},T_P{i+1}] = dichromatCorrection('gray','Deuteranopia',0,'linTransform',1,'M',lambda(i),585,T{i},T_P{i},1);
                      end
%}
%   T_prev_P      - Initial transformation matrix for the RGB image(isochromatic plate version)  
%                   Same thing as above, but you have to get another start
%                   point for the isochromatic plate version of the image.
%
% Outputs:
%   triRGBImgFormatCorrected:  - Transformed RGB image after PCA and scaling. Also replaced missing cone as done in other code
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
% Loop over lambdas:
for i = 1:30
lambda = linspace(0,1,30)
lambdaval = lambda(i)
T{1} = eye(3,3);
T_P{1} = eye(3,3);
[RGBImage_dichromat,s_raw_P(i), v_raw_P(i), s_bal_P(i), v_bal_P(i),T{i+1},T_P{i+1}] = dichromatCorrection('gray','Deuteranopia',0,'linTransform',1,'M',lambda(i),585,T{i},T_P{i});
end

[RGBImage_dichromat,s_raw_P, v_raw_P, s_bal_P, v_bal_P] = dichromatCorrection('gray','Deuteranopia',0,'linTransform',1,'M',0,585,eye(3,3),eye(3,3));
[RGBImage_dichromat,s_raw_P, v_raw_P, s_bal_P, v_bal_P] = dichromatCorrection('ishihara','Deuteranopia',0,'linTransform',1,'M',0,585,eye(3,3),eye(3,3));

%}

% Close out any stray figures
% close all;

if strcmp(img,'ishihara')
    % Display
    % Wavelengths for display
    imgSize = 128;
    wls = (400:10:700)';
    d = displayCreate('LCD-Apple');
    P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);
    load T_cones_ss2;
    T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);


    Disp.wls = wls;
    Disp.d = d;
    Disp.P_monitor = P_monitor;
    Disp.m = imgSize;
    Disp.n = imgSize;
    Disp.T_cones = T_cones;
    M_rgb2cones = Disp.T_cones*Disp.P_monitor;
    M_cones2rgb = inv(M_rgb2cones);
    grayRGB = [0.5 0.5 0.5]';
    grayLMS = M_rgb2cones*grayRGB;

    % Make initial inside (number) and outside (background) colors the
    % same... then add to the inside colors by only adding the M cone color
%     insideColors = [
%     0.95, 0.6, 0.4;   % light orange
%     0.85, 0.4, 0.3;   % muted red-orange
%     0.9,  0.5, 0.2    % pumpkin tone
% ];
    insideColors = [
    0.35  0.3   0.3;  
    0.3   0.45  0.4;  
    0.55  0.25  0.6  
    ];
    insideColors = [
    0.5  0.5   0.5;  
    0.5   0.5  0.5;  
    0.5  0.5  0.5  
    ];

    outsideColors = insideColors;

    % % See if it's working... 
    % imgRGB = generateIshiharaPlate('74', insideColors, outsideColors,imgSize);
    % imgRGB = im2double(imgRGB);
    % figure();imagesc(imgRGB)
    % axis square;

    switch (modType)
        case 'rand'
            % % This code case allows the squares to be random colors
            % modulationDirection_LMS = modulationDirection_LMS;
        case 'M' % m cone deficiency
            modulationDirection_LMS  = [0 1 0]';
            modulationDirection_back = [1 .5 0]';
        case 'L'   % l cone deficiency
            modulationDirection_LMS = [1 0 0]';
        case 'S'   % s cone deficiency
            modulationDirection_LMS = [0 0 1]';
    end

    % Direction of color (e.g., move in M cone direction to make a plate
    % that deuteranopes cannot see
    modulationDirection_rgb      = M_cones2rgb*modulationDirection_LMS;
    modulationDirection_rgb_back = M_cones2rgb*modulationDirection_back;

    % NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
    % is different from scaleFactor that goes into rgb->LMS conversions (and
    % vice versa) which scales the image values to a sensible number
    % Background is the value of each pixel: this determines cone contrast separately for each pixel
    for i = 1:size(insideColors,1)
        scaleFactor_rgb(i)      = MaximizeGamutContrast(modulationDirection_rgb,insideColors(i,:)'); % bg is in rgb cal format
        scaleFactor_rgb_back(i) = MaximizeGamutContrast(modulationDirection_rgb_back,outsideColors(i,:)'); % bg is in rgb cal format

        % Stay away from the very edge
        toleranceFactor = 0.7;
        % Scale modulation direction by scale factor to get modulation=
        modulation_rgb(i,:)      = scaleFactor_rgb(i).*modulationDirection_rgb;
        modulation_rgb_back(i,:) = scaleFactor_rgb_back(i).*toleranceFactor.*modulationDirection_rgb_back;
    end

    % New colors. Outside colors stay the same. Inside colors simply add
    % missing cone information to existing background colors.
    % insideColorsMod  = insideColors + modulation_rgb.*[.9 .7 .5]';
    insideColorsMod  = insideColors + modulation_rgb;

    outsideColorsMod = outsideColors;

    % Generate plate now that you have the correct colors
    imgRGBmod = generateIshiharaPlate('74', insideColorsMod, outsideColorsMod,imgSize);
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

else
    % Get trichromatic (LMS) image
    [triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(img,renderType,0,bScale,nSquares,modType);
end

% This is also happening in colorCorrect, which is called by the block
% processing function
switch (method)
    case 'linTransform'
        % decolorOptimize does mean subtraction, then maximizes variance fmincon
        % expects x y z dimensions in rows and measurements in columns ie. [3 x 1000]
        [triLMScalFormatCorrected,s_raw, v_raw, s_bal, v_bal, T] = colorCorrectionOptimize(triLMSCalFormat,renderType,lambda_var,constraintWL,T_prev,Disp);
        [triLMScalFormatCorrected_plate,s_raw_P, v_raw_P, s_bal_P, v_bal_P, T_P] = colorCorrectionOptimize(triLMSCalFormat_plate,renderType,lambda_var,constraintWL,T_prev_P, Disp);
    case 'easyPCA'
        triLMScalFormatCorrected = colorCorrectionEasyPCA(triLMSCalFormat,renderType,Disp,bScale);
        triLMScalFormatCorrected_plate = colorCorrectionEasyPCA(triLMSCalFormat_plate,renderType,Disp,bScale);
    case 'hardPCA'
        numPCs = 2;
        triLMScalFormatCorrected = colorCorrectionHardPCA(triLMSCalFormat,numPCs,Disp);
        triLMScalFormatCorrected_plate = colorCorrectionHardPCA(triLMSCalFormat_plate,numPCs,Disp);
end


%%%% Old colorCorrection function %%%%
% [triLMScalFormatCorrected       triLMSmeans]                   = colorCorrection(triLMSCalFormat,renderType,Disp,bScale);   % Original image
% [triLMScalFormatCorrected_plate triLMSmeans_plate]             = colorCorrection(triLMSCalFormat_plate,renderType,Disp,bScale); % Image with plate
% correctedLMS = K_opt_plate * D_mnew + T_mean_plate;

% Create RGB image from LMS
%%%%% ORIGINAL IMAGE AND PLATE %%%%%
% Dichromat simulation of original image
[diRGBCalFormatOrig]        = LMS2RGBCalFormat(diLMSCalFormat, Disp,bScale);
[diRGBCalFormatOrig_plate]  = LMS2RGBCalFormat(diLMSCalFormat_plate, Disp,bScale);         % isochromatic plate

% Trichromat simulation of original image
[triRGBcalFormatOrig]       = LMS2RGBCalFormat(triLMSCalFormat, Disp,bScale);
[triRGBcalFormatOrig_plate] = LMS2RGBCalFormat(triLMSCalFormat_plate, Disp,bScale);  % isochromatic plate

%%%%% CORRECTED IMAGE AND PLATE %%%%%
% Corrected trichromat image
[triRGBcalFormatCorrected]        = LMS2RGBCalFormat(triLMScalFormatCorrected, Disp,bScale);
[triRGBcalFormatCorrected_plate]  = LMS2RGBCalFormat(triLMScalFormatCorrected_plate, Disp,bScale);       % isochromatic plate

% Corrected dichromat image
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);
grayRGB = [0.5 0.5 0.5]';
grayLMS = M_rgb2cones*grayRGB;

[diLMSCalFormatCorrected,~]        = DichromSimulateLinear(triLMScalFormatCorrected, grayLMS,  constraintWL, renderType, Disp);
diRGBCalFormatCorrected            = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp,bScale);
[diLMSCalFormatCorrected_plate,~]  = DichromSimulateLinear(triLMScalFormatCorrected_plate, grayLMS,  constraintWL, renderType, Disp);
diRGBCalFormatCorrected_plate      = LMS2RGBCalFormat(diLMSCalFormatCorrected_plate, Disp,bScale); % isochromatic plate

% diLMSCalFormatCorrected                  = tri2dichromatLMSCalFormat(triLMScalFormatCorrected,renderType,Disp,bScale);
% diRGBCalFormatCorrected                  = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp,bScale);
% diLMSCalFormatCorrected_plate            = tri2dichromatLMSCalFormat(triLMScalFormatCorrected_plate,renderType,Disp,bScale);      % isochromatic plate
% diRGBCalFormatCorrected_plate            = LMS2RGBCalFormat(diLMSCalFormatCorrected_plate, Disp,bScale); % isochromatic plate


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

sgtitle(['lambdavar = ' num2str(lambda_var)])
end

