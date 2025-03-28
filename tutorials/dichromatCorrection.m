
function [triRGBImgFormatCorrected,s_raw_P, v_raw_P, s_bal_P, v_bal_P] = dichromatCorrection(img,renderType,bScale,method,nSquares,modType,lambda_var,var)
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
lambdas = linspace(0,.1,10)
for i = 1:length(lambdas)
[RGBImage_dichromat] = dichromatCorrection('gray','Deuteranopia',0,'linTransform',1,'M',lambdas(i));
end

[RGBImage_dichromat] = dichromatCorrection('gray','Deuteranopia',0,'linTransform',1,'M',0.1);
[RGBImage_dichromat] = dichromatCorrection('74','Deuteranopia',0,'linTransform',1,'M',0.1);
[RGBImage_dichromat] = dichromatCorrection('scene2.mat','Deuteranopia',1,'linTransform',10,'M',0.1);
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

    % insideColors = [
    %     1.0 0.6 0.3;
    %     0.9 0.7 0.2;
    %     1.0 0.5 0.5
    %     ];
    % outsideColors = [
    %     1.0 0.6 0.3;
    %     0.9 0.7 0.2;
    %     1.0 0.5 0.5
    %     ];
    insideColors = [
        0.5 0.5 0.5;
        0.5 0.5 0.5;
        0.5 0.5 0.5
        ];
    outsideColors = [
        0.5 0.5 0.5;
        0.5 0.5 0.5;
        0.5 0.5 0.5
        ];

    imgRGB = generateIshiharaPlate('74', insideColors, outsideColors,imgSize);
    imgRGB = im2double(imgRGB);
    figure();imagesc(imgRGB)
    axis square;

    rgbImageCalFormat = ImageToCalFormat(imgRGB);

    switch (modType)
        case 'rand'
            % % This code case allows the squares to be random colors
            % modulationDirection_LMS = modulationDirection_LMS;
        case 'M' % m cone deficiency
            modulationDirection_LMS = [0 1 0]';
        case 'L'   % l cone deficiency
            modulationDirection_LMS = [1 0 0]';
        case 'S'   % s cone deficiency
            modulationDirection_LMS = [0 0 1]';
    end
    modulationDirection_rgb = M_cones2rgb*modulationDirection_LMS;


    % NOTE: this scaleFactor_rgb is for scaling the modulation direction. This
    % is different from scaleFactor that goes into rgb->LMS conversions (and
    % vice versa) which scales the image values to a sensible number
    % Background is the value of each pixel: this determines cone contrast separately for each pixel
    for i = 1:size(insideColors,1)
        scaleFactor_rgb(i) = MaximizeGamutContrast(modulationDirection_rgb,insideColors(i,:)'); % bg is in rgb cal format
        % Stay away from the very edge
        toleranceFactor = 0.8;
        % Scale modulation direction by scale factor to get modulation=
        modulation_rgb(i,:) = scaleFactor_rgb(i).*toleranceFactor.*modulationDirection_rgb;
    end

    % M_rgb2cones * modulation_rgb(1,:)'
    insideColorsMod = insideColors + modulation_rgb;

    imgRGBmod = generateIshiharaPlate('74', insideColorsMod, outsideColors,imgSize);
    imgRGBmod = im2double(imgRGBmod);
    figure();imagesc(imgRGBmod)
    axis square;
    imgRGBmodCalFormat = ImageToCalFormat(imgRGBmod);
    triLMSCalFormat = M_rgb2cones * imgRGBmodCalFormat;
    triLMSCalFormat_plate = triLMSCalFormat;
    constraintWl=585;

    [diLMSCalFormat,M_triToDi] = DichromSimulateLinear(triLMSCalFormat, grayLMS,  constraintWl, renderType, Disp);
    [diLMSCalFormat_plate,M_triToDi] = DichromSimulateLinear(triLMSCalFormat_plate, grayLMS,  constraintWl, renderType, Disp);

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
        [triLMScalFormatCorrected,s_raw, v_raw, s_bal, v_bal] = colorCorrectionOptimize(var, triLMSCalFormat,renderType,lambda_var,Disp);
        [triLMScalFormatCorrected_plate,s_raw_P, v_raw_P, s_bal_P, v_bal_P] = colorCorrectionOptimize(var, triLMSCalFormat_plate,renderType,lambda_var,Disp);
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

[diLMSCalFormatCorrected,~]        = DichromSimulateLinear(triLMScalFormatCorrected, grayLMS,  585, renderType, Disp);
diRGBCalFormatCorrected            = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp,bScale);
[diLMSCalFormatCorrected_plate,~]  = DichromSimulateLinear(triLMScalFormatCorrected_plate, grayLMS,  585, renderType, Disp);
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

