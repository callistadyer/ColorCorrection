
function [triRGBImgFormatCorrected] = dichromatCorrection(img,renderType,bScale,method,nSquares,modType,lambda_var)
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
close all;

% Get trichromatic (LMS) image
[triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(img,renderType,0,bScale,nSquares,modType);    

%%% Currently trying block processing - Callista, evenetually come back to
%%% this to make it a toggle (ie block or not)

% % 8x8 = 64
% % 4x4 = 16 blocks
% blockNum  = [4, 4];
% % Block size [m x n] ... currently assumes divisible by blockNum
% blockSize = [Disp.m/blockNum(1) Disp.n/blockNum(2)];
% 
% % Get image into image format for block processing
% triLMSImgFormat       = CalFormatToImage(triLMSCalFormat,Disp.m,Disp.n);
% triLMSImgFormat_plate = CalFormatToImage(triLMSCalFormat_plate,Disp.m,Disp.n);
% 
% % Function that calls color correction algorithm (this gets inserted into
% % blockproc)
% fun = @(blockStruct) processBlock(blockStruct.data, method, renderType, lambda_var, Disp, bScale);
% 
% % Apply block processing for both the original image and the one with the
% % isochromatic plate
% triLMSImgCorrected = blockproc(triLMSImgFormat, blockSize, fun);
% triLMSImgFormatCorrected_plate = blockproc(triLMSImgFormat_plate, blockSize, fun);
% 
% % Convert back to cal format
% triLMScalFormatCorrected = ImageToCalFormat(double(triLMSImgCorrected));
% triLMScalFormatCorrected_plate = ImageToCalFormat(double(triLMSImgFormatCorrected_plate));


% This is also happening in colorCorrect, which is called by the block
% processing function
switch (method)
    case 'linTransform'
        % decolorOptimize does mean subtraction, then maximizes variance fmincon 
        % expects x y z dimensions in rows and measurements in columns ie. [3 x 1000]  
        % lambda_var = 0.1;
        [triLMScalFormatCorrected] = colorCorrectionOptimize(triLMSCalFormat,renderType,lambda_var,Disp);
        [triLMScalFormatCorrected_plate] = colorCorrectionOptimize(triLMSCalFormat_plate,renderType,lambda_var,Disp);
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
diLMSCalFormatCorrected                  = tri2dichromatLMSCalFormat(triLMScalFormatCorrected,renderType,Disp,bScale); 
diRGBCalFormatCorrected                  = LMS2RGBCalFormat(diLMSCalFormatCorrected, Disp,bScale);
diLMSCalFormatCorrected_plate            = tri2dichromatLMSCalFormat(triLMScalFormatCorrected_plate,renderType,Disp,bScale);      % isochromatic plate 
diRGBCalFormatCorrected_plate            = LMS2RGBCalFormat(diLMSCalFormatCorrected_plate, Disp,bScale); % isochromatic plate 


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


figure('Position',[214         261         501        1076]);

% Create a tiled layout
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add images to the tiles
nexttile
imshow(triRGBImgFormatOrig);
title('trichromat');

nexttile
imshow(diRGBImgFormatOrig);
title('dichromat');

nexttile
imshow(triRGBImgFormatCorrected);
title('trichromat corrected');

nexttile
imshow(diRGBImgFormatCorrected);
title('dichromat corrected');

nexttile
imshow(triRGBImgFormatOrig_plate);
title('trichromat - plate');

nexttile
imshow(diRGBImgFormatOrig_plate);
title('dichromat - plate');

nexttile
imshow(triRGBImgFormatCorrected_plate);
title('trichromat corrected - plate');

nexttile
imshow(diRGBImgFormatCorrected_plate);
title('dichromat corrected - plate');

sgtitle('Left: Trichromat, Right: Dichromat')
end


% Process block by block. Gotta use an extra function because blockproc
% seems to want only one input at a time
% function [triLMSImgFormatCorrected] = processBlock(block, method, renderType, lambda_var, Disp, bScale)
% 
%     triLMSCalFormat = ImageToCalFormat(double(block)); % Convert to double and cal format
% 
%     % Size of current block
%     newSize = size(block);
%     blockDisp = Disp;
%     blockDisp.m = newSize(1);
%     blockDisp.n = newSize(2);
% 
%     % Call appropriate color correction function
%     [triLMSCalFormatCorrected] = colorCorrect(method, triLMSCalFormat, renderType, lambda_var, blockDisp, bScale);
% 
%     % Get image in image format
%     triLMSImgFormatCorrected = CalFormatToImage(triLMSCalFormatCorrected,newSize(1),newSize(2));
% 
% 
% end
