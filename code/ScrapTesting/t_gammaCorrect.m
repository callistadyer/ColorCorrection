
%% Gamma correction: how to use RGB2rgbLin and rgbLin2RGB
close all
clear
clc

% Load display parameters
Disp = loadDisplay();

% Get gamma and inverse gamma tables
gammaTable      = displayGet(Disp.d,'gamma table');   % gammaTable used for reverse gamma correction 
invGammaTable   = displayGet(Disp.d,'inverse gamma'); % invGammaTable used for gamma correction 

%% Plot the gamma table and inverse gamma table
figure();
subplot(1,2,1)
plot(gammaTable,LineWidth=2)
axis square
xlabel('Gamma corrected value (input)')
ylabel('Linear rgb value (output)')
title('Gamma Table (DAC to Linear)')

subplot(1,2,2)
plot(invGammaTable,LineWidth=2)
axis square
xlabel('Linear rgb value (input)')
ylabel('Gamma corrected value (output)')
title('Inverse Gamma Table (Linear to DAC)')

%% Let's take a gamma corrected image, then un-correct it, and re-correct it to see what we get

%% Load figure
Disp = loadDisplay();
imgParams = buildSetParameters('ishihara',1,128*3,128*3);
renderType = 'Deuteranopia';
[insideColors, outsideColors] = chooseIshiharaColors(renderType,imgParams.plateType,Disp);
ishiharaRGB = generateIshiharaPlate('74', insideColors, outsideColors,imgParams.m);

%% Gamma corrected image (original)
RGBImage     = im2double(ishiharaRGB); % this transforms from [0 255] to [0 1]
RGBCalFormat = ImageToCalFormat(RGBImage);

%% Load in the gamma tables
gammaTable      = displayGet(Disp.d,'gamma table');
invGammaTable   = displayGet(Disp.d,'inverse gamma');
% invGammaTable is not in the right format, so we normalize to get it into
% [0 1] range.
dacSize = displayGet(Disp.d, 'dacsize');   % should return 10 for 10-bit
maxDAC  = 2^dacSize - 1;                   % 1023 for 10-bit
invGammaTable = invGammaTable / maxDAC;

%% Do the gamma stuff
% RGB gamma corrected -> linear rgb
rgbLinCalFormat = RGB2rgbLin(RGBCalFormat, Disp, imgParams);
rgbLinImage     = CalFormatToImage(rgbLinCalFormat,imgParams.m,imgParams.n);

% linear rgb -> RGB gamma corrected
RGBCalFormat      = rgbLin2RGB(rgbLinCalFormat, Disp, imgParams);
RGBImageRecovered = CalFormatToImage(RGBCalFormat,imgParams.m,imgParams.n);

%% Visualize
% Look at original image, uncorrect gamma, then get back gamma corrected
figure('Position',[474         558        1011         264]);

subplot(1,3,1)
imagesc(RGBImage);
title('original RGB')
axis square

subplot(1,3,2)
imagesc(rgbLinImage);
title('linear rgb')
axis square

subplot(1,3,3)
imagesc(RGBImageRecovered);
title('recovered RGB')
axis square


%% Check if original RGB and recovered RGB are the same within tolerance  
tolerance = .1;  

maxError = abs(max(RGBImage(:) - RGBImageRecovered(:)));
meanA = mean(RGBImage(:));
relativeError = maxError / meanA;

fprintf('Maximum absolute error: %.6f\n', maxError);
fprintf('Relative error: %.6f\n', relativeError);

if relativeError > tolerance
    error('Recovered image differs from original beyond tolerance (%.2f%%)', tolerance*100);
else
    disp('Recovered image matches original within specified tolerance.');
end

