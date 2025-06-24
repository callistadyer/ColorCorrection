function [triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues(img,renderType,modType,nSquares,constraintWL,plateType,Disp)
% loadLMSvalues  Loads or generates an image and converts to LMS for dichromat simulation
%
% Syntax:
%   [triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues(img,renderType,modType,nSquares,constraintWL,plateType,Disp)
%
% Inputs:
%   img:              Either 'ishihara' or a filename ('.png', '.jpg') or a hyperspectral identifier
%   renderType:       Type of dichromacy to simulate
%                         'Deuteranopia', 'Protanopia', or 'Tritanopia'
%   modType:          Type of cone modulation for hyperspectral cases
%                         'L', 'M', 'S', or 'rand'
%   nSquares:         Number of tiles (used only in hyperspectral generation)
%   constraintWL:     Wavelength that defines projection plane for dichromat simulation
%                         Typically 575 for Deuteranopia
%   Disp:             Structure containing display calibration, cone sensitivities, and image dimensions
%
% Outputs:
%   triLMSCalFormat:  LMS image in calibration format for trichromatic viewer
%   diLMSCalFormat:   LMS image simulated for dichromatic viewer (same format)
%   Disp:             Updated display struct with possibly new dimensions
%
% Description:
%   Depending on the input image, this function performs one of several operations:
%   (1) For 'ishihara', it generates a synthetic Ishihara plate with embedded digit
%   (2) For file inputs (.png, .jpg), it loads and converts the image to LMS
%   (3) For hyperspectral inputs, it calls a renderer to generate simulated LMS scenes
%   In all cases, it produces LMS versions for both trichromat and dichromat viewers.
%
%
% History:
%   04/14/2025  cmd  Initial draft and structure
%
% Examples:
%{
Disp = loadDisplay('ishihara');
testLMS = loadLMSvalues('ishihara','Deuteranopia','M',[],585,1,Disp);

% Check that behavior has not changed since we declared it good.
if (abs(sum(testLMS(:)) - 525.8024)/525.8024 > 1e-4)
    error('No longer get same LMS image returned by loadLMSValues');
end
%}

if strcmp(img,'ishihara')

    % 1 -> gray with missing cone mod
    % 2 -> background random inside with missing cone mod
    % 3 -> LS background, M inside
    % 4 -> like 2 but constrained between .3 and .7 colors so more room for
    %      modulation

    [insideColors, outsideColors] = chooseIshiharaColors(renderType,plateType,Disp);

    % Generate plate now that you have the correct colors
    ishiharaRGB = generateIshiharaPlate('74', insideColors, outsideColors,Disp.m);
    ishiharaRGB = im2double(ishiharaRGB);

    % Plot modified RGB Image 
    figure();imagesc(ishiharaRGB)
    axis square;

    % Put modified image into LMS 
    triRGBCalFormat    = ImageToCalFormat(ishiharaRGB);
    triLMSCalFormat    = Disp.M_rgb2cones * triRGBCalFormat;
    % 
    % M_plane = triLMSCalFormat(2,:);
    % % Assuming M_plane is 1 x N
    % N = length(M_plane);
    % 
    % % Fill in L and S planes with mean of M (or another estimate)
    % L_plane = mean(M_plane) * zeros(1, N);
    % S_plane = mean(M_plane) * zeros(1, N);
    % LMS_full = [L_plane; M_plane; S_plane];
    % M_cones2rgb = inv(Disp.M_rgb2cones);
    % RGB = M_cones2rgb * LMS_full;
    % 
    % rgbImage = reshape(RGB', Disp.m, Disp.n, 3);
    % figure();
    % imshow(rgbImage);



    % Run the modulated image through the linear dichromat simulation
    [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, Disp.grayLMS,  constraintWL, renderType, Disp);

    % Check
    rgb1 = inv(Disp.M_rgb2cones) * diLMSCalFormat;
    image = CalFormatToImage(rgb1,Disp.m,Disp.n);
    figure();imagesc(image); axis square;
    
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)

    % img_rgb = im2double(imread(img));           % Load and convert image to double
    % if Disp.m > 128 || Disp.n > 128
    %     scaleFactor = 0.6;                      % Downsample
    %     img_rgb = imresize(img_rgb, scaleFactor);
    %     [rows, cols, ~] = size(img_rgb);
    %     Disp.m         = cols;
    %     Disp.n         = rows;

    img_rgb = im2double(imread(img));                % Load and convert image to double
    img_rgb = imresize(img_rgb, [128*2, 128*2]);         % Resize to 128x128 pixels

    [rows, cols, ~] = size(img_rgb);
    Disp.m = cols;
    Disp.n = rows;

    % Cal format RGB 
    triRGBCalFormat = ImageToCalFormat(img_rgb); 
    triRGBCalFormat(triRGBCalFormat>1) = 1;
    triRGBCalFormat(triRGBCalFormat<0) = 0;
    % Convert to LMS 
    triLMSCalFormat = Disp.M_rgb2cones * triRGBCalFormat;

    [diLMSCalFormat,M_triToDi]       = DichromSimulateLinear(triLMSCalFormat, Disp.grayLMS,  constraintWL, renderType, Disp);
    
    % check
    rgb1 = inv(Disp.M_rgb2cones) * diLMSCalFormat;
    imageDi = CalFormatToImage(rgb1, Disp.m, Disp.n);
    figure();imagesc(imageDi);

else
    % Get trichromatic (LMS) image
    Disp.m = 32;
    Disp.n = 32;

    % I think t_renderHyperspectralImage is only used in order to create
    % LMS values for the gray image with square isochromatic plates added
    % on. Maybe you can simplify this? Not sure. 
    [triLMSCalFormat,triLMSCalFormat_plate,diLMSCalFormat,diLMSCalFormat_plate] = t_renderHyperspectralImage(img,renderType,constraintWL,nSquares,modType,Disp);
    clear triLMSCalFormat;
    clear diLMSCalFormat;
    triLMSCalFormat = triLMSCalFormat_plate; % do this when you just want to see the isochromatic plate square version (other is just gray)
    diLMSCalFormat  = diLMSCalFormat_plate;
end


end