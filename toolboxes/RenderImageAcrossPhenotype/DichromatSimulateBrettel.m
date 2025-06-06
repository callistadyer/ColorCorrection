function [diXYZ] = DichromatSimulateBrettel(xyzImage, cbTypes, rgbImage)
% Simulates color vision for dichromatic viewers using the Brettel, Vienot, Mollon (1997) method.
%
% Syntax:
%   [diXYZ] = DichromatSimulateBrettel(xyzImage, cbTypes,rgbImage)
%
% Inputs:
%   xyzImage: Input XYZ image
%   cbTypes:  Array of integers specifying the types of color blindness to simulate:
%             1 = Protanopia, 2 = Deuteranopia, 3 = Tritanopia.
%             Example: [1, 2] will simulate Protanopia and Deuteranopia.
%   rgbImage: Input RGB image (optional). If empty, defaults to '74.jpg' on the local path.
%
% Outputs:
%   diXYZ:    dichromat XYZ values
%
% Description:
%   This function renders images to simulate their appearance to color-blind viewers,
%   specifically using the Brettel, Vienot, Mollon (1997) algorithm. It handles the
%   transformation from XYZ color space to LMS for both trichromatic and dichromatic
%   vision, and converts back to sRGB for visualization.
%
% History:
%   01/10/2025 - Updated for all cbTypes by cmd
%
% Examples:
%{
% Simulate Protanopia and Deuteranopia:
[srgb_dichromat, lmsDichromat, lmsTrichromat] = Brettel([], [1, 2]);

% Simulate all three types:
[srgb_dichromat, lmsDichromat, lmsTrichromat] = Brettel([], [1, 2, 3]);
%}

%% Input Handling and Setup
% Basically, if xyz values are provided, make sure to use that instead of rgb 
if ~isempty(xyzImage)
    rgbImage = [];
    % disp("Note: XYZ values are being used, NOT rgb values")
end
if nargin < 1 || (isempty(rgbImage) && isempty(xyzImage))
    % Default to '74.jpg' if no input is provided
    rgbImage = imread('macbeth.tif');
    rgbImage = imresize(rgbImage, [128, 128]);
    rgbImage = double(rgbImage)/255;
end

% Turn the image we read into gray (test! gray should be gray in this algorithm)
forceGray = false;
if (forceGray)
    rgbImage = 0.7*ones(size(rgbImage));
end

if nargin < 2 || isempty(cbTypes)
    % Default to all three types if no cbTypes are specified
    cbTypes = [1, 2, 3];
end

% Ensure cbTypes is a row vector
cbTypes = unique(cbTypes(:))';

%%%% Can uncomment this stuff if you wanna do some plotting %%%%
% Create a regular MATLAB figure
% f = figure('Position',[863         899        1289         284]);
% set(f, 'Name', 'Color Blindness Simulations', 'NumberTitle', 'off');
% Display the original RGB image
% subplot(1, length(cbTypes) + 2, 1); % Add the original as the first subplot
% imshow(rgbImage);
% title('Original Image');

% Wavelength range for display 
wls = (400:10:700)';
% Display parameters
d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);

%% Load and Process the Image
% Convert the RGB image into a scene with calibrated display

% Get white point - need to do this manually or matlab does something funky
% that makes our gray image white-ish
rgbWhite(1,1,:) = [1, 1, 1];
sceneWhite = sceneFromFile(rgbWhite, 'rgb', [], d, wls);
whiteXYZ = squeeze(sceneGet(sceneWhite, 'xyz'))';

% Routine for when rgb values are used instead of XYZ (do NOT want to do
% this if we are importing image through Callista's decolorization routine... 
% need to use XYZ values to keep conversions consistent)
if ~ isempty(rgbImage)
    scene = sceneFromFile(rgbImage, 'rgb', [], d, wls);
    xyzImage = sceneGet(scene, 'xyz');
    xyzImage = srgb2xyz(rgbImage);
    whiteXYZ = srgb2xyz(rgbWhite);

    % Convert the image from XYZ to LMS
    lmsTrichromat = xyz2lms(xyzImage, [], whiteXYZ);
    xyzTrichromat = imageLinearTransform(lmsTrichromat, colorTransformMatrix('lms2xyz'));
    srgbTrichromat = xyz2srgb(xyzTrichromat);
    % subplot(1, length(cbTypes) + 2, 2); % Add each cbType simulation
    % imshow(srgbTrichromat);
end



%% Simulate Color Blindness for Specified Types
% srgb_dichromat = cell(length(cbTypes), 1);
lmsDichromat = cell(length(cbTypes), 1);
    
% Simulate for each specified type of color blindness
cbTypeNames = {'Protanopia', 'Deuteranopia', 'Tritanopia'};
for idx = 1:length(cbTypes)
    cbType = cbTypes(idx);
    
    % Convert to LMS for the current dichromatic viewer
    lmsDichromat{idx} = xyz2lms(xyzImage, cbType, whiteXYZ);
    
    % Convert the dichromatic LMS values back to XYZ
    diXYZ = imageLinearTransform(lmsDichromat{idx}, colorTransformMatrix('lms2xyz'));
    
    % % Can uncomment this if you are not using decolorization routine (aka using rgb vals instead of xyz)  
    % Convert the color-blind XYZ values to sRGB
    % srgb_dichromat{idx} = xyz2srgb(cbXYZ);
    % subplot(1, length(cbTypes) + 2, idx + 2); % Add each cbType simulation
    % imshow(srgb_dichromat{idx});
    % title(['Simulated ' cbTypeNames{cbType}]);
end

end
