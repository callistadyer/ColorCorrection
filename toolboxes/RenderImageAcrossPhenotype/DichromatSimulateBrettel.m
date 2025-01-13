function [srgb_dichromat, lmsDichromat, lmsTrichromat] = DichromatSimulateBrettel(rgbImage, cbTypes)
% Simulates color vision for dichromatic viewers using the Brettel, Vienot, Mollon (1997) method.
%
% Syntax:
%   [srgb_dichromat, lmsDichromat, lmsTrichromat] = DichromatSimulateBrettel(rgbImage, cbTypes)
%
% Inputs:
%   rgbImage: Input RGB image (optional). If empty, defaults to '74.jpg' on the local path.
%   cbTypes:  Array of integers specifying the types of color blindness to simulate:
%             1 = Protanopia, 2 = Deuteranopia, 3 = Tritanopia.
%             Example: [1, 2] will simulate Protanopia and Deuteranopia.
%
% Outputs:
%   srgb_dichromat: Simulated sRGB images for specified dichromatic viewers (cell array).
%   lmsDichromat:   LMS values for the specified dichromatic viewers (cell array).
%   lmsTrichromat:  LMS values for trichromatic viewers (normal vision).
%
% Description:
%   This function renders images to simulate their appearance to color-blind viewers,
%   specifically using the Brettel, Vienot, Mollon (1997) algorithm. It handles the
%   transformation from XYZ color space to LMS for both trichromatic and dichromatic
%   vision, and converts back to sRGB for visualization.
%
% History:
%   01/10/2025 - Updated to dynamically use subplots for cbTypes by cmd.
%
% Examples:
%{
% Simulate Protanopia and Deuteranopia:
[srgb_dichromat, lmsDichromat, lmsTrichromat] = Brettel([], [1, 2]);

% Simulate all three types:
[srgb_dichromat, lmsDichromat, lmsTrichromat] = Brettel([], [1, 2, 3]);
%}

%% Input Handling and Setup
if nargin < 1 || isempty(rgbImage)
    % Default to '74.jpg' if no input is provided
    rgbImage = imread('74.jpg');
    rgbImage = imresize(rgbImage, [128, 128]);

end

if nargin < 2 || isempty(cbTypes)
    % Default to all three types if no cbTypes are specified
    cbTypes = [1, 2, 3];
end

% Ensure cbTypes is a row vector
cbTypes = unique(cbTypes(:))';

% Create a regular MATLAB figure

f = figure('Position',[863         899        1289         284]);
set(f, 'Name', 'Color Blindness Simulations', 'NumberTitle', 'off');

% Display the original RGB image
subplot(1, length(cbTypes) + 1, 1); % Add the original as the first subplot
imshow(rgbImage);
title('Original Image');

% Wavelength range for display 
wls = (400:10:700)';

d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);

%% Load and Process the Image
% Convert the RGB image into a scene with calibrated display
scene = sceneFromFile(rgbImage, 'rgb', [], d, wls);

% Extract XYZ values from the scene
imgXYZ = sceneGet(scene, 'xyz');
whiteXYZ = sceneGet(scene, 'illuminant xyz');

% Convert the image from XYZ to LMS 
lmsTrichromat = xyz2lms(imgXYZ, [], whiteXYZ);

%% Simulate Color Blindness for Specified Types
srgb_dichromat = cell(length(cbTypes), 1);
lmsDichromat = cell(length(cbTypes), 1);

% Simulate for each specified type of color blindness
cbTypeNames = {'Protanopia', 'Deuteranopia', 'Tritanopia'};
for idx = 1:length(cbTypes)
    cbType = cbTypes(idx);
    
    % Convert to LMS for the current dichromatic viewer
    lmsDichromat{idx} = xyz2lms(imgXYZ, cbType, whiteXYZ);
    
    % Convert the dichromatic LMS values back to XYZ
    cbXYZ = imageLinearTransform(lmsDichromat{idx}, colorTransformMatrix('lms2xyz'));
    
    % Convert the color-blind XYZ values to sRGB
    srgb_dichromat{idx} = xyz2srgb(cbXYZ);
    
    subplot(1, length(cbTypes) + 1, idx + 1); % Add each cbType simulation
    imshow(srgb_dichromat{idx});
    title(['Simulated ' cbTypeNames{cbType}]);
end

end
