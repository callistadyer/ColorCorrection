function imgRGB = generateIshiharaPlate(textStr, insideColors, outsideColors, outputSize)
% Renders an Ishihara-style dot image with a number embedded
%
% Syntax:
%   imgRGB = generateIshiharaPlate(textStr, insideColors, outsideColors, outputSize)
%
% Inputs:
%   textStr       - string. The numeric string to embed in the image (e.g., '74')
%   insideColors  - 3x3 matrix (3 possible colors). RGB colors for dots inside the number
%   outsideColors - 3x3 matrix (3 possible colors). RGB colors for dots outside the number
%   outputSize    - Final image resolution (e.g., 128)
%
% Output:
%   imgRGB        - outputSize x outputSize x 3 RGB IMAGE. Final rendered Ishihara plate
%
% Example:
%{ 
insideColors = [
    0.1 0.8 0.5;
    0.2 0.6 0.7;
    0.3 0.9 0.4
];
outsideColors = [
    1.0 0.6 0.3;
    0.9 0.7 0.2;
    1.0 0.5 0.5
];
img = generateIshiharaPlate('74', insideColors, outsideColors, 128);
figure();
imagesc(img)
axis equal;
title(['Isochromatic plate', textStr]);
%}

% Random seed for reproducibility
rng(0);

% Size of working image (not end result image)
imgSize = 512;

% Radius of the circular plate
plateRadius = imgSize / 2 - 5;

% Center coordinate of the plate
plateCenter = [imgSize/2, imgSize/2];

%%%%%%%%%%%%%%%%%%%%%%% Plate parameters %%%%%%%%%%%%%%%%%%%%%%%
numDots = [80, 100, 550, 800, 400]; % Number of dots for each size
radii   = [12, 10,  6,   5,   3];   % Radii for each dot 
buffer = 1.2;                       % Some space between dots
fontSize = 310;                     % Font size of embedded number

% Plate parameters
% numDots = [10, 10, 150, 250, 450]*.8;
% radii   = [12, 10,  6,   5,   3]*.2;
% buffer = 1.2*.2;
% fontSize = 310*.3; % Font size of embedded number


% This part turns the grayscale image into a binary mask (numberMask)
% where the number appears as ones (white) and the background 
% as zeros (black). This mask is later used to determine whether 
% each dot in the plate lies inside or outside the number
% This helps us do the coloring correctly 
fig = figure('Visible','off');
axes('Position',[0 0 1 1]);
% Write the number 
% text(0.5, 0.5, textStr, 'FontSize', fontSize, 'FontWeight','bold', ...
%     'HorizontalAlignment','center', 'VerticalAlignment','middle');
text(0.5, 0.5, textStr, 'FontSize', fontSize, 'FontWeight','bold', ...
    'FontName','Arial Rounded MT Bold', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
xlim([0 1]); ylim([0 1]);
axis off
frame = getframe(gca);
img = rgb2gray(frame.cdata);
numberMask = imbinarize(imresize(img, [imgSize imgSize]));
close(fig);

% % Uncomment to look at mask
% figure();
% imshow(numberMask);
% title(['Binary Mask for "', textStr, '"']);

f = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, imgSize, imgSize]);

% create axes within the figure that match the image dimensions
ax = axes('Units', 'pixels', 'Position', [0, 0, imgSize, imgSize]);

hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');

% set the axes limits to match the canvas size
xlim(ax, [0 imgSize]);
ylim(ax, [0 imgSize]);

% initialize an empty array to store placed dot positions and radii
positions = [];

% Check if dot can be placed
function isValid = isValidPlacement(x_center, y_center, radius, positions)
    % Compute distance from the center of the plate
    distFromCenter = sqrt((x_center - plateCenter(1))^2 + (y_center - plateCenter(2))^2);
    % Inside the cirlce plate boundary?
    if distFromCenter > (plateRadius - radius)
        isValid = false;
        return;
    end
    
    % Overlapping? 
    for i = 1:size(positions, 1)
        % Get the distance between the new and all of the old circle centers 
        dx = x_center - positions(i,1);
        dy = y_center - positions(i,2);
        dist = sqrt(dx^2 + dy^2);
        % If the distance is too close, don't place the new circle
        if dist < (radius + positions(i,3) + buffer)
            isValid = false;
            return;
        end
    end
        isValid = true;
end

% Function that places circles
function positions = placeCircles(positions, numCircles, radius)
    count = 0; attempts = 0; attemptsLimit = 10050;
    
    % Place dots so long as you have more dots to place... or if you've
    % tried for a long time, then give up
    while count < numCircles && attempts < attemptsLimit
        % pick a random (x, y) location
        x = rand * imgSize;
        y = rand * imgSize;
        xi = round(x); yi = round(y);
        
        % reject if coordinates are out of bounds (between 1 and image size)
        if xi < 1 || xi > imgSize || yi < 1 || yi > imgSize
            attempts = attempts + 1;
            continue;
        end
        
        % check if placement is valid
        if isValidPlacement(x, y, radius, positions)
            % Determine if the dot is inside the number shape
            % Plotting coordinates start from the bottom-left BUT 
            % image matrices like numberMask start from the top-left, 
            % so the y-coordinate must be flipped when indexing into the mask. 
            % So we compute imgSize - yi + 1; for example, if yi = 50 in a 
            % 512×512 image, the corresponding row in the image is 512 - 50 + 1 = 463
            % "Give me the pixel at row y (height), column x (width)"

            % numberMask has 0s in the number and 1s elsewhere... so you
            % have to use the not ~ when checking if it is inside (~0 = ~FALSE = TRUE)
            inside = ~numberMask(imgSize - yi + 1, xi);

            % Color based on inside or outside
            if inside == 1
                color = insideColors(randi(size(insideColors,1)), :);
            else
                color = outsideColors(randi(size(outsideColors,1)), :);
            end
            
            % Draw this dot
            theta = linspace(0, 2*pi, 30);
            xCirc = x + radius * cos(theta); % x coordinates
            yCirc = y + radius * sin(theta); % y coordinates
            fill(xCirc, yCirc, color, 'EdgeColor', color, 'LineWidth', 1);
            
            % Record position and radius
            positions = [positions; x, y, radius];
            
            % Successful dot counter
            count = count + 1;
        end
        attempts = attempts + 1;
    end
end

% Loop over each dot size and place the specified number of dots
for i = 1:length(radii)
    positions = placeCircles(positions, numDots(i), radii(i));
end

frame = getframe(ax);
% cdata has the rgb information
% resize the image so it's not sooo slow in the correction algorithm
imgRGB = imresize(frame.cdata, [outputSize outputSize]);
close(f);

bPLOT = 0;
if bPLOT == 1
    figure();
    imagesc(imgRGB);
    axis square;
    title(['Isochromatic plate', textStr]);
end

end
