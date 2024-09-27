function image = plateO(imgSz,modulation,outerRadius)


image = zeros(imgSz);

% Define the center and base radius of the circle
centerX = imgSz(1) / 2;
centerY = imgSz(2) / 2;
noiseLevel = 0; % Noise level (try 30 for a fuzzy looking version)

% Loop through each pixel in the image
for x = 1:imgSz(1)
    for y = 1:imgSz(2)
        % Calculate distance from the center
        distance = sqrt((x - centerX)^2 + (y - centerY)^2);
        
        % Add noise to the outer and inner radius
        noisyOuterRadius = outerRadius + noiseLevel * (rand - 0.5); % Outer radius with noise
        noisyInnerRadius = outerRadius - 30 + noiseLevel * (rand - 0.5); % Inner radius with noise
        
        % If the distance is within the noisy radius range, set the pixel
        % to modulation value
        if distance < noisyOuterRadius && distance > noisyInnerRadius
            image(x, y) = modulation;
        end
    end
end



end
