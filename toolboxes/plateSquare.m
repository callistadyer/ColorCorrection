function image = plateSquare(imgSz, modulation, sideLength)

% Initialize the image to zero
image = zeros(imgSz);

% Define the center and half side length of the square
centerX = imgSz(1) / 2;
centerY = imgSz(2) / 2;
halfSideLength = sideLength / 2;

for x = 1:imgSz(1)
    for y = 1:imgSz(2)
        % Calculate the distance from the center in x and y directions
        distanceX = abs(x - centerX);
        distanceY = abs(y - centerY);
        
        % If the distance is within the noisy radius range, set the pixel
        % to modulation value
        if distanceX < halfSideLength && distanceY < halfSideLength
            if size(modulation,1) > 1 || size(modulation,2) > 1
                image(x, y) = modulation(x,y);
            else
                image(x, y) = modulation;
            end
        end
    end
end

end
