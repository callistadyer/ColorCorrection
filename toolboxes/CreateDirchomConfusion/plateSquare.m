function image = plateSquare(imgSz, modulation, nSquares)

    image = zeros(imgSz);
    
    % CENTER OF CIRCLE
    centerX = imgSz(1) / 2;
    centerY = imgSz(2) / 2;
    
    % RADIUS OF CIRCLE
    radius = min(imgSz(1:2)) * 0.4;
    
    % STEPS TO TAKE AROUND THE CIRCLE, DEPENDS ON HOW MANY SQUARES
    angleStep = 2 * pi / nSquares;
    
    % LOOP FOR EACH SQUARE
    for k = 1:nSquares
        % ANGLE FOR SQUARE K (FIRST ONE IS AT TOP SO ANGLE = 0)
        angle = (k - 1) * angleStep;
        
        % POSITION OF SQUARE (CENTER) IN X AND Y 
        squareCenterX = round(centerX + radius * cos(angle));
        squareCenterY = round(centerY - radius * sin(angle)); 
        
        % HOW BIG SHOULD THE SQUARE BE?
        halfSideLength = round(imgSz(1) * 0.05); 
        
        % SQUARE X AND Y VALUES
        xMin = round(squareCenterX - halfSideLength);
        xMax = round(squareCenterX + halfSideLength);
        yMin = round(squareCenterY - halfSideLength);
        yMax = round(squareCenterY + halfSideLength);
        
        % INSERT MODULATION AT DESIGNATED SQUARE LOCATIONS
        for z = 1:imgSz(3)
            for x = xMin:xMax
                for y = yMin:yMax
                    image(x, y, z) = modulation(x,y,z,k);
                end
            end
        end
    end
end




% function image = plateSquare(imgSz, modulations, nSquares)
% 
% % function image = plateSquare(imgSz, modulation, sideLength)
% 
% % Initialize the image to zero
% image = zeros(imgSz);
% 
% % Define the center and half side length of the square
% centerX = imgSz(1) / 2;
% centerY = imgSz(2) / 2;
% halfSideLength = sideLength / 2;
% 
% for z = 1:imgSz(3)
%     for x = 1:imgSz(1)
%         for y = 1:imgSz(2)
%             % Calculate the distance from the center in x and y directions
%             distanceX = abs(x - centerX);
%             distanceY = abs(y - centerY);
% 
%             % If the distance is within the noisy radius range, set the pixel
%             % to modulation value
%             if distanceX < halfSideLength && distanceY < halfSideLength
%                 if size(modulation,1) > 1 || size(modulation,2) > 1
%                     image(x, y, z) = modulation(x,y,z);
%                 else
%                     image(x, y, z) = modulation;
%                 end
%             end
%         end
%     end
% end
% 
% end
