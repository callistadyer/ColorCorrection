function scaleFactor = getModulation(rgbImageCalFormat,renderType,takeMin,m,n)

% function that calculates the modulation in the L M or S cone based on the
% input image and the gamut limitations

% Which cone do you want to modulate
switch (renderType)
    case 'Deuteranopia' % m cone deficiency
        modulationDirection = [0 1 0]';
    case 'Protanopia'   % l cone deficiency
        modulationDirection = [1 0 0]';
    case 'Tritanopia'   % s cone deficiency
        modulationDirection = [0 0 1]';
end

% Background is the value of each pixel: this determines cone contrast separately for each pixel
for i = 1:size(rgbImageCalFormat,2)
    scaleFactor(i) = MaximizeGamutContrast(modulationDirection,rgbImageCalFormat(:,i)); % bg is in rgb cal format
end

if takeMin == 1
% Use these two lines to use the same modulation for all pixels... will need to take the minimum mod of all 
nonzeroMod = scaleFactor(scaleFactor>0);  % modulations above 0 
scaleFactor = min(nonzeroMod(:));        % minimum modulation of all pixels
end

% If there is more than 1 value for the modulation (ie. different for each pixel) 
if size(scaleFactor,1) > 1 || size(scaleFactor,2) > 1
    scaleFactor = CalFormatToImage(scaleFactor,m,n);
end

end