function modulation = getModulation(rgbImageCalFormat,renderType,modulationDirection,m,n)

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
    modulation(i) = MaximizeGamutContrast(modulationDirection,rgbImageCalFormat(:,i)); % bg is in rgb cal format
end

% Use these two lines to use the same modulation for all pixels... will need to take the minimum mod of all 
nonzeroMod = modulation(modulation>0);
modulation = min(nonzeroMod(:));

% If there is more than 1 value for the modulation (ie. different for each pixel) 
if size(modulation,1) > 1 || size(modulation,2) > 1
    modulation = CalFormatToImage(modulation,m,n);
end



end