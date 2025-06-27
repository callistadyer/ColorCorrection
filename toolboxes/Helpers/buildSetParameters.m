function [setParams] = buildSetParameters(img,setType)

setParams = struct();

if strcmp(img, 'gray')
    switch setType
        case 1
            setParams.nSquares = 1;
        case 2
            setParams.nSquares = 2;
        otherwise
            error('ERROR: undefined setType. Make sure your setType is compatible with the img type')
    end

elseif strcmp(img, 'ishihara')
    %   plateType:    - Double. Only relevant for ishihara plates. 
%                        1 -> gray with missing cone mod
%                        2 -> background random inside with missing cone mod
%                        3 -> LS background, M inside
%                        4 -> like 2 but constrained between .3 and .7 colors so more room for modulation

    switch setType
        case 1
            setParams.plateType = 1;
        case 2
            setParams.plateType = 2;    
        case 3
            setParams.plateType = 3; 
        case 4
            setParams.plateType = 4;  
        otherwise
            error('ERROR: undefined setType. Make sure your setType is compatible with the img type')
    end
elseif endsWith(img, '.png', 'IgnoreCase', true) || endsWith(img, '.jpg', 'IgnoreCase', true)

end
