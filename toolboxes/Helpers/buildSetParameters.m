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
