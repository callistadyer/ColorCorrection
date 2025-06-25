function [setParams] = buildSetParameters(img,setType)

setParams = struct();

if strcmp(img, 'gray')
    switch setType
        case 'gray1'
            setParams.nSquares = 1;
        case 'gray2'
            setParams.nSquares = 2;
        otherwise
            error('ERROR: undefined setType. Make sure your setType is compatible with the img type')
    end

elseif strcmp(img, 'ishihara')
    switch setType
        case 'ishihara1'
            setParams.plateType = 1;
        case 'ishihara2'
            setParams.plateType = 2;    
        otherwise
            error('ERROR: undefined setType. Make sure your setType is compatible with the img type')
    end
end
