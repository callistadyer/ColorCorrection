function imgParams = buildSetParameters(img,setType,m,n)

% buildSetParameters  Define stimulus-specific parameters, including image size, based on image type and setType.
%
% Syntax:
%   setParams = buildSetParameters(img, setType, m, n)
%
% Inputs:
%   img:      String specifying image type or filename. Can be:
%               - 'gray'      : gray squares stimulus
%               - 'ishihara'  : Ishihara plates
%               - filename    : external image file (.png, .jpg)
%
%   setType:  Numeric identifier specifying the stimulus variation, e.g.:
%               - Number of squares for 'gray'
%               - Plate type for 'ishihara'
%
%   m:        Numeric. Desired image height in pixels. If empty or not provided, defaults to 128.
%
%   n:        Numeric. Desired image width in pixels. If empty or not provided, defaults to 128.
%
% Outputs:
%   setParams: Struct with stimulus-specific parameters:
%                - nSquares  : number of squares (for 'gray')
%                - plateType : Ishihara plate type (for 'ishihara')
%                - m, n      : image size (height & width in pixels)
%
% Examples:
%{
    % Gray squares with 2 squares, default size 128×128
    setParams = buildSetParameters('gray', 2, [], []);
    
    % Ishihara plate 3, custom size 256×256
    setParams = buildSetParameters('ishihara', 3, 256, 256);
    
    % External image 'flower.png', default size 128×128
    setParams = buildSetParameters('flower.png', [], [], []);
%}

% Initialize empty struct
imgParams = struct();

if isempty(m)
    m = 128;
end

if isempty(n)
    n = 128;
end

if strcmpi(img, 'gray')
    % For gray stimuli
    if isempty(setType), setType = 1; end
    switch setType
        case {1, 2, 3, 4, 5}
            imgParams.nSquares = setType;
        otherwise
            error('ERROR: undefined setType for gray. Choose setType between 1–5.');
    end

elseif strcmpi(img, 'ishihara')
    % For Ishihara plates
    if isempty(setType), setType = 1; end
    switch setType
        case {1, 2, 3, 4}
            imgParams.plateType = setType;
        otherwise
            error('ERROR: undefined setType for ishihara. Choose setType between 1–4.');
    end


elseif endsWith(img, {'.png', '.jpg'}, 'IgnoreCase', true)
    % For external images
    if isempty(setType), setType = 1; end
    switch setType
        case 1
            imgParams.plateType = setType;
        case 2
            imgParams.plateType = setType;
        otherwise
            error('ERROR: undefined setType for ishihara. Choose setType between 1–4.');
    end

else
    error('ERROR: Unrecognized img type "%s".', img);
end

% Size of image
imgParams.m = m;
imgParams.n = n;

imgParams.img = img;
imgParams.setType = setType;


end
