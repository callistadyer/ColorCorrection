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
%   setType:  Numeric identifier specifying the stimulus variation, e.g.:
%               - Number of squares for 'gray'
%               - Plate type for 'ishihara'
%   m:        Image dimension 1
%   n:        Image dimension 2
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

% Initialize
imgParams = struct();
describe = struct();

% Default size
if isempty(m), m = 128; end
if isempty(n), n = 128; end

% Set image description depending on image type
if strcmpi(img, 'gray')
    if isempty(setType), setType = 1; end
    if ~ismember(setType, 1:10)
        error('ERROR: undefined setType for gray. Choose setType between 1–10.');
    end
    describe.imgType   = 'gray';
    describe.nSquares  = setType;

elseif strcmpi(img, 'ishihara')
    if isempty(setType), setType = 1; end
    if ~ismember(setType, 1:4)
        error('ERROR: undefined setType for ishihara. Choose setType between 1–4.');
    end
    describe.imgType   = 'ishihara';
    describe.plateType = setType;

elseif endsWith(img, {'.png', '.jpg'}, 'IgnoreCase', true)
    if isempty(setType), setType = 1; end
    describe.imgType   = 'natural image (.png or .jpg)';
    describe.filename  = img;

else
    error('ERROR: Unrecognized img type "%s".', img);
end

imgParams.m = m;
imgParams.n = n;
imgParams.describe = describe;

end
