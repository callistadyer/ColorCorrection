function [triLMScalFormat,triLMSCalFormat_plate] = generateGrayImage(nSquares,modType,Disp,imgParams)    
% Demonstrate how to read, add cone directed info, and render for tri- and dichromat
%
% Syntax:
%   [triLMScalFormat,triLMSCalFormat_plate] = generateGrayImage(nSquares,modType,Disp,imgParams)  
%
% Description:
%
% Inputs:
%   nSquares:     - Double. Number of squares in the isochromatic plate.
%   modType:      - String. Type of isochromatic plate modulation. If you want to
%                   modulate in the direction of a cone, then write which 
%                   cone (e.g., "M"). "rand" just gives random colors in
%                   the squares.
%                       'rand'
%                       'M'
%                       'L'
%                       'S'
%  Disp:          - Struct. Display parameters. If you dont have it, leave
%                   it empty [] and the code will handle it. 
%
% Outputs:
%   None
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code
%
% History
%   07/30/2024  dhb, cmd  Initial go.
%   08/27/2024  dhb, cmd  It's working, cleaning up.
%   06/16/2025  cmd, cleaning up after lots of project changes
%
% Examples:
%{
Disp = loadDisplay();
imgParams = buildSetParameters('gray',10,64,64);
[triLMScalFormat,triLMSCalFormat_plate] = generateGrayImage(10,'M',Disp,imgParams);
triRGBCalFormat_plate = LMS2RGBCalFormat(triLMSCalFormat_plate,Disp,imgParams);
GrayImage = CalFormatToImage(triRGBCalFormat_plate,64,64);
figure(); imagesc(GrayImage); axis square;
%}
% Create gray hyperspectral image
[whitergbLinCalFormat] = ImageToCalFormat(ones(imgParams.m,imgParams.n));

% Gray 0.5 rgb at each pixel in image
grayrgbLinCalFormat       = (0.5.*(repmat(whitergbLinCalFormat,3,1)));

% Make hyperspectral img by multiplying primaries * rgb values at each pixel
% This is a weighted sum of primaries
% hyperspecGrayCalFormat = Disp.P_monitor * grayImgCalFormat;

% Image format
% grayrgbImage    = CalFormatToImage(grayrgbLinCalFormat,imgParams.m,imgParams.n);

% Render lms image that a trichromat sees
triLMScalFormat = rgbLin2LMSCalFormat(grayrgbLinCalFormat,Disp);
triLMSImage = CalFormatToImage(triLMScalFormat,imgParams.m,imgParams.n);
trirgbLinCalFormat = Disp.M_cones2rgb * triLMScalFormat;

% Get modulation for isochromatic plate modulation.
% This function is taking triRGBCalFormat and using MaximizeGamutContrast to determine how much we can move
% in a given cone direction (if specified, the direction of the missing cone) without going out of gamut
for i = 1:nSquares
    [lmsModulationImgFormat(:,:,:,i)] = getDichromatConfusionModulation(trirgbLinCalFormat,modType,Disp,imgParams);
end

% Create isochromatic plates
% This is taking in the original hyperspectralImage and then adding the
% modulation in squares to the it
[~, triLMSCalFormat_plate] = isochromaticPlates(triLMSImage,lmsModulationImgFormat,Disp,imgParams, ...
    'verbose',true);


end