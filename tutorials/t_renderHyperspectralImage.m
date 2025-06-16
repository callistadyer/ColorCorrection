function [triLMScalFormat,triLMSCalFormat_plate,diLMScalFormat,diLMScalFormat_plate] = t_renderHyperspectralImage(image,renderType,constraintWl,nSquares,modType,Disp)    
% Demonstrate how to read, add cone directed info, and render for tri- and dichromat
%
% Syntax:
%   [triLMScalFormat,triLMSCalFormat_plate,diLMScalFormat,diLMScalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(image,renderType,bPLOTscatter,nSquares,modType)
%
% Description:
%
% Inputs:
%   image:        - String. Name of image to be rendered. If passed as the empty matrix, you get a
%                   hyperspectral image of some stuffed animals. Some other options are
%                       'sceneN.mat' - N is 1 to 5. One of our hypespectral scenes.
%                       'gray'       - Gray spatially uniform field.                       
%   renderType:   - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   constraintWl: - Double. Constraint wavelength for dichromat simulation.
%                   Default this to 585.
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
[triLMScalFormat,triLMSCalFormat_plate,diLMScalFormat,diLMScalFormat_plate] = t_renderHyperspectralImage('gray','Deuteranopia',585,10,'rand',[]); 
% Visualize:
Disp = loadDisplay('gray');
RBG = LMS2RGBCalFormat(triLMSCalFormat_plate,Disp);
figure(); imagesc(CalFormatToImage(RBG,Disp.m,Disp.n)); axis square;       
%}

if isempty(Disp)
    % Load display 
Disp = loadDisplay(image);
end

% Load hyperspectral image data
% No longer using loadImage. Basically only use it for getting gray image.
% [hyperspectralImage] = loadImage(image);

% Create gray hyperspectral image
[grayImgCalFormat,m,n] = ImageToCalFormat(ones(Disp.m,Disp.n));
% Gray 0.5 rgb at each pixel in image
grayImgCalFormat       = (0.5.*(repmat(grayImgCalFormat,3,1)));

d = displayCreate('LCD-Apple');
P_monitor = SplineSrf(displayGet(d,'wave'),displayGet(d,'spd'),Disp.wls);
% Make hyperspectral img by multiplying primaries * rgb values at each pixel
% This is a weighted sum of primaries
hyperspecGrayCalFormat = P_monitor * grayImgCalFormat;

% Image format
hyperspectralImage     = CalFormatToImage(hyperspecGrayCalFormat,m,n);

% Get cone responses for every pixel of the hyperspectral image
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Render lms image that a trichromat sees
triLMScalFormat = Disp.T_cones*hyperspectralImageCalFormat;
triRGBCalFormat = Disp.M_cones2rgb * triLMScalFormat;

% Get modulation for isochromatic plate modulation.
% This function is taking triRGBCalFormat and using MaximizeGamutContrast to determine how much we can move
% in a given cone direction (if specified, the direction of the missing cone) without going out of gamut
for i = 1:nSquares
    [lmsModulationImgFormat(:,:,:,i)] = getDichromatConfusionModulation(triRGBCalFormat,modType,Disp,0);
end

% Create isochromatic plates
% This is taking in the original hyperspectralImage and then adding the
% modulation in squares to the it
[triRGBCalFormat_plate, triLMSCalFormat_plate] = isochromaticPlates(image,renderType,hyperspectralImage,lmsModulationImgFormat,Disp,nSquares, ...
    'verbose',true);
triRGBImgFormat_plate = CalFormatToImage(triRGBCalFormat_plate,Disp.m,Disp.n);

% Dichromat manipulation (push trichromat LMS image into this function to get out LMS of dichromat)
[diLMScalFormat,M_triToDi]       = DichromSimulateLinear(triLMScalFormat, Disp.grayLMS,  constraintWl, renderType, Disp);
[diLMScalFormat_plate,M_triToDi] = DichromSimulateLinear(triLMSCalFormat_plate, Disp.grayLMS,  constraintWl, renderType, Disp);

end