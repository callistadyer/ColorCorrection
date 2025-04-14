function [triLMScalFormat,triLMSCalFormat_plate,diLMScalFormat,diLMScalFormat_plate,Disp,modDirection] = t_renderHyperspectralImage(image,renderType,bPLOTscatter,nSquares,modType)
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
%   renderType    - String. Type of dichromat.  Options are:
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   bPLOTscatter  - Boolean. Plot scatter plots of lms values (1 or 0)
%
% Outputs:
%   None
%
% Optional key/value pairs:
%   None
%
% Examples are included within the code

% History
%   07/30/2024  dhb, cmd  Initial go.
%   08/27/2024  dhb, cmd  It's working, cleaning up.
% Examples:
%{
t_renderHyperspectralImage([],'Deuteranopia'  ,0,0,10,'rand')
t_renderHyperspectralImage('scene1.mat','Protanopia'  ,0,0,10)    
[lmsImageCalFormat,lmsModuledCalFormat] = t_renderHyperspectralImage('scene2.mat','Deuteranopia',0,1,0,10,'rand');    
t_renderHyperspectralImage('scene3.mat','Tritanopia'  ,0,0,10,'rand')    
t_renderHyperspectralImage('scene4.mat','Deuteranopia',0,0,10,'rand')  
t_renderHyperspectralImage('scene5.mat','Deuteranopia',0,0,10,'rand')    
t_renderHyperspectralImage('gray','Deuteranopia',0,0,10,'rand')           
%}

%% Load hyperspectral image data
[hyperspectralImage Disp] = loadImage(image);

% Get cone responses for every pixel of the hyperspectral image
[hyperspectralImageCalFormat,m,n] = ImageToCalFormat(hyperspectralImage);

% Render lms image that a trichromat sees.  Convert from cal format to image as well.
triLMScalFormat = Disp.T_cones*hyperspectralImageCalFormat;

M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);
triRGBCalFormat = M_cones2rgb * triLMScalFormat;

% [triRGBCalFormat,triRGBLinCalFormat,scaleFactor] = LMS2RGBCalFormat(triLMScalFormat,Disp);
triRGBImgFormat                                  = CalFormatToImage(triRGBCalFormat,m,n);


% Note: if you change modType to be "Deuteranopia" etc, then the following
% code will be overwritten inside of getDichromatConfusionModulation.m ...
% the vectors will be changed to the missing cone direction. otherwise,
% random color directions will be used

vectors = rand(nSquares, 3); % Uniform random values between 0 and 1
magnitudes = sqrt(sum(vectors.^2, 2)); % Compute the magnitudes
modDirection = vectors ./ magnitudes;  % Normalize to unit vectors

% Get modulation for isochromatic plate modulation.
% This function is taking rgbLinImageCalFormat2 and using MaximizeGamutContrast to determine how much we can move
% in a given cone direction (specifically, the direction of the missing cone) without going out of gamut
for i = 1:nSquares
    [lmsModulationImgFormat(:,:,:,i)] = getDichromatConfusionModulation(triRGBCalFormat,modDirection(i,:)',modType,Disp,0);
end

% Create isochromatic plates
[triRGBCalFormat_plate, triLMSCalFormat_plate] = isochromaticPlates(image,renderType,lmsModulationImgFormat,Disp,nSquares, ...
    'verbose',true);
triRGBImgFormat_plate = CalFormatToImage(triRGBCalFormat_plate,Disp.m,Disp.n);

% Params
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);
grayRGB = [0.5 0.5 0.5]';
grayLMS = M_rgb2cones*grayRGB;
constraintWl = 585;

% Dichromat manipulation (push trichromat LMS image into this function to get out LMS of dichromat)
[diLMScalFormat,M_triToDi] = DichromSimulateLinear(triLMScalFormat, grayLMS,  constraintWl, renderType, Disp);
[diLMScalFormat_plate,M_triToDi] = DichromSimulateLinear(triLMSCalFormat_plate, grayLMS,  constraintWl, renderType, Disp);

end