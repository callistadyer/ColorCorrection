
function totalVariance = whiteNoiseVariance(varianceType,renderType,Disp)
% Calculates total variance to normalize variance of given image size 
%
% Syntax:
%   totalVariance = whiteNoiseVariance(Disp)
%
% Inputs:
%   varianceType:  Type of variance computation
%                  "newConeVar"
%   Disp:          Display parameters
% Outputs:
%   totalVariance: Variance normalizer term 
%
% Optional key/value pairs:
%   None
%

% Random seed to calculate total variance 
rng(2);
% s = rng
% Create a white noise image of size mxn
whiteNoiseImage = rand(Disp.m, Disp.n, 3);
whiteNoiseCalFormat = ImageToCalFormat(whiteNoiseImage);

% 1... MANUAL
% Make hyperspectral img by multiplying primaries * rgb values at each pixel
hyperspectralNoiseCalFormat = Disp.P_monitor * whiteNoiseCalFormat;

% LMS image (1)
triLMSNoiseCalFormat       = Disp.T_cones*hyperspectralNoiseCalFormat;

% Compute total variance 
totalVariance = varianceLMS(varianceType,renderType,[],triLMSNoiseCalFormat);

end

% ALTERNATE IMAGE CREATION
% 2... CREATE SCENE VERSION
% Image format
%hyperspectralNoiseImage     = CalFormatToImage(hyperspectralNoiseCalFormat,Disp.m,Disp.n);
%
% NOTE: whiteNoiseImage needs to be ungammacorrected before sceneToFile will do
% what you want. For example, whiteNoiseImage.^2.
% scene = sceneFromFile(whiteNoiseImage, 'rgb', [], Disp.d, Disp.wls);
% % Get XYZ values (needed for dichromat simulation algorithm from Brettel)
% imgXYZ = sceneGet(scene,'xyz');
% % Hyperspectral image (needed for calculating cone responses)
% hyperspectralImage2     = double(sceneGet(scene,'energy'));
% hyperspectralCalFormat2 = ImageToCalFormat(hyperspectralImage2);
% 
% % LMS image (2)
% triLMSNoiseCalFormat2   = Disp.T_cones*hyperspectralCalFormat2;
