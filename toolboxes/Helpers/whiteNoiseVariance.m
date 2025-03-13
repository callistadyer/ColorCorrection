
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

% Create a white noise image of size mxn
whiteNoiseImage = rand(Disp.m, Disp.n, 3);
whiteNoiseCalFormat = ImageToCalFormat(whiteNoiseImage);

% 1... MANUAL
% Make hyperspectral img by multiplying primaries * rgb values at each pixel
hyperspectralNoiseCalFormat = Disp.P_monitor * whiteNoiseCalFormat;

% LMS image
triLMSNoiseCalFormat       = Disp.T_cones*hyperspectralNoiseCalFormat;

% Compute total variance 
totalVariance = varianceLMS(varianceType,renderType,[],triLMSNoiseCalFormat);

end
