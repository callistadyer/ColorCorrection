
function totalVariance = whiteNoiseVariance(varianceType,renderType,T,Disp)
% Calculates total variance to normalize variance of given image size 
%
% Syntax:
%   totalVariance = whiteNoiseVariance(varianceType,renderType,Disp)
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
triLMSNoiseCalFormat         = Disp.T_cones*hyperspectralNoiseCalFormat;

% Corrected dichromat image
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
M_cones2rgb = inv(M_rgb2cones);
grayRGB = [0.5 0.5 0.5]';
grayLMS = M_rgb2cones*grayRGB;


triLMSNoiseCalFormatContrast = (triLMSNoiseCalFormat - grayLMS)./grayLMS;
triLMSNoiseCalFormatContrastnew = T * triLMSNoiseCalFormatContrast;
% Compute total variance

switch (varianceType)
    case 'LMdifferenceContrast'
        LMSold = triLMSNoiseCalFormatContrast;
        LMSnew = triLMSNoiseCalFormatContrastnew;
    case 'delta'
        LMSold = triLMSNoiseCalFormatContrast;
        LMSnew = triLMSNoiseCalFormatContrastnew;
    case 'newConeVar'
        LMSold = triLMSNoiseCalFormat;
        LMSnew = (triLMSNoiseCalFormatContrastnew.*grayLMS)+grayLMS;
end
totalVariance = varianceLMS(varianceType,renderType,LMSold,LMSnew);

end
