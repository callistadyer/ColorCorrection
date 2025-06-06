
function totalSimilarity = whiteNoiseSimilarity(similarityType,Disp)
% Calculates total similarity metric to normalize squared error similarity of given image size 
%
% Syntax:
%   totalSimilarity = whiteNoiseSimilarity(similarityType,Disp)
%
% Inputs:
%   similarityType:  Type of similarity computation
%                     "angle"
%                     "squared"
%   Disp:          Display parameters
% Outputs:
%   totalSimilarity: Similarity normalizer term 
%
% Optional key/value pairs:
%   None
%

% Random seed to calculate total variance 
rng(2);
% Create a white noise image of size mxn
whiteNoiseImage = rand(Disp.m, Disp.n, 3);
whiteNoiseCalFormat = ImageToCalFormat(whiteNoiseImage);

rng(3);
whiteNoiseImage2 = rand(Disp.m, Disp.n, 3);
whiteNoiseCalFormat2 = ImageToCalFormat(whiteNoiseImage2);


% 1... MANUAL
% Make hyperspectral img by multiplying primaries * rgb values at each pixel
hyperspectralNoiseCalFormat  = Disp.P_monitor * whiteNoiseCalFormat;
hyperspectralNoiseCalFormat2 = Disp.P_monitor * whiteNoiseCalFormat2;


% LMS image
triLMSNoiseCalFormat       = Disp.T_cones*hyperspectralNoiseCalFormat;
triLMSNoiseCalFormat2      = Disp.T_cones*hyperspectralNoiseCalFormat2;

% % Uncomment this to normalize for RGB values if you want to use RGB in
% % similarity metric
% triLMSNoiseCalFormat= whiteNoiseCalFormat;
% triLMSNoiseCalFormat2= whiteNoiseCalFormat2;

% Compute similarity 
totalSimilarity = similarityLMS(similarityType,triLMSNoiseCalFormat,triLMSNoiseCalFormat2);

end
