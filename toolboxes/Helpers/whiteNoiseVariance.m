
function totalVariance = whiteNoiseVariance(Disp)

% Create a white noise image of size mxn
% m = 256; n = 256;
% rng(2);
s = rng
whiteNoiseImage = rand(Disp.m, Disp.n, 3);
whiteNoiseCalFormat = ImageToCalFormat(whiteNoiseImage);

% MANUAL
% Make hyperspectral img by multiplying primaries * rgb values at each pixel
hyperspectralNoiseCalFormat = Disp.P_monitor * whiteNoiseCalFormat;
% Image format
hyperspectralNoiseImage     = CalFormatToImage(hyperspectralNoiseCalFormat,Disp.m,Disp.n);
hyperspectralNoiseCalFormat = ImageToCalFormat(hyperspectralNoiseImage);


% CREATE SCENE VERSION
scene = sceneFromFile(whiteNoiseImage, 'rgb', [], Disp.d, Disp.wls);
% Get XYZ values (needed for dichromat simulation algorithm from Brettel)
imgXYZ = sceneGet(scene,'xyz');
% Hyperspectral image (needed for calculating cone responses)
hyperspectralImage2 = double(sceneGet(scene,'energy'));
hyperspectralCalFormat2 = ImageToCalFormat(hyperspectralImage2);
triLMSNoiseCalFormat = Disp.T_cones*hyperspectralCalFormat2;

% Compute total variance 
totalVariance = var(triLMSNoiseCalFormat(:));

end