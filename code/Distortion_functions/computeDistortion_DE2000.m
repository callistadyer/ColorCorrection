
function [distortion, distortionNormalized] = computeDistortion_DE2000(LMSCalFormatOld, LMSCalFormatNew, imgParams, normalizingValue, Disp, paramsStruct)
% computeDistortion_squared  Computes sum of squared error between LMS images
%
% Syntax:
%   [distortion, distortionNormalized] = computeDistortion_DE2000(LMS_old, LMS_new, imgParams, normalizingValue, Disp, paramsStruct)
%
% Inputs:
%   LMS_old:            3 x N matrix of original LMS values
%   LMS_new:            3 x N matrix of transformed LMS values
%   imgParams:          image parameters
%   normalizingValue:   Value for normalizing the output
%   paramsStruct: 
%
% Outputs:
%   distortion:  Scalar – distortion in image via how much total LMS has changed
%                         (large value = high distortion)
%
% Example:
%   distortion = ...

% Number of pixels in image
nPixels = imgParams.m * imgParams.n;

% Convert LMS to XYZ (old/original LMS)
XYZ1CalFormat = Disp.M_lms2xyz*LMSCalFormatOld;

% Convert LMS to XYZ (transformed LMS)
XYZ2CalFormat = Disp.M_lms2xyz*LMSCalFormatNew;

% LAB values from XYZ (original LMS)
Lab1CalFormat = XYZToLab(XYZ1CalFormat,Disp.labWhiteXYZ);

% LAB values from XYZ (transformed LMS)
Lab2CalFormat = XYZToLab(XYZ2CalFormat,Disp.labWhiteXYZ);

% Compute deltaE values for color pairs given in CIELAB coordinates
% Essentially a difference or error score.
DE_2000 = ComputeDE2000_Lab(Lab1CalFormat',Lab2CalFormat');

% Calculate the distortion values based on the deltaE results
% DE_2000 is length of number of pixels
distortion = sum(DE_2000.^2)/nPixels;

distortionNormalized = distortion/normalizingValue;

end
