function [inGamut, rgbImageCalFormat] = checkGamut(LMSimageCalFormat, Disp)
% checkGamut  Verifies whether LMS image is within RGB monitor gamut
%
% Syntax:
%   [inGamut, rgbImageCalFormat] = checkGamut(LMSimageCalFormat, Disp)
%
% Inputs:
%   LMSimageCalFormat: 3×N matrix of LMS values in calibration format
%   Disp:              Display struct containing colorimetric transforms
%
% Outputs:
%   inGamut:           Logical flag (1 if all RGB values ∈ [0,1], 0 otherwise)
%   rgbImageCalFormat: 3×N matrix of linear RGB values (same format as input)
%
% Description:
%   This function converts an LMS-calibrated image into the display's
%   linear RGB space and checks whether all resulting RGB values lie
%   within the unit cube [0,1]. Values outside the gamut will result
%   in `inGamut = 0`. Gamma correction is not applied here.
%
% History:
%   06/24/2025  cmd cleaned up
%
% Examples:
%{
Disp = loadDisplay('ishihara');
[triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues('ishihara','Deuteranopia',1,Disp);
[inGamut, rgbOut] = checkGamut(triLMSCalFormat, Disp);
%}

% Convert LMS to RGB
rgbImageCalFormat = LMS2rgbLinCalFormat(LMSimageCalFormat, Disp);
rgbImage = CalFormatToImage(rgbImageCalFormat, Disp.m, Disp.n);

% Tolerance... maybe out of gamut by small amount. should clip this if so
tol = 1e-5;
% Check for out-of-gamut values
if any(rgbImage(:) < -tol | rgbImage(:) > 1+tol)
    inGamut = 0;
else
    inGamut = 1;
end

end
