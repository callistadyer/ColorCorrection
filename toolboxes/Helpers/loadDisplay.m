function Disp = loadDisplay()

% loadDisplay  Sets display parameters and color space transforms
%
% Syntax:
%   Disp = loadDisplay()
%
% Inputs:
%   none.
%
% Outputs:
%   Disp: Struct with display parameters including:
%         - Disp.wls                 : Wavelength sampling (400–700 nm)
%         - Disp.d                   : Display structure (from ISET)
%         - Disp.P_monitor           : Monitor spectral power distribution
%         - Disp.T_cones             : Cone fundamentals (SS2 2°)
%         - Disp.M_rgb2cones         : Matrix from RGB to LMS
%         - Disp.M_cones2rgb         : Matrix from LMS to RGB
%         - Disp.grayRGB             : RGB vector for neutral gray [0.5 0.5 0.5]
%         - Disp.grayLMS             : LMS value corresponding to gray
%
% Description:
%   Loads and defines display colorimetric transformations using a predefined
%   LCD display model (Apple LCD from ISETBio). Converts spectral monitor
%   primaries to cone fundamentals to compute RGB ↔ LMS mappings. 
%   Image size is inferred based on the type or path provided.
%
% History:
%   06/17/2025  cmd commented 
%   07/01/2025  cmd updated. made it picture independent
%
% Examples:
%{
Disp = loadDisplay()
%}

% Display parameters

% Wavelengths to be used for display
wls = (400:10:700)';
% Display
d = displayCreate('LCD-Apple');
% Monitor primaries
P_monitor = SplineSrf(displayGet(d, 'wave'), displayGet(d, 'spd'), wls);
% Spectral sensitivities (cone fundamentals)
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wls);
% Conversions from rgb to cones 
M_rgb2cones = T_cones*P_monitor;

% Conversions from rgb to xyz
load T_xyz1931.mat
T_xyz = SplineCmf(S_xyz1931,T_xyz1931,wls);
M_rgb2xyz = T_xyz * P_monitor;
M_lms2xyz = M_rgb2xyz * inv(M_rgb2cones);

labWhiteXYZ = M_rgb2xyz * [1 1 1]';

% Put parameters into Disp structure
Disp.wls       = wls;
Disp.d         = d;
Disp.P_monitor = P_monitor;
Disp.T_cones   = T_cones;

Disp.M_rgb2cones = Disp.T_cones*Disp.P_monitor;
Disp.M_cones2rgb = inv(M_rgb2cones);
Disp.M_rgb2xyz   = M_rgb2xyz;
Disp.M_lms2xyz   = M_lms2xyz;

Disp.labWhiteXYZ = labWhiteXYZ;
Disp.grayRGB     = [0.5 0.5 0.5]';
Disp.grayLMS     = M_rgb2cones*Disp.grayRGB;
end
