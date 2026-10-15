function diLMSCalFormat = DichromatSimulateBrettelLMS(LMSCalFormat, dichromatType, Disp, m, n)
% DichromatSimulateBrettelLMS
% Syntax:
%   diLMSCalFormat = DichromatSimulateBrettelLMS(LMSCalFormat, dichromatType, Disp, m, n)
%
% Inputs:
%   LMSCalFormat   - LMS CalFormat 3 x N 
%   dichromatType  - Color blindness type:
%                       'Protanopia' 
%                       'Deuteranopia'
%                       'Tritanopia'
%   Disp           - Contains
%                     - Disp.M_lms2xyz
%                     - Disp.labWhiteXYZ
%                     - Disp.M_rgb2xyz
%                     - Disp.M_rgb2cones
%   m,n            - dimensions of image
%
% Outputs:
%   diLMSCalFormat  - 3 x N dichromat LMS 

switch lower(string(dichromatType))
    case "protanopia"
        cbType = 1;
    case "deuteranopia"
        cbType = 2;
    case "tritanopia"
        cbType = 3;
    otherwise
        error('Unknown dichromatType: %s (use Protanopia, Deuteranopia, or Tritanopia)', string(dichromatType));
end

% Grab everything from Disp 
whiteXYZ  = Disp.labWhiteXYZ(:);
M_lms2xyz = Disp.M_lms2xyz;
M_xyz2lms = Disp.M_rgb2cones * inv(Disp.M_rgb2xyz);

M_xyz2lms_stock = colorTransformMatrix('xyz2lms');     % XYZ to Stockman LMS
M_lms2xyz_stock = colorTransformMatrix('lms2xyz');     % Stockman LMS to XYZ


% LMS to XYZ 
xyzCalFormat = M_lms2xyz * LMSCalFormat;

% xyz2lms expects image shaped input
xyzImage = CalFormatToImage(xyzCalFormat, m, n);

% Brettel XYZ image to dichromat LMS image 
% The Brettel algorithm is implemented inside xyz2lms, which I think outputs LMS
% values in the Stockman LMS coordinate system (colorTransformMatrix('xyz2lms')).
% I think we must convert to xyz, then go from xyz back to our LMS
%%%%%% This is where Brettel is implemented!! %%%%%%
lmsDiImage_stock = xyz2lms(xyzImage, cbType, whiteXYZ); % is it ok to use our whiteXYZ here?

% Stockman LMS to XYZ 
xyzDiImage = imageLinearTransform(lmsDiImage_stock, M_lms2xyz_stock);

% XYZ to our LMS
xyzDiCalFormat = ImageToCalFormat(xyzDiImage);            
diLMSCalFormat = M_xyz2lms_disp * xyzDiCalFormat;         
end
