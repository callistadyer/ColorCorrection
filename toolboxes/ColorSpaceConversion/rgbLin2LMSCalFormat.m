function lmsCalFormat = rgbLin2LMSCalFormat(rgbLinCalFormat,Disp,scaleFactor)

% function takes rgb in cal format and outputs LMS image in cal format 
%
% Syntax:
%   lmsCalFormat = rgbLin2LMSCalFormat(rgbLinCalFormat,Disp,scaleFactor,bScale)
%
% Description:
%
% Inputs:
%   rgbLinCalFormat       - linear rgb image
%   Disp includes:
%       Disp.T_cones               - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor             - [nWlx3]. Display primaries
%       Disp.m                     - Scalar.  Row dimension of image     
%       Disp.n                     - Scalar.  Column dimension of image  
%   scaleFactor           - scale factor used to scale img values   
%
% Outputs:
%   lmsCalFormat          - LMS image
%
% Optional key/value pairs:
%   None
%

% LMS image
M_rgb2cones = Disp.T_cones*Disp.P_monitor;
lmsCalFormat = M_rgb2cones * rgbLinCalFormat;

end