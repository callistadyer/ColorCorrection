function lmsCalFormat = rgbLin2LMSCalFormat(rgbLinCalFormat,Disp)

% function takes linear rgb in cal format and outputs LMS image in cal format 
%
% Syntax:
%   lmsCalFormat = rgbLin2LMSCalFormat(rgbLinCalFormat,Disp,scaleFactor,bScale)
%
% Description:
%
% Inputs:
%   rgbLinCalFormat       - linear rgb image
%   Disp includes:
%       Disp.T_cones      - [3xnWl]. Cone spectral sensitivities
%       Disp.P_monitor    - [nWlx3]. Display primaries
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