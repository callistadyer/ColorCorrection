function [calFormatDiLMS] = DichromatSimulateLinear(calFormatLMSContrast, grayLMS,  constraintWl, cbType, Disp)
% Simulates color vision for dichromatic viewers by projecting onto a plane
% in LMS space
%
% Syntax:
%    [calFormatDiLMS] = DichromatSimulateLinear(calFormatLMS, grayLMS,  constraintWl, cbType, Disp)
%
% Inputs:
%   xyzImage: Input XYZ image
%   cbTypes:  Array of integers specifying the types of color blindness to simulate:
%             1 = Protanopia, 2 = Deuteranopia, 3 = Tritanopia.
%             Example: [1, 2] will simulate Protanopia and Deuteranopia.
%   rgbImage: Input RGB image (optional). If empty, defaults to '74.jpg' on the local path.
%
% Outputs:
%   diXYZ:    dichromat XYZ values
%
% Description:
%   This function renders images to simulate their appearance to
%   color-blind viewers.
%
%   It is a linear algorithm, in that a linear transformation is applied to
%   cone contrast to get the answer.  It operates by projecting cone
%   contrast onto a plane, where the plane contains the achromatic
%   direction and the cone contrast of a specified wavelength.  This
%   ensures that achromatic contrast in comes out as achromatic contrast
%   and also that light at the specified wavelength comes back as itself.
%
%   Basically, this is a linear version of the Brettel, Vienot, Mollon
%   (1997) algorithm. We only specify one contraint that applies to all
%   stimuli, rather than two constraints each of which applies to one half
%   of the color space.
%
%   For some purposes, a linear algorithm has computational advantages
%   (e.g. in enforcing gamut contraints with linear constraints).
%
% History:
%   03/24/2025  cmd, dhb    Wrote it.  Or at least tried to.
%
% Examples:
%{

%}

% Get some key information out of the passed (enhanced) display structure.
wls = Disp.wls;
T_cones = Disp.T_cones;
P_monitor = Disp.P_monitor;

% Convert cone excitations to cone contrast
calFormatLMSContrast = ExcitationsToContrast(calFormatLMS, grayLMS);

% Define achromatic constraint vector
constraint1LMSContrast = [1 1 1]';

% Define monochromatic constraint vector
[~,index] = min(abs(wls-constraintWl)));
index = index(1);
constraint2LMS = T_cones(:,index);
constraint2LMSContrast = ExcitationToContrast(constraint2LMS,grayLMS);

% Now we want to find the best least squares approximation to the trichromatic
% LMS contrast, in the plane spanned by the two constraint vectors
constraintMatrix = [constraint1LMSContrast  constraint2LMSContrast];
calFormatDiLMSContrast = constraintMatrix*(constraintMatrix\calFormatLMSContrast);

% Convert back to excitations
calFormatDiLMS = ContrastToExcitation(calFormatDiLMSContrast,grayLMS);


end

end
