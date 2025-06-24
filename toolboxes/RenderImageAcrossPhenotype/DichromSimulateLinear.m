function [calFormatDiLMS,M_triToDi] = DichromSimulateLinear(calFormatLMS, grayLMS,  constraintWl, cbType, Disp)
% Simulates color vision for dichromatic viewers by projecting onto a plane
% in LMS space
%
% Syntax:
%    [calFormatDiLMS] = DichromSimulateLinear(calFormatLMS, grayLMS,  constraintWl, cbType, Disp)
%
% Inputs:
%   calFormatLMS:   Input LMS image
%   grayLMS:        gray point in LMS
%   constraintWl:   which wavelength defines the projection plane for
%                   dichromats.
%   cbType:         which type of dichromacy
%                       "Protanopia"
%                       "Deuteranopia"
%                       "Tritanopia"
%   Disp:           display parameters
%
% Outputs:
%   calFormatDiLMS: dichromat LMS values
%   M_triToDi:      transformation matrix that left multiplies on trichromat 
%                   contrast LMS and outputs dichromat contrast LMS values
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

Disp = loadDisplay('ishihara');
[triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues('ishihara','Deuteranopia','M',[],585,3,Disp);
[calFormatDiLMS,M_triToDi] = DichromSimulateLinear(triLMSCalFormat, Disp.grayLMS,  585, "Deuteranopia", Disp)

%}

% Get some key information out of the passed (enhanced) display structure.
wls = Disp.wls;
T_cones = Disp.T_cones;
P_monitor = Disp.P_monitor;

% Convert cone excitations to cone contrast
calFormatLMSContrast = (calFormatLMS - grayLMS)./grayLMS;

% Define monochromatic constraint vector
[~,index] = min(abs(wls-constraintWl));
index = index(1);
constraint2LMS = T_cones(:,index);

% Define achromatic constraint vector
constraint1LMScontrast = [1 1 1]';
constraint2LMScontrast = (constraint2LMS - grayLMS)./grayLMS;

% Define achromatic constraint vector (excitations)
% constraint1LMS = grayLMS;
% constraint2LMS = constraint2LMS;

% Now we want to find the best least squares approximation to the trichromatic
% LMS contrast, in the plane spanned by the two constraint vectors
constraintMatrix = [constraint1LMScontrast  constraint2LMScontrast];

switch (cbType)
    case 'Deuteranopia' % m cone deficiency
        missingConeIdx = 2;
        availableConeIdx = [1 3];
    case 'Protanopia'   % l cone deficiency
        missingConeIdx = 1;
        availableConeIdx = [2 3];
    case 'Tritanopia'   % s cone deficiency
        missingConeIdx = 3;
        availableConeIdx = [1 2];
end

% Calculate missing cone value. See where on the plane the missing cone
% should land when the plane is defined by the two constraint vectors.
A = constraintMatrix(availableConeIdx,:);     % 2x2
B = constraintMatrix(missingConeIdx,:);       % 1x2
% missingCone = B * inv(A) * calFormatLMSContrast(availableConeIdx,:);

% % Replace missing cone row with the value achieved through projection (above)
% calFormatDiLMSContrast = calFormatLMSContrast;
% calFormatDiLMSContrast(missingConeIdx,:) = missingCone;
% % Convert back to excitations
% calFormatDiLMS = (calFormatDiLMSContrast.*grayLMS)+grayLMS;

% Potential way to get M matrix:
M_triToDi = eye(3);
M_triToDi(missingConeIdx, :) = 0;
M_triToDi(missingConeIdx, availableConeIdx) = B * inv(A);

% Compute dichromat from transformation matrix M
% NOTE: M_triToDi operates on contrast LMS 
calFormatDiLMSContast = M_triToDi * calFormatLMSContrast;

% Convert back to LMS excitations
calFormatDiLMS = (calFormatDiLMSContast .* grayLMS) + grayLMS;

calFormatDiRGB = inv(Disp.M_rgb2cones) * calFormatDiLMS;
calFormatDiRGB(calFormatDiRGB>1) = 1;
calFormatDiRGB(calFormatDiRGB<0) = 0;

calFormatDiLMS = Disp.M_rgb2cones * calFormatDiRGB;

end
