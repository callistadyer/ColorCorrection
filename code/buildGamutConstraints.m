function [A_total, b_total] = buildGamutConstraints(triLMSCalFormat, renderType, Disp)
% buildGamutConstraints
% Constructs linear constraints for color transformation matrix T
% such that both trichromat and dichromat RGB renderings stay within gamut [0,1].
%
% Inputs:
%   triLMSCalFormat   3 x N matrix of LMS values
%   renderType        type of dichromat 
%   Disp              display struct with fields:
%                         M_cones2rgb 
%                         M_rgb2cones
%                         grayRGB   
%                         grayLMS    
%
% Outputs:
%   A_total : Combined constraint matrix (12N x 9)
%   b_total : Combined b vector (12N x 1)

% Number of pixels
nPix = size(triLMSCalFormat, 2);

% Step 1: Convert LMS to RGB
triRGBCalFormat = Disp.M_cones2rgb * triLMSCalFormat;

% Step 2: Convert RGB to contrast RGB
triRGBContrastCalFormat = (triRGBCalFormat - Disp.grayRGB) ./ Disp.grayRGB;

% Step 3: Build trichromat constraint matrix (A_tri)
constraintA = zeros(3 * nPix, 9);
for i = 1:nPix
    r = triRGBContrastCalFormat(1,i);
    g = triRGBContrastCalFormat(2,i);
    b = triRGBContrastCalFormat(3,i);
    block = [eye(3)*r, eye(3)*g, eye(3)*b];
    constraintA((3*i-2):(3*i), :) = block;
end

A_tri_upper = constraintA;
A_tri_lower = -A_tri_upper;
A_tri = [A_tri_upper; A_tri_lower];
b_tri = [ones(nPix * 3, 1); ones(nPix * 3, 1)];

% Step 4: Simulate dichromat LMS and get M_triToDi (conversion from
% trichromat LMS constrat to dichromat LMS contrast)
[~, ~, M_triToDi] = DichromSimulateLinear(triLMSCalFormat, renderType, Disp);

% Step 5: Build M_triRGBc2diRGBc that transforms tri RGB contrast to di RGB contrast
M_conesC2rgbC = diag(1 ./ Disp.grayRGB) * inv(Disp.M_rgb2cones) * diag(Disp.grayLMS);
M_triRGBc2diRGBc = M_conesC2rgbC * M_triToDi * inv(M_conesC2rgbC);

% Step 6: Extend M_all across pixels using Kronecker product
M_triRGBc2diRGBc_allpixels = kron(speye(nPix), M_triRGBc2diRGBc);

% Step 7: Apply M_triRGBc2diRGBc_allpixels to tri RGB contrast constraint matrix
A_di_upper = M_triRGBc2diRGBc_allpixels * constraintA;
A_di_lower = -A_di_upper;
A_di = [A_di_upper; A_di_lower];
b_di = [ones(nPix * 3, 1); ones(nPix * 3, 1)];

% Step 8: Combine both sets of constraints
A_total = [A_tri; A_di];
b_total = [b_tri; b_di];

end
