function [A_total, b_total, M_triRGBc2diRGBc] = buildGamutConstraints(triLMSCalFormat, renderType, Disp)
% buildGamutConstraints
% Constructs linear constraints for color transformation matrix T
% such that both trichromat and dichromat RGB renderings stay within gamut [0,1].
%
% Description:
%   T (3x3, passed to fmincon as a 9x1 vector) multiplies each pixel's RGB
%   contrast vector c: new contrast = T*c. Staying in gamut means every
%   channel of the new RGB stays in [0,1], which in contrast units is
%   -1 <= (T*c)_k <= 1 (and the same for the dichromat rendering, which is a
%   fixed matrix times T*c). Each of these is linear in T, so they can be
%   written as A_total * T(:) <= b_total and handed to fmincon.
%
%   Each pixel gives 12 rows: 3 channels x (upper, lower) x (trichromat,
%   dichromat).
%
% CONVEX-HULL REDUCTION (added 10/28)
%   Building rows for every pixel gave 12*N rows (120,000 at 100x100), which
%   made each fmincon run slow and made fmincon's 'sqp' algorithm unusable
%   (it needed a 107 GB matrix).
%
%   Only a few pixels can ever be the limiting ones. For a fixed T, every
%   row is a linear function of one pixel's contrast vector c (for example
%   (T*c)_1 <= 1). The largest value of a linear function over a set of
%   points is always reached at a corner of the set's convex hull (the
%   smallest convex shape containing all the points). So if the hull
%   pixels satisfy a row, every other pixel does too. Keeping only the hull
%   pixels' rows therefore gives EXACTLY the same allowed region for T, not
%   an approximation. The dichromat rows are also linear in c (a fixed
%   matrix times T*c), so the same hull pixels cover them.
%
%   Tested 10/29 on a full sweep (Gaugin 100x100, L,M,S): 8.4x faster, but
%   steps 3-6 got worse answers (distortion 1.3-2.5x higher). The allowed
%   region is the same; fmincon's interior-point method just takes a
%   different path with fewer constraints and lands in a worse local minimum.
%
%   Set useConvexHull = false below to build rows for every pixel (the old
%   behavior), e.g. for comparisons.
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
%   A_total : Combined constraint matrix (12*nUsed x 9), where nUsed is the
%             number of convex-hull pixels (or N if the hull is not used)
%   b_total : Combined b vector (12*nUsed x 1)
%   M_triRGBc2diRGBc: matrix the converts trichromat contrast RGB to
%                     dichromat contrast RGB

% Use only convex-hull pixels (same constraints, far fewer rows). See above.
% true for updates 2 and 3; false for updates 0 and 1.
% CURRENTLY true: update 2 run, flower1 100x100 (10/29).
useConvexHull = true;

% Step 1: Convert LMS to RGB
triRGBCalFormat = Disp.M_cones2rgb * triLMSCalFormat;

% Step 2: Convert RGB to contrast RGB
triRGBContrastCalFormat = (triRGBCalFormat - Disp.grayRGB) ./ Disp.grayRGB;

% Step 2b: Keep only the pixels that can be limiting: the corners of the
% convex hull of the pixels' RGB contrast vectors (see CONVEX-HULL REDUCTION).
% convhulln needs at least 4 points that are not all on one plane (a gray or
% two-color image is flat); in those cases, fall back to using every pixel.
nPixAll = size(triRGBContrastCalFormat, 2);
usePix  = 1:nPixAll;
if useConvexHull && nPixAll >= 4
    try
        hullFacets = convhulln(triRGBContrastCalFormat');   % rows = triangles of pixel indices
        usePix     = unique(hullFacets(:))';                % pixels at the hull's corners
    catch
        % Degenerate point set (e.g. all colors on a plane): keep all pixels
        usePix = 1:nPixAll;
    end
end
triRGBContrastCalFormat = triRGBContrastCalFormat(:, usePix);

% Number of pixels whose constraints are built
nPix = size(triRGBContrastCalFormat, 2);

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
% trichromat LMS constrat to dichromat LMS contrast).
% (Uses the full image; M_triToDi depends only on the dichromat type and
% display, not on which pixels are kept.)
[~, ~, M_triToDi] = DichromRenderLinear(triLMSCalFormat, renderType, Disp);

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
