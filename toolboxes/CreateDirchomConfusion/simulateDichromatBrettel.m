
function diLMSCalFormat = simulateDichromatBrettel(lmsImage,renderType,Disp)


error("This function doesn't work properly... (use DichromatSimulateBrettel.m instead)")


% Using Brettel et al. algorithm to try and simulate dichromat rendering
%
% Syntax:
%   function deuterLMSCalFormat = simulateDichromatBrettel(lmsImage,renderType,Disp)
%
% Description:
%
% Inputs:
%   lmsImage:       - LMS representation of the input image
%   renderType:     - Type of dichromacy
%                       'Deuteranopia'
%                       'Protanopia'
%                       'Tritanopia'
%   Disp            - contains spectral sensitivities and wls necessary
%
% Outputs:
%   deuterLMSCalFormat - LMS values for simulating what deuteronope sees
%   when input LMS image is shown
%
% Optional key/value pairs:
%   None
%
% check if in image format (if not, reshape it)
if length(size(lmsImage)) ~= 3
    lmsImage = CalFormatToImage(lmsImage,Disp.m,Disp.n);
end
cone_fundamentals = Disp.T_cones';
wavelengths       = Disp.wls;

% Extract cone sensitivity functions
L_bar = cone_fundamentals(:, 1); % L-cone sensitivity
M_bar = cone_fundamentals(:, 2); % M-cone sensitivity
S_bar = cone_fundamentals(:, 3); % S-cone sensitivity

% Calculate the Neutral Point (E) for Equal-Energy Spectrum
w_white = ones(size(wavelengths)); % Equal-energy spectrum
L_E = trapz(wavelengths, w_white .* L_bar);
M_E = trapz(wavelengths, w_white .* M_bar);
S_E = trapz(wavelengths, w_white .* S_bar);
E = [L_E, M_E, S_E]; % Neutral Point

% Define anchor wavelengths for each dichromacy type
switch renderType
    case 'Protanopia'
        anchor_wavelengths = [575, 475]; % nm
        missing = 'L';
    case 'Deuteranopia'
        anchor_wavelengths = [575, 475]; % nm
        missing = 'M';
    case 'Tritanopia'
        anchor_wavelengths = [660, 485]; % nm
        missing = 'S';
    otherwise
        error('Unknown dichromat type');
end

% Calculate Anchor Points (A1 and A2)
[~, idx_A1] = min(abs(wavelengths - anchor_wavelengths(1))); % Closest to λ1
[~, idx_A2] = min(abs(wavelengths - anchor_wavelengths(2))); % Closest to λ2
A1 = [L_bar(idx_A1), M_bar(idx_A1), S_bar(idx_A1)];
A2 = [L_bar(idx_A2), M_bar(idx_A2), S_bar(idx_A2)];

% Call projection logic with dynamically calculated E, A1, and A2
simulatedLMSimg          = project_to_plane(lmsImage, E, A1, A2, missing);
diLMSCalFormat       = ImageToCalFormat(simulatedLMSimg);
% Note: this arbitrary scaling is necessary because nothing inherent to
% Brettel's algorithm prevents the resulting (projected) RGB values from
% going out of gamut
% [correctedLMS_scaled, k] = scaleInGamut(deuterLMSCalFormat,Disp,0);
% deuterLMSCalFormat       = correctedLMS_scaled;
% Uncomment these to get conversions to RGB images
% RGBImgCalFormatSimulated = LMS2RGBCal     Format(deuterLMSCalFormat,Disp,0);
% RGBImgFormatSimulated    = CalFormatToImage(RGBImgCalFormatSimulated,Disp.m,Disp.n);
end


function Qp = project_to_plane(Q, E, A1, A2, missing)
% Determine which anchor point to use
% Compare the ratio of LMS components for Q with the ratios for A1 and A2

switch missing
    case 'L'
        % Determine projection plane based on S/M ratios
        ratio_Q = Q(:,:,3) ./ Q(:,:,2); % S/M for Q
        ratio_A1 = A1(3) / A1(2); % S/M for A1
        ratio_A2 = A2(3) / A2(2); % S/M for A2

        if abs(ratio_Q - ratio_A1) < abs(ratio_Q - ratio_A2)
            anchor = A1; % Closer to A1
        else
            anchor = A2; % Closer to A2
        end

        % Calculate plane coefficients (a, b, c)
        a = E(2) * anchor(3) - E(3) * anchor(2);
        b = E(3) * anchor(1) - E(1) * anchor(3);
        c = E(1) * anchor(2) - E(2) * anchor(1);
        % Calculate projection onto the chosen plane
        LQp = -(b * Q(:,:,2) + c * Q(:,:,3)) / a;
        Qp = cat(3, Q(:,:,2), LQp, Q(:,:,3));
        % Qp = [LQp; Q(:,:,2); Q(:,:,3)];

    case 'M'
        % Determine projection plane based on S/L ratios
        ratio_Q = Q(:,:,3) ./ Q(:,:,1); % S/L for Q
        ratio_A1 = A1(3) / A1(1); % S/L for A1
        ratio_A2 = A2(3) / A2(1); % S/L for A2

        if abs(ratio_Q - ratio_A1) < abs(ratio_Q - ratio_A2)
            anchor = A1; % Closer to A1
        else
            anchor = A2; % Closer to A2
        end

        % Calculate plane coefficients (a, b, c)
        a = E(2) * anchor(3) - E(3) * anchor(2);
        b = E(3) * anchor(1) - E(1) * anchor(3);
        c = E(1) * anchor(2) - E(2) * anchor(1);
        % Calculate projection onto the chosen plane
        MQp = -(a .* Q(:,:,1) + c .* Q(:,:,3)) / b;
        Qp = cat(3, Q(:,:,1), MQp, Q(:,:,3));
        % Qp = [Q(:,:,1); MQp; Q(:,:,3)];

    case 'S'
        % Determine projection plane based on M/L ratios
        ratio_Q = Q(:,:,2) ./ Q(:,:,1); % M/L for Q
        ratio_A1 = A1(2) / A1(1); % M/L for A1
        ratio_A2 = A2(2) / A2(1); % M/L for A2

        if abs(ratio_Q - ratio_A1) < abs(ratio_Q - ratio_A2)
            anchor = A1; % Closer to A1
        else
            anchor = A2; % Closer to A2
        end
        % Calculate plane coefficients (a, b, c)
        a = E(2) * anchor(3) - E(3) * anchor(2);
        b = E(3) * anchor(1) - E(1) * anchor(3);
        c = E(1) * anchor(2) - E(2) * anchor(1);
        % Calculate projection onto the chosen plane
        SQp = -(a * Q(:,:,1) + b * Q(:,:,2)) / c;
        Qp = cat(3, Q(:,:,1), SQp, Q(:,:,2));
        % Qp = [Q(:,:,1); Q(:,:,2); SQp];
end
end
