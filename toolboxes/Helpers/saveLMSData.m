function saveLMSData(outputSubdir, triLMSCalFormat, trirgbLinCalFormat, triRGBImage, diLMSCalFormat, dirgbLinCalFormat, diRGBImage, Disp, imgParams)
% saveLMSData  Saves LMS and RGB data for trichromat and dichromat simulations
%
% Syntax:
%   saveLMSData(outputSubdir, triLMSCalFormat, trirgbLinCalFormat, triRGBImage, ...
%               diLMSCalFormat, dirgbLinCalFormat, diRGBImage, Disp, imgParams)
%
% Inputs:
%   outputSubdir:         Directory where output .mat and .png files should be saved
%   triLMSCalFormat:      LMS image (calibration format) for trichromatic observer
%   trirgbLinCalFormat:   Linear RGB image (calibration format) for trichromatic observer
%   triRGBImage:          Gamma-corrected RGB image for trichromat (for visualization)
%   diLMSCalFormat:       LMS image (calibration format) for dichromatic observer
%   dirgbLinCalFormat:    Linear RGB image (calibration format) for dichromatic observer
%   diRGBImage:           Gamma-corrected RGB image for dichromat (for visualization)
%   Disp:                 Display structure (e.g., color conversion matrices)
%   imgParams:            Image parameters used during generation or rendering
%
% Description:
%   This function saves all relevant LMS and RGB data used in the color vision simulation
%   pipeline. It includes both trichromatic and dichromatic versions of the image in
%   calibration and gamma-corrected formats, as well as the associated display and image
%   parameter structures.
%
% History:
%   04/14/2025  cmd  Initial draft and documentation update

% Construct full paths to output files 
% These files will be saved into the specified subdirectory
triLMSPath    = fullfile(outputSubdir, 'triLMSCalFormat.mat');       % LMS values for trichromat
trirgbLinPath = fullfile(outputSubdir, 'trirgbLinCalFormat.mat');    % Linear RGB for trichromat
triRGBPath    = fullfile(outputSubdir, 'triRGBCalFormat.mat');       % Gamma-corrected RGB for trichromat

diLMSPath     = fullfile(outputSubdir, 'diLMSCalFormat.mat');        % LMS values for dichromat
dirgbLinPath  = fullfile(outputSubdir, 'dirgbLinCalFormat.mat');     % Linear RGB for dichromat
diRGBPath     = fullfile(outputSubdir, 'diRGBCalFormat.mat');        % Gamma-corrected RGB for dichromat

dispPath      = fullfile(outputSubdir, 'Disp.mat');                  % Display structure
imgParamsPath = fullfile(outputSubdir, 'imgParams.mat');             % Parameters used for image generation

% Save trichromat data
save(triLMSPath,    'triLMSCalFormat');       % Save trichromat LMS in calibration format
save(trirgbLinPath, 'trirgbLinCalFormat');    % Save trichromat linear RGB
save(triRGBPath,    'triRGBImage');           % Save trichromat gamma-corrected RGB image
imwrite(triRGBImage, fullfile(outputSubdir, 'trichromat_render.png'));  % Save image file for easy viewing

% Save dichromat data 
save(diLMSPath,     'diLMSCalFormat');        % Save dichromat LMS in calibration format
save(dirgbLinPath,  'dirgbLinCalFormat');     % Save dichromat linear RGB
save(diRGBPath,     'diRGBImage');            % Save dichromat gamma-corrected RGB image
imwrite(diRGBImage, fullfile(outputSubdir, 'dichromat_render.png'));    % Save image file for easy viewing

% Save display and rendering parameters 
save(dispPath,      'Disp');                  % Save display calibration information
save(imgParamsPath, 'imgParams');             % Save image generation parameters

% Print a confirmation message to the console
fprintf('[saveLMSData] Saved LMS and RGB data to: %s\n', outputSubdir);

end
