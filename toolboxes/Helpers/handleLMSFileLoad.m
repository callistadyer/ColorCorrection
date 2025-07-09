function [didLoad, triLMSCalFormat, trirgbLinCalFormat, triRGBImage, diLMSCalFormat, dirgbLinCalFormat, diRGBImage, pathName, outputSubdir] = handleLMSFileLoad(img, renderType, imgParams, Disp, clearTestImages)
% handleLMSFileLoad  Handles file path setup, directory clearing, loading, and cache validation
%
% Syntax:
%   [didLoad, triLMSCalFormat, trirgbLinCalFormat, triRGBImage, ...
%    diLMSCalFormat, dirgbLinCalFormat, diRGBImage, ...
%    pathName, outputSubdir] = handleLMSFileLoad(img, renderType, imgParams, Disp, clearTestImages)
%
% Inputs:
%   img:               Image identifier or filename 
%                           (e.g., 'ishihara' or 'flower.png')
%   renderType:        Type of dichromacy simulation    
%                           ('Deuteranopia', 'Protanopia', 'Tritanopia')
%   imgParams:         Structure containing parameters such as setType, m, n
%   Disp:              Display structure
%   clearTestImages:   whether to delete the testImages folder and regenerate
%
% Outputs:
%   didLoad:           True if cached data was successfully loaded and matched `imgParams`
%   triLMSCalFormat:   Trichromat LMS (calibration format)
%   trirgbLinCalFormat: Trichromat linear RGB (calibration format)
%   triRGBImage:       Trichromat gamma-corrected image
%   diLMSCalFormat:    Dichromat LMS (calibration format)
%   dirgbLinCalFormat: Dichromat linear RGB (calibration format)
%   diRGBImage:        Dichromat gamma-corrected image
%   pathName:          Truncated subfolder path for organizational use
%   outputSubdir:      Full path to output directory
%
% Description:
%   This function handles output directory naming, folder clearing, and cache checking.
%   If all expected files exist and match the input `imgParams`, it loads and returns them.
%   Otherwise, it returns `didLoad = false`, signaling that regeneration is needed.

% Base output path
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');
testImagesDir = fullfile(outputDir, 'testImages');

%  Clear all contents in 'testImages' folder if requested
if clearTestImages
    if exist(testImagesDir, 'dir')
        fprintf('Clearing entire folder: %s\n', testImagesDir);
        rmdir(testImagesDir, 's');   % Recursively delete all files and subfolders
    end
    mkdir(testImagesDir);            % Recreate the folder after deletion
end

%  Build unique parameter-based subfolder name (e.g., 's1_m128_n128') 
paramFields = {'setType', 'm', 'n'};
paramSubdirParts = {};
for i = 1:length(paramFields)
    field = paramFields{i};
    if isfield(imgParams, field)
        prefix = lower(field(1));  % e.g., 's', 'm', 'n'
        paramSubdirParts{end+1} = [prefix num2str(imgParams.(field))]; %#ok<AGROW>
    end
end
paramSubdir = strjoin(paramSubdirParts, '_');  % Join parts with underscores

%  Determine full output subdirectory path 
outputSubdir = fullfile(outputDir, 'testImages', renderType, img, paramSubdir);
if ~exist(outputSubdir, "dir")
    mkdir(outputSubdir);  % Create folder if it doesn't already exist
end

%  Extract path relative to 'testImages' 
idx = strfind(outputSubdir, ['testImages' filesep]);
pathName = outputSubdir(idx + length('testImages') + 1 : end);  % Remove prefix to shorten

%  Define expected file paths 
triLMSPath    = fullfile(outputSubdir, 'triLMSCalFormat.mat');
trirgbLinPath = fullfile(outputSubdir, 'trirgbLinCalFormat.mat');
triRGBPath    = fullfile(outputSubdir, 'triRGBCalFormat.mat');

diLMSPath     = fullfile(outputSubdir, 'diLMSCalFormat.mat');
dirgbLinPath  = fullfile(outputSubdir, 'dirgbLinCalFormat.mat');
diRGBPath     = fullfile(outputSubdir, 'diRGBCalFormat.mat');

dispPath      = fullfile(outputSubdir, 'Disp.mat');
imgParamsPath = fullfile(outputSubdir, 'imgParams.mat');

%  Check if all expected files exist 
if all(cellfun(@exist, ...
    {triLMSPath, trirgbLinPath, triRGBPath, ...
     diLMSPath, dirgbLinPath, diRGBPath, ...
     dispPath, imgParamsPath}, ...
     repmat({'file'}, 1, 8)))

    % Load and compare imgParams
    tmp = load(imgParamsPath, 'imgParams');
    if isequaln(tmp.imgParams, imgParams)
        % All files found and parameters match — load all
        load(triLMSPath,    'triLMSCalFormat');
        load(trirgbLinPath, 'trirgbLinCalFormat');
        load(triRGBPath,    'triRGBImage');

        load(diLMSPath,     'diLMSCalFormat');
        load(dirgbLinPath,  'dirgbLinCalFormat');
        load(diRGBPath,     'diRGBImage');

        load(dispPath,      'Disp');

        didLoad = true;
        return;
    else
        fprintf('[handleLMSFileLoad] Cached imgParams mismatch – regenerating data.\n');
    end
end

%  If missing files or mismatched parameters, return empty 
didLoad = false;
triLMSCalFormat   = [];
trirgbLinCalFormat = [];
triRGBImage       = [];

diLMSCalFormat    = [];
dirgbLinCalFormat = [];
diRGBImage        = [];
end
