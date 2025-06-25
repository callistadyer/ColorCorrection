
% Script for running colorCorrection project code

% Define output directory.
% The key preference gets set up by the TbTb local hook.
projectName = 'ColorCorrection';
myFullPath = mfilename('fullpath');
[myPath,myName] = fileparts(myFullPath);
% myFullPath = matlab.desktop.editor.getActiveFilename();
% [myPath, myName] = fileparts(myFullPath);


imageType = 'gray';
setType   = 'gray1';
dichromatType = 'Deuteranopia';

outputDir = getpref(projectName,'outputDir');
outputSubdir = fullfile(outputDir,'testImages',imageType,setType,dichromatType);
if (~exist(outputSubdir,"dir"))
    mkdir(outputSubdir);
end


%% Loading in images. Want this to be completely separate from correction routine

img = imageType;
renderType = dichromatType;

% make function that determines set variables
[setParams] = buildSetParameters(img,setType);

% Load display 
Disp = loadDisplay(img);

% Load LMS values for this image
[triLMSCalFormat,diLMSCalFormat,Disp] = loadLMSvalues(img,renderType,setParams,Disp);

save(fullfile(outputSubdir, 'triLMSCalFormat.mat'), 'triLMSCalFormat');


%%