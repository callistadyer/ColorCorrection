
%% Script for running colorCorrection project code

% Parameters to set. This will determine file location as well
imageTypes = {'gray','ishihara','flower1.png','flower2.png'};
imageType = imageTypes{4};
setType   = 1;
% Set types depend on the image type
%       gray
%           number of squares
%       ishihara
%           plate type
%       .png
%           none so far?
%           maybe could do size of the image. consider downsampling 

dichromatType = 'Deuteranopia';

%%%% you might want to do image size somewhere? maybe another subfolder? 

%% Define output directory.
% The key preference gets set up by the TbTb local hook.
projectName = 'ColorCorrection';
myFullPath = mfilename('fullpath');
[myPath,myName] = fileparts(myFullPath);
% myFullPath = matlab.desktop.editor.getActiveFilename();
% [myPath, myName] = fileparts(myFullPath);

outputDir = getpref(projectName,'outputDir');
if endsWith(imageType, {'.png', '.jpg'}, 'IgnoreCase', true)
    outputSubdir = fullfile(outputDir, 'testImages', dichromatType, imageType);
else
    outputSubdir = fullfile(outputDir, 'testImages', dichromatType, imageType, num2str(setType));
end
if (~exist(outputSubdir,"dir"))
    mkdir(outputSubdir);
end


% Check if image already saved; if so, load it and skip processing
if endsWith(imageType, {'.png', '.jpg'}, 'IgnoreCase', true)
    imageBaseName = imageType;
else
    imageBaseName = [imageType, '.png'];
end

imageOutputPath = fullfile(outputSubdir, imageBaseName);

if exist(imageOutputPath, 'file')
    fprintf('Found precomputed image: %s\n', imageOutputPath);
    triRGBImage = im2double(imread(imageOutputPath));

    % Optionally, load other saved .mat files
    triLMSPath = fullfile(outputSubdir, 'triLMSCalFormat.mat');
    triRGBPath = fullfile(outputSubdir, 'triRGBCalFormat.mat');
    dispPath   = fullfile(outputSubdir, 'Disp.mat');
    
    if exist(triLMSPath, 'file'), load(triLMSPath, 'triLMSCalFormat'); end
    if exist(triRGBPath, 'file'), load(triRGBPath, 'triRGBCalFormat'); end
    if exist(dispPath,   'file'), load(dispPath, 'Disp'); end
else
    % Type of image
    img         = imageType;
    % Type of dichromat
    renderType  = dichromatType;
    % Make function that determines set variables
    [setParams] = buildSetParameters(img,setType);

    % Load display
    Disp = loadDisplay(img);

    % Load LMS values for this image
    [triLMSCalFormat,triRGBCalFormat,Disp] = loadLMSvalues(img,renderType,setParams,Disp);
    triRGBImage = CalFormatToImage(triRGBCalFormat,Disp.m,Disp.n);

    % Save out the values
    save(fullfile(outputSubdir, 'triLMSCalFormat.mat'), 'triLMSCalFormat');
    save(fullfile(outputSubdir, 'triRGBCalFormat.mat'), 'triRGBCalFormat');
    save(fullfile(outputSubdir, 'Disp.mat'), 'Disp');

    % Save out the image
    imwrite(triRGBImage, imageOutputPath);
end
