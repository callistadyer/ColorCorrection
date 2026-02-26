function outputsGrid = computeInfoDistortionGrid(obj, LMSCalFormat, imgParams, dichromatType, nInfoSteps, nDistSteps, pathName)
% computeInfoDistortionGrid  Build an (info × distortion) grid of optimized color corrections.
%
% Syntax:
%   outputsGrid = computeInfoDistortionGrid(obj, LMSCalFormat, imgParams, dichromatType, nInfoSteps, nDistSteps, pathName)
%
% Description:
% Build a 2-D grid of transformed image where:
%   1) Each row corresponds to a target INFO value
%   2) Each column corresponds to a target DISTORTION value that is achievable while maintaining that row's TARGET INFO
%
%   computeInfoSweep.m minimizes the distortion for each info target.
%   computeInfoDistortionGrid.m gives a range of distortion values at each info target
%
%   For each target info level:
%     (1) Find the min achievable distortion at that info
%     (2) Find the max achievable distortion at that info
%     (3) Interpolate a set of target distortion values between [min, max].
%     (4) For each target distortion, solve for a transform that satisfies
%           info = target info
%           distortion = close to tgt info (in loss function)
% Inputs:
%   obj              - Daltonize/optimizer object. Must provide:
%   LMSCalFormat     - 3×N LMS values of the original image
%   imgParams        - Struct of image-related parameters used by metric functions
%   dichromatType    - type of dichromacy
%                            'Protaniopia'
%                            'Deuteranopia'
%                            'Tritanopia'
%   nInfoSteps       - Number of info targets (rows). Default 10 if []
%   nDistSteps       - Number of distortion targets per row (cols). Default 7 if []
%   pathName         - String describing the where the image is
%                        e.g., 'Deuteranopia/flower2.png/s1_m32_n32'
%
% Outputs:
%   outputsGrid      - nInfoSteps × nDistSteps struct array. Each cell (i,j) has:
%                        .LMSDaltonizedCalFormat            (3×N)   % trichromat LMS, transformed
%                        .rgbLinDaltonizedCalFormat         (3×N)   % trichromat linear RGB, transformed
%                        .LMSDaltonizedRenderedCalFormat    (3×N)   % dichromat-rendered LMS of transformed
%                        .rgbLinDaltonizedRenderedCalFormat (3×N)   % dichromat-rendered linear RGB
%                        .transformRGBmatrix                (3×3)   % optimizing RGB transform for (i,j)
%                        .targetInfoNormalized              (1×1)   % target info for this row i
%                        .targetDistortionNormalized        (1×1)   % target distortion for this col j at row i
%                        .infoNormalized                    (1×1)   % achieved info (normalized)
%                        .distortionNormalized              (1×1)   % achieved distortion (normalized)
%                        .imgParams                        
%                        .Disp                           
%

% Defaults
if isempty(nInfoSteps);  nInfoSteps = 10; end         % number of rows (info targets)
if isempty(nDistSteps);  nDistSteps = 7;  end         % number of cols (distortion targets)

% These epsilons are the half-widths of the "stay near this target" bands.
% epsInfo is small: we want the achieved info to be as close to the row's
% target info as possible.
epsInfo = 1e-6;                                       

% Output paths
projectName = 'ColorCorrection';
outputDir   = getpref(projectName, 'outputDir');
saveBase    = fullfile(outputDir, 'testImagesTransformed');

% Use the exact same metric folder naming as your sweep, so results sit together
metricFolder = buildMetricFolderName(obj.infoFcn, obj.infoParams, obj.distortionFcn);

% Name this "grid run" folder by its size; file will store the whole grid.
runFolder = sprintf('%dinfo_x_%ddist', nInfoSteps, nDistSteps);
saveSubdir = fullfile(saveBase, pathName, metricFolder, runFolder);
saveFile   = fullfile(saveSubdir, 'sweepOutputs_GRID.mat');

% If we've done this already, load 
if exist(saveFile, 'file')
    fprintf('[computeInfoDistortionGrid] Loading old grid from: %s\n', saveFile);
    S = load(saveFile);
    outputsGrid = S.outputsGrid;   % struct array (nInfo x nDist)
    return;
end


% Normalizing
Disp = obj.Disp;
% Convert the original LMS to contrast 
LMSContrastCalFormat = (LMSCalFormat - Disp.grayLMS) ./ Disp.grayLMS;

% Compute how the original image would be seen by three dichromats. Then
% make a composite LMS that uses L from Protan, M from Deutan, S from Tritan
[calFormatLMS_prot, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Protanopia',   Disp);
[calFormatLMS_deut, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Deuteranopia', Disp);
[calFormatLMS_trit, ~, ~] = DichromRenderLinear(LMSCalFormat, 'Tritanopia',   Disp);
LMSCalFormat_new          = [calFormatLMS_prot(1,:); calFormatLMS_deut(2,:); calFormatLMS_trit(3,:)];

% Put that composite into contrast
LMSContrastCalFormat_new  = (LMSCalFormat_new - Disp.grayLMS) ./ Disp.grayLMS;
% Normalizers
normalizerValueToGetRawValue = 1;
infoNormalizer       = obj.infoFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, ...
                           imgParams, dichromatType, normalizerValueToGetRawValue, Disp, obj.infoParams);
distortionNormalizer = obj.distortionFcn(LMSContrastCalFormat, LMSContrastCalFormat_new, ...
                           imgParams, normalizerValueToGetRawValue,                      obj.distortionParams);

% Define the row targets (info values)
% We reproduce your sweep's idea: use lambda endpoints to learn the info range.
%   lambda=0  prioritize "distortion" (i.e., keep image natural) 
%   lambda=1  prioritize "info"       
% Sample linearly in (infoNormalized) between those endpoints
[~,~,~,~,infoNormalized_0,~,~] = colorCorrectionOptimize("lambda", 0, [], [], ...
    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
    infoNormalizer, distortionNormalizer, Disp);

[~,~,~,~,infoNormalized_1,~,~] = colorCorrectionOptimize("lambda", 1, [], [], ...
    LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
    infoNormalizer, distortionNormalizer, Disp);

targetInfoVec = linspace(infoNormalized_0, infoNormalized_1, nInfoSteps);

% Preallocate
LMSDaltonizedCalFormatGrid            = cell(nInfoSteps, nDistSteps);
rgbLinDaltonizedCalFormatGrid         = cell(nInfoSteps, nDistSteps);
LMSDaltonizedRenderedCalFormatGrid    = cell(nInfoSteps, nDistSteps);
rgbLinDaltonizedRenderedCalFormatGrid = cell(nInfoSteps, nDistSteps);
transformRGBmatrixGrid                = cell(nInfoSteps, nDistSteps);
infoNormalizedGrid                    = cell(nInfoSteps, nDistSteps);
distortionNormalizedGrid              = cell(nInfoSteps, nDistSteps);
targetInfoNormalizedGrid              = zeros(nInfoSteps, nDistSteps);
targetDistortionNormalizedGrid        = zeros(nInfoSteps, nDistSteps);

% Prepare gamut constraints for max-distortion
% We will do one extra optimization per row to find the max distortion attainable
% while holding info fixed (that's how we learn the [min,max] distortion range)
[A_total, b_total, ~] = buildGamutConstraints(LMSCalFormat, dichromatType, Disp);
opts = optimoptions('fmincon','Algorithm','interior-point','Display','none', ...
    'ConstraintTolerance',1e-10,'StepTolerance',1e-10,'MaxIterations',200);

%%%%%%%%%%%% Main loop info targets) %%%%%%%%%%%%
for i = 1:nInfoSteps
    thisInfoTarget = targetInfoVec(i);

    % Find the min distortion at this info
    %   objective = minimize distortion (lambda==0)
    %   constraints = sta at target info + gamut constraints
    [~,~,T_min, ~, ~, ~, distN_min] = colorCorrectionOptimize( ...
        'targetinfo', [], thisInfoTarget, [], ...
        LMSCalFormat, imgParams, dichromatType, ...
        obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer,...
        Disp,'T_init', eye(3), 'epsInfo', epsInfo);

    % Find the max distortion at this info 
    % Build the same info constraint as above so info stays near thisInfoTarget
    nonl_info = @(t_vec) infoBandConstraint( ...
        t_vec, LMSCalFormat, imgParams, dichromatType, obj.infoFcn, obj.distortionFcn, ...
        infoNormalizer, distortionNormalizer, Disp, thisInfoTarget, epsInfo);

    % Define a SCALAR objective that returns "normalized distortion"
    % for a given transform vector t_vec
    %
    %   distOnly_norm(...) is a tiny file-scope helper at the bottom of this file.
    %   It calls your lossFunction and returns ONLY the normalized distortion (5th output).
    %
    %   Conceptually:
    %     distOnly(t_vec)  = distortionNormalized(T),  with T reshaped from t_vec
    %     We then solve:  maximize  distOnly(t_vec)
    %                     subject to info band + gamut
    %     which we implement as:  minimize -distOnly(t_vec)
    distOnly = @(t_vec) distOnly_norm( ...
        t_vec, LMSCalFormat, imgParams, dichromatType, ...
        obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp);

    % Maximize distortion by minimizing its negative
    %         Start from T_min(:) — it is already feasible under the info band,
    %         which helps the solver stay in the feasible set from the first step.
    [T_max_vec, ~] = fmincon(@(t) -distOnly(t), T_min(:), A_total, b_total, [], [], [], [], nonl_info, opts);
    distN_max = distOnly(T_max_vec);

    % Make sure [min, max] order is right
    if distN_max < distN_min
        tmp = distN_min; distN_min = distN_max; distN_max = tmp;
    end

    % Build the distortion targets
    distTargets = linspace(distN_min, distN_max, nDistSteps);

    % Optimize with every (info, distortion) combination
    T_seed = eye(3);  
    for j = 1:nDistSteps
        thisDistTarget = distTargets(j);

        [LMS_ij, rgbLin_ij, T_ij, info_ij, infoN_ij, dist_ij, distN_ij] = colorCorrectionOptimize( ...
            'targetdist', [], thisInfoTarget, thisDistTarget, ...
            LMSCalFormat, imgParams, dichromatType, ...
            obj.infoFcn, obj.distortionFcn, infoNormalizer, distortionNormalizer, Disp, ...
            'T_init', T_seed, 'epsInfo', epsInfo);

        % Store the trichromat images and transform
        LMSDaltonizedCalFormatGrid{i,j}    = LMS_ij;
        rgbLinDaltonizedCalFormatGrid{i,j} = rgbLin_ij;
        transformRGBmatrixGrid{i,j}        = T_ij;

        % Store achieved info and distortion (to compare with desired target)
        infoNormalizedGrid{i,j}            = infoN_ij;
        distortionNormalizedGrid{i,j}      = distN_ij;
        targetInfoNormalizedGrid(i,j)      = thisInfoTarget;
        targetDistortionNormalizedGrid(i,j)= thisDistTarget;

        % Dichromat rendering for each
        [LMSDaltonizedRenderedCalFormatGrid{i,j}, rgbLinDaltonizedRenderedCalFormatGrid{i,j}] = ...
            DichromRenderLinear(LMS_ij, dichromatType, Disp);

        % Start next column with current solution
        T_seed = T_ij;
    end
end

% Save the whole grid 
nInfo = size(targetInfoNormalizedGrid,1);
nDist = size(targetDistortionNormalizedGrid,2);

for i = 1:nInfo
    for j = 1:nDist
        outputsGrid(i,j).LMSDaltonizedCalFormat            = LMSDaltonizedCalFormatGrid{i,j};
        outputsGrid(i,j).rgbLinDaltonizedCalFormat         = rgbLinDaltonizedCalFormatGrid{i,j};
        outputsGrid(i,j).LMSDaltonizedRenderedCalFormat    = LMSDaltonizedRenderedCalFormatGrid{i,j};
        outputsGrid(i,j).rgbLinDaltonizedRenderedCalFormat = rgbLinDaltonizedRenderedCalFormatGrid{i,j};
        outputsGrid(i,j).transformRGBmatrix                = transformRGBmatrixGrid{i,j};

        outputsGrid(i,j).targetInfoNormalized              = targetInfoNormalizedGrid(i,1);
        outputsGrid(i,j).targetDistortionNormalized        = targetDistortionNormalizedGrid(i,j);

        outputsGrid(i,j).infoNormalized                    = infoNormalizedGrid{i,j};
        outputsGrid(i,j).distortionNormalized              = distortionNormalizedGrid{i,j};

        outputsGrid(i,j).imgParams                         = imgParams;
        outputsGrid(i,j).Disp                              = Disp;
    end
end

if ~exist(saveSubdir,'dir'); mkdir(saveSubdir); end
save(saveFile, 'outputsGrid', '-v7.3');
return;
end


function dN = distOnly_norm(t_vec, triLMSCalFormat, imgParams, dichromatType, ...
        infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp)
% We need to maximize normalized distortion at a fixed info level to get a
% range of distortion values
% fmincon can only minimize a scalar function, so we define a simple scalar 
% function, dN = distortionNormalized(T), then ask fmincon to minimize the negative
% minimize  -dN under the info constraint (and gamut), which is equivalent to
% maximizing dN while keeping info = target.

[~, ~, ~, ~, dN] = lossFunction('lambda', 0.0, t_vec, ...
    triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
    infoNormalizer, distortionNormalizer, Disp);
end


function [c,ceq] = infoBandConstraint(t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
    infoNormalizer, distortionNormalizer, Disp, targetInfo, epsInfo)

[~, ~, infoNorm_here] = lossFunction('lambda', 0.0, t_vec, ...
    triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
    infoNormalizer, distortionNormalizer, Disp);

c   = (infoNorm_here - targetInfo).^2 - (epsInfo.^2);
ceq = [];
end
