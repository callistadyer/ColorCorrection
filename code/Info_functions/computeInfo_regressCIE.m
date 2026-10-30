function [info,infoNormalized,LMSCalFormat_pred,beta] = computeInfo_regressCIE(LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, dichromatType, normalizingValue, Disp, paramsStruct)
% computeInfo_regressCIE  How much of the original color can be worked out from what the dichromat sees
%
% Syntax:
%   [info,infoNormalized,LMSCalFormat_pred,beta] = computeInfo_regressCIE( ...
%       LMSContrastCalFormat_old, LMSContrastCalFormat_new, imgParams, ...
%       dichromatType, normalizingValue, Disp, paramsStruct)
%
% Description:
%   A linear regression predicts something about the ORIGINAL image (the
%   target) from the TRANSFORMED image's cones that the dichromat has (the
%   predictors). The better the prediction, the higher the info.
%
%   paramsStruct.predictingWhat (target) and how it is scored:
%     'M'     : original M. Rebuilt image = original L, predicted M,
%               original S. Info = -DE2000(original, rebuilt image).
%     'L,M,S' : original L, M and S. Rebuilt image = predicted L, M, S.
%               Info = -DE2000(original, rebuilt image).
%     'L-M'   : original L-M (the red-green signal). Rebuilt image =
%               original L, predicted M (= original L - predicted L-M),
%               original S. Info = -DE2000(original, rebuilt image).
%               Rewards transforms that make L-M stand out in what the
%               dichromat sees.
%     ('L' and 'S' work like 'M'.)
%
%   paramsStruct.predictingFromWhat (predictors), e.g.
%     'L and S'  : new L, new S (cone contrast)
%     'opponent' : what the dichromat PERCEIVES: the transformed image as
%                  the dichromat sees it, in CIELAB [L*, a*, b*] (DE units).
%                  Set weightPenaltySD in DE units (e.g. 1), so differences
%                  below about 1 DE count as not visible.
%
%   paramsStruct.whichWeightPenalty: WHICH method stops the regression
%   from using huge weights to blow up a tiny signal.
%     'none'  : plain regression (default)
%     'noise' : add noise (SD = weightPenaltySD) to the predictors
%     'ridge' : ridge regression with the same strength as that noise
%   paramsStruct.weightPenaltySD: HOW STRONG that method is, in cone
%   contrast (default 0.02). Signals smaller than about this count as
%   unreadable.
%
% Inputs:
%   LMSContrastCalFormat_old   3 x N original LMS contrast
%   LMSContrastCalFormat_new   3 x N transformed LMS contrast
%   imgParams                  image parameters
%   dichromatType              'Protanopia', 'Deuteranopia', or 'Tritanopia'
%                              (only used for defaults if paramsStruct is empty)
%   normalizingValue           infoNormalized = info / normalizingValue
%   Disp                       display parameters
%   paramsStruct               settings above
%
% Outputs:
%   info                       info score (higher = better)
%   infoNormalized             info / normalizingValue
%   LMSCalFormat_pred          3 x N rebuilt ("predicted") LMS image. For
%                              'L-M': original L, predicted M, original S
%                              (used for the DE2000 score and the picture).
%   beta                       regression weights, rows = predictors,
%                              columns = targets. Huge weights = a tiny
%                              signal is being blown up.


%% 1. Default target and predictors (if none given): predict the missing cone
if nargin < 7 || isempty(paramsStruct) || ...
   ~isfield(paramsStruct,'predictingWhat') || isempty(paramsStruct.predictingWhat) || ...
   ~isfield(paramsStruct,'predictingFromWhat') || isempty(paramsStruct.predictingFromWhat)

    switch dichromatType
        case 'Protanopia'    % Missing L; predict L from [M S]
            paramsStruct.predictingWhat     = 'L';
            paramsStruct.predictingFromWhat = 'M and S';
        case 'Deuteranopia'  % Missing M; predict M from [L S]
            paramsStruct.predictingWhat     = 'M';
            paramsStruct.predictingFromWhat = 'L and S';
        case 'Tritanopia'    % Missing S; predict S from [L M]
            paramsStruct.predictingWhat     = 'S';
            paramsStruct.predictingFromWhat = 'L and M';
        otherwise
            error('Unknown DichromatType: %s. Use ''Protanopia'', ''Deuteranopia'', or ''Tritanopia''.', dichromatType);
    end
end

% Weight penalty: two settings, both set in the runner
%   whichWeightPenalty = WHICH method ('none', 'noise', 'ridge')
%   weightPenaltySD    = HOW STRONG it is (e.g. 0.02)
% For each: use the runner's value if given, otherwise the default.
whichWeightPenalty = 'none';                                   % default: no penalty
if isfield(paramsStruct, 'whichWeightPenalty')                 % set in the runner?
    whichWeightPenalty = paramsStruct.whichWeightPenalty;      % yes -> use it
end
weightPenaltySD = 0.02;                                        % default strength
if isfield(paramsStruct, 'weightPenaltySD')                    % set in the runner?
    weightPenaltySD = paramsStruct.weightPenaltySD;            % yes -> use it
end


%% 2. Target (y) and predictors (X). Rows = pixels.
L_old = LMSContrastCalFormat_old(1,:);  M_old = LMSContrastCalFormat_old(2,:);  S_old = LMSContrastCalFormat_old(3,:);
L_new = LMSContrastCalFormat_new(1,:);  M_new = LMSContrastCalFormat_new(2,:);  S_new = LMSContrastCalFormat_new(3,:);

y = buildTarget(paramsStruct.predictingWhat, L_old, M_old, S_old);                      % from the ORIGINAL
X = buildPredictors(paramsStruct.predictingFromWhat, L_old, M_old, S_old, L_new, M_new, S_new, ...
                    LMSContrastCalFormat_new, dichromatType, Disp);                       % from the TRANSFORMED

if isrow(y)
    y = y(:);                       % one target: make it a column (N x 1)
end
N = size(X,1);                      % number of pixels


%% 3. Clean up the predictors
% Set tiny contrasts to exactly 0
smallestInterestingContrast = 1e-5;
X(abs(X) < smallestInterestingContrast) = 0;

% 'noise' penalty: add the same noise on every call (rng(1)), so fmincon
% sees a fixed function
if strcmp(whichWeightPenalty, 'noise')
    rng(1);
    X = X + weightPenaltySD * randn(size(X));
end

% 'L-M' only: remove the averages. Only the PATTERN of L-M across the
% image matters, not its overall level (same as adding a constant term).
yMean = 0;
if strcmp(paramsStruct.predictingWhat, 'L-M')
    yMean = mean(y);                % added back in step 5, for the picture
    y     = y - yMean;
    X     = X - mean(X, 1);
end


%% 4. Regression: find weights (beta) so that X * beta is as close to y as possible
if all(X == 0, 'all')
    % No signal at all in the predictors: nothing can be predicted
    beta = zeros(size(X,2), size(y,2));
else
    switch whichWeightPenalty
        case {'none', 'noise'}
            beta = X \ y;                                              % least squares
        case 'ridge'
            lambdaRidge = N * weightPenaltySD^2;                       % same strength as the noise
            beta = (X'*X + lambdaRidge*eye(size(X,2))) \ (X'*y);       % ridge
        otherwise
            error('Unsupported whichWeightPenalty: %s. Use ''none'', ''noise'', or ''ridge''.', whichWeightPenalty);
    end
end

y_hat    = X * beta;                % prediction
residual = y - y_hat;               % the part of the target that could NOT be predicted


%% 5. Rebuilt ("predicted") image, in LMS contrast
% Start from the original; replace what was predicted
LMSContrast_pred = LMSContrastCalFormat_old;

switch paramsStruct.predictingWhat
    case 'L'
        LMSContrast_pred(1,:) = y_hat(:).';
    case 'M'
        LMSContrast_pred(2,:) = y_hat(:).';
    case 'S'
        LMSContrast_pred(3,:) = y_hat(:).';
    case 'L,M,S'
        LMSContrast_pred = y_hat';                       % y_hat is N x 3 here
    case 'L-M'
        % The regression predicted only L-M, one number per pixel. To make
        % an image we need L, M and S, so:
        % (Sizes: N = number of pixels.)

        % Step 1: predicted L-M (add back the average removed in step 3)
        %   y_hat: [N x 1]
        %   .'   : [1 x N]
        %   yMean: [1 x 1]
        predictedLminusM = y_hat(:).' + yMean;          % [1 x N]

        % Step 2: assume the dichromat has the original L, so keep it
        predictedL = L_old;                             % [1 x N]

        % Step 3: solve for M using the predicted L-M:
        %       predicted L-M = original L - predicted M
        %       predicted M   = original L - predicted L-M
        predictedM = predictedL - predictedLminusM;     % 1 x N  =  1 x N - 1 x N

        % Step 4: predicted image = [original L, predicted M, original S]
        LMSContrast_pred = [predictedL; predictedM; S_old];   % [3 x N]
    otherwise
        error('Unsupported predictingWhat: %s', paramsStruct.predictingWhat);
end

% LMS contrast -> LMS excitations
LMSCalFormat_old  = (LMSContrastCalFormat_old .* Disp.grayLMS) + Disp.grayLMS;
LMSCalFormat_pred = (LMSContrast_pred         .* Disp.grayLMS) + Disp.grayLMS;


%% 6. Info score (higher = better)
% All cases (including 'L-M'): -(perceptual error, DE2000) between the
% original and the predicted image.
% For 'L-M' the predicted image is [original L, predicted M, original S],
% so it differs from the original only where L-M was not predicted.
[distortion, ~] = computeDistortion_DE2000(LMSCalFormat_old, LMSCalFormat_pred, imgParams, 1, Disp, paramsStruct);
info = -distortion;

% R^2 VERSION FOR 'L-M' (turned off for now; scores with DE2000 above instead).
% To turn back on, uncomment, and uncomment the matching normalizer lines
% in computeSweep.m.
% if strcmp(paramsStruct.predictingWhat, 'L-M')
%     % R^2 = 1 - (L-M left over) / (all of the L-M pattern)
%     %   1 = all of L-M predicted, 0 = none
%     info = 1 - sum(residual.^2) / max(sum(y.^2), eps);   % eps: avoid dividing by 0
% end


%% 7. Normalize
infoNormalized = info / normalizingValue;

end


function y = buildTarget(spec, L, M, S)
% Target from the ORIGINAL image
spec = strtrim(spec);
switch spec
    case 'L',     y = L;
    case 'M',     y = M;
    case 'S',     y = S;
    case 'L-M',   y = L - M;
    case 'L,M,S', y = [L(:), M(:), S(:)];
    otherwise
        error('Unsupported predictingWhat: %s', spec);
end
end

function X = buildPredictors(spec, L_old, M_old, S_old, L_new, M_new, S_new, LMSContrastCalFormat_new, dichromatType, Disp)
% Predictors from the TRANSFORMED image. Returns N x P (rows = pixels).
spec = strtrim(spec);
switch spec
    case 'opponent'
        % What the dichromat PERCEIVES, in CIELAB units (1 unit ~ 1 DE):
        %   L* = lightness, a* = red-green, b* = yellow-blue
        % Use with weightPenaltySD in DE units (e.g. 1): differences
        % smaller than about that count as not visible.
        % Step 1: transformed image, LMS contrast -> LMS excitations     3 x N
        newLMS = LMSContrastCalFormat_new .* Disp.grayLMS + Disp.grayLMS;
        % Step 2: how the dichromat sees it                               3 x N
        diLMS  = DichromRenderLinear(newLMS, dichromatType, Disp);
        % Step 3: LMS -> XYZ -> CIELAB                                     3 x N
        diLab  = XYZToLab(Disp.M_lms2xyz * diLMS, Disp.labWhiteXYZ);
        % Step 4: one row per pixel: [L*, a*, b*]                          N x 3
        X = diLab';
    case 'L and S'
        X = [L_new(:), S_new(:)];
    case 'L and M'
        X = [L_new(:), M_new(:)];
    case 'M and S'
        X = [M_new(:), S_new(:)];
    case 'L,M,S'
        X = [L_new(:), M_new(:), S_new(:)];
    case 'deltaL and deltaS'
        X = [(L_new(:)-L_old(:)), (S_new(:)-S_old(:))];
    case 'deltaL+M and deltaS'
        X = [(L_new(:)+M_new(:) - (L_old(:)+M_old(:))), (S_new(:)-S_old(:))];
    otherwise
        error('Unsupported predictingFromWhat: %s', spec);
end
end
