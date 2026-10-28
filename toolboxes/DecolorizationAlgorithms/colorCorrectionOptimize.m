function [LMSDaltonizedCalFormat, rgbLinDaltonizedCalFormat, transformRGBmatrix, info, infoNormalized, distortion, distortionNormalized, nlconSatisfied, fminconInfo] = ...
    colorCorrectionOptimize(target, lambda, targetInfo, targetDist, ...
    triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, infoParams, distortionParams, infoNormalizer, distortionNormalizer, Disp, varargin)
% Optimizes a linear transformation to enhance color contrast for dichromats
%
% % Syntax:
% [triLMSCalFormatOpt, trirgbLinCalFormatOpt, transformRGBmatrix, info, infoNormalized, distortion, distortionNormalized] = ...
%     colorCorrectionOptimize(useLambdaOrTargetInfo, lambdaOrTargetInfo, ...
%     triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp, varargin)
%
% Inputs:
%   target                     Are you looking for a specific 'lambda','info', or 'distortion'?
%   lambda                     Tradeoff weight in [0,1] when target='lambda'.
%                                   lambda=1 means maximize information
%                                   lambda=0 means minimize distortion
%                              
%   targetInfo                 Desired normalized information when target='info'.
%   targetDist                 Desired normalized distortion value when target='distortion'.
%
%   triLMSCalFormat            Trichromatic LMS data in CalFormat
%   imgParams                  Image/metric params used by infoFnc and distortionFcn.
%   dichromatType              Dichromat type
%                                   'Deuteranopia'
%                                   'Protanopia'
%                                   'Tritanopia'.
%   infoFnc                    Handle for the information metric.
%   distortionFcn              Handle for similarity metric to preserve naturalness.
%                                   'squared' (sum of squared error),
%   infoParams
%   distortionParams
%   infoNormalizer             Scalar to normalize the info metric.
%   distortionNormalizer       Scalar to normalize the distortion metric.
%   Disp:                      Struct with display parameters. Must include:
%                                   .M_cones2rgb  – Matrix to convert LMS to RGB
%                                   .T_cones      – Cone fundamentals
%                                   .P_monitor    – Monitor spectral power distribution
%                                   .wls          – Wavelength sampling
%
% Outputs:
%   LMSDaltonizedCalFormat:    daltonized LMS rendering for trichromat
%   rgbLinDaltonizedCalFormat: daltonized rgb Linear rendering for trichromat
%   transformRGBmatrix:        transformation matrix
%   info:                      info from infoFcn
%   infoNormalized:            info from infoFcn, normalized
%   distortion:                distortion from distortionFcn
%   distortionNormalized:      distortion from distortionFcn, normalized
%   nlconSatisfied:            true if the returned transform satisfies the band constraint
%   fminconInfo:               struct describing how fmincon stopped:
%                                   .exitflag   final fmincon exitflag (NaN if skipFmincon)
%                                   .reason     exitflag in words
%                                   .iterations total iterations, including restarts
%                                   .funcCount  total function evaluations, including restarts
%                                   .nRestarts  warm restarts after hitting MaxIterations
%
% Constraints:
%   - RGB values must be between 0 and 1
%
% Optional key/value pairs:
%   'T_init'                   3x3 starting transformation (default eye(3)).
%   'skipFmincon'              If true, return T_init without optimizing (default false).
%   'bandTol'                  Half-width of the allowed band around targetInfo
%                              or targetDist. Required when target='info' or
%                              'distortion'. computeSweep sets it (1% of the
%                              endpoint range) and chooseStartPoints passes the
%                              same value to its violation check, so there is
%                              one band definition.
%
% Examples are included within the code

% History
%   11/03?/2025  cmd  Initial go.
%
% Examples:
%{
%}

% Sample data to see how code workswe
if isempty(triLMSCalFormat)
    % Pretend data to visualize what the function does
    N = 100;                         % Number of points
    x = randn(1, N);                 % Random data for the first dimension
    y = 2 * x + randn(1, N) * 0.2;   % Strong correlation with the first dimension
    z = 0.5 * x + 0.3 * randn(1, N); % Weaker correlation with the first dimension
    triLMSCalFormat = [x; y; z];
end

%% Key value pairs
parser = inputParser;
parser.addParameter('T_init', eye(3), @(x) isnumeric(x) && isequal(size(x),[3 3]));

parser.addParameter('skipFmincon', false, @(x) islogical(x) && isscalar(x));

% Half-width of the band around the target (no default: the caller must set it)
parser.addParameter('bandTol', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));

parser.parse(varargin{:});

T_init      = parser.Results.T_init;
skipFmincon = parser.Results.skipFmincon;
bandTol     = parser.Results.bandTol;

switch lower(target)
    case 'lambda'
        if isempty(lambda), error('target ''lambda'' requires a lambda.'); end
    case 'info'
        if isempty(targetInfo), error('target ''info'' requires target info.'); end
        if isempty(bandTol), error('target ''info'' requires ''bandTol''.'); end
    case 'distortion'
        if isempty(targetDist)
            error('target ''distortion'' requires target distortion');
        end
        if isempty(bandTol), error('target ''distortion'' requires ''bandTol''.'); end
    otherwise
        error('target must be ''lambda'', ''targetinfo'', or ''targetdist''.');
end
 
rng(1);

% Get linear constraints
[A_total, b_total, M_triRGBc2diRGBc] = buildGamutConstraints(triLMSCalFormat, dichromatType, Disp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Just reached optimization')

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...  
    'ConstraintTolerance',  1e-10, ...  
    'StepTolerance',        1e-10, ...
    'FiniteDifferenceType', 'central', ...   % see note below options
    'Display','iter', ...
    'MaxIterations',200);

% FiniteDifferenceType = 'central' (fmincon's default is 'forward'):
% central differences estimate gradients more accurately, especially near
% the solution where the true slope is small. Tested 10/28 on Gaugin, 100x100,
% L,M,S info sweep vs the same run with forward differences:
%   - fmincon runs that stalled (exitflag 2): 35 -> 9; converged: 14 -> 42;
%     failed / no feasible point: 2 -> 0
%   - a little slower: +13% run time (8700 s vs 7679 s), despite +84%
%     function evaluations
%   - did NOT change the results: identical transforms and images at every
%     step (max pixel difference 1e-5), because another start point usually
%     converged anyway in the forward-difference run
% Kept for now for more reliable individual fmincon runs. Revisit if run
% time becomes a problem.

% paramsStruct = struct();  % default empty
if exist('infoParams','var') && ~isempty(infoParams)
    paramsStruct = infoParams;

    if isfield(infoParams,'predictingWhat')
        paramsStruct.predictingWhat = infoParams.predictingWhat;
    end
    if isfield(infoParams,'predictingFromWhat')
        paramsStruct.predictingFromWhat = infoParams.predictingFromWhat;
    end
end


% Attempt to add in nonlinear constraint where we keep info a constant
% value and search over distortion values (want to keep info as the
% target info and then get the transformation that minimizes distortion, or
% keep distortion fixed and maximize info)
switch lower(target)

    case 'lambda'
        % Trade off info vs distortion via lambda
        fun     = @(t_vec) lossFunction('lambda', lambda, t_vec, ...
                        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp,paramsStruct);
        nonlcon = [];

    case 'info'
        % Minimize distortion subject to info constraint (=targetInfo)

        % Minimize distortion only (lambda=0)
        fun = @(t_vec) lossFunction('lambda', 0.0, t_vec, ...
                        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp,paramsStruct);
        % Nonlinear constraint on info
        nonlcon = @(t_vec) constraintBand( ...
                        t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, paramsStruct, Disp, ...
                        'info', targetInfo, bandTol);

    case 'distortion'
        % Maximize info subject to distortion constraint (=targetDist)

        % Maximize info only (lambda=1)
        fun = @(t_vec) lossFunction('lambda', 1.0, t_vec, ...
                        triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, Disp,paramsStruct);
        % Nonlinear constraint on distortion
        nonlcon = @(t_vec) constraintBand( ...
                        t_vec, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
                        infoNormalizer, distortionNormalizer, paramsStruct, Disp, ...
                        'distortion', targetDist, bandTol);
end

% Now do the minimization with that 
% [transformRGB_opt, fval, exitflag, output] = fmincon(fun, T_init(:), ...            
%     A_total, b_total, ...           % Linear inequality constraints (gamut)
%     [], [], [], [], ...             % Aeq, beq, lb, ub 
%     nonlcon, ...                    % Nonlinear constraint
%     options);

if skipFmincon
    % This is to make sure that the first step (in the forward pass, last
    % step in backwards pass) always returns identity. Basically, don't
    % even give the algorithm a chance to deviate from identity when it is
    % the first step.
    transformRGB_opt = T_init(:);   % return identity (or whatever T_init was)

    % fmincon never runs here, so there is no fval/exitflag/output from it.
    % Fill them with placeholders so later code can still refer to them.
    fval = NaN; exitflag = NaN; output = struct();
    nlconSatisfied = true;   % forced endpoint; do not test nonlinear band

    % Record that fmincon was skipped (exitflag NaN, zero iterations), so
    % anyone reading the saved candidate files can tell this step apart
    % from a step where fmincon actually ran.
    fminconInfo = struct('exitflag', NaN, 'reason', 'skipped (forced T_init)', ...
        'iterations', 0, 'funcCount', 0, 'nRestarts', 0);
else
    % Run fmincon from the starting transformation T_init.
    %   fun      - objective (distortion for an info sweep, -info for a
    %              distortion sweep, or the lambda tradeoff)
    %   A_total, b_total - linear gamut constraints (RGB must stay in [0,1])
    %   nonlcon  - band constraint keeping info (or distortion) within
    %              bandTol of the target ([] for a lambda search)
    % Outputs:
    %   transformRGB_opt - the 3x3 transformation found, as a 9x1 vector
    %   exitflag         - why fmincon stopped (see exitflagReason below)
    %   output           - details, including .iterations and .funcCount
    [transformRGB_opt, fval, exitflag, output] = fmincon(fun, T_init(:), ...
        A_total, b_total, ...
        [], [], [], [], ...
        nonlcon, ...
        options);

    % Keep running totals of the work done, so that if fmincon is
    % restarted below, the counts cover all of the runs combined.
    totalIterations = output.iterations;   % optimizer steps taken
    totalFuncCount  = output.funcCount;    % objective evaluations (includes the
                                           % extra ones for finite-difference gradients)

    % WARM RESTARTS
    % exitflag == 0 means fmincon hit MaxIterations (or MaxFunctionEvaluations)
    % before converging: it simply ran out of iterations, it did not finish.
    % In that case, start fmincon again from the point where it stopped
    % (transformRGB_opt), not from T_init. Allow at most maxRestarts extra
    % runs, so a search that never settles cannot loop forever.
    maxRestarts = 3;
    nRestarts   = 0;
    while exitflag == 0 && nRestarts < maxRestarts
        nRestarts = nRestarts + 1;
        fprintf('fmincon hit its iteration limit; warm restart %d of %d\n', nRestarts, maxRestarts);
        [transformRGB_opt, fval, exitflag, output] = fmincon(fun, transformRGB_opt(:), ...
            A_total, b_total, ...
            [], [], [], [], ...
            nonlcon, ...
            options);
        % Add this run's work to the running totals
        totalIterations = totalIterations + output.iterations;
        totalFuncCount  = totalFuncCount  + output.funcCount;
    end

    % Summarize how fmincon finished. exitflag here is from the LAST run
    % (after any restarts). This struct is returned as an output and
    % chooseStartPoints saves it in each candidate's .mat file.
    fminconInfo = struct('exitflag', exitflag, 'reason', exitflagReason(exitflag), ...
        'iterations', totalIterations, 'funcCount', totalFuncCount, 'nRestarts', nRestarts);

    % Print a one-line summary for the command window log
    fprintf('fmincon exitflag %d (%s), %d iterations, %d restarts\n', ...
        exitflag, fminconInfo.reason, totalIterations, nRestarts);

    % exitflag == -2: fmincon could not find ANY transformation that meets
    % all the constraints (gamut + band). The returned transformRGB_opt is
    % then not a valid answer, so warn loudly. Note the band check below
    % only re-tests the band (nonlcon), not the gamut constraints, so it
    % will not catch an answer that is infeasible only because of gamut.
    if exitflag == -2
        warning('colorCorrectionOptimize:noFeasiblePoint', ...
            'fmincon found no feasible point (exitflag -2).');
    end
end


if skipFmincon
    % forced endpoint
    nlconSatisfied = true;

elseif isempty(nonlcon)
    % no nonlinear constraint to satisfy
    nlconSatisfied = true;

% Fallback if nonlinear constraint violated
else

    % We want cchk to be less than 0, that would mean that the solution is
    % within the constraint band
    [cchk, ~] = nonlcon(transformRGB_opt(:));

    % If maxC <= 0, all inequality constraints are satisfied exactly
    % If maxC > 0, that value is the magnitude of the worst violation
    maxC = max(cchk(:));

    % Determine whether the nonlinear constraints are satisfied:
    %   All inequality constraints must be <= tol
    %   All equality constraints must have |ceq| <= tol
    % nlconSatisfied = (maxC <= tol) && (maxAbsCeq <= tol);
    % slack = 1e-12;  % tiny numerical wiggle room
    slack = options.ConstraintTolerance;
    nlconSatisfied = (maxC <= slack);

    % If the nonlinear constraints are not satisfied:
    %   Reject the optimized solution
    %   Fall back to "Start T"
    if ~nlconSatisfied
        transformRGB_opt = T_init(:);   % For the identity start, this will resort to using identity, for other it will resort to T_prev
        fprintf('Constraint violated: maxC=%.3g (<=0 is feasible). Returning T_init.\n', maxC);
    end
end

disp('Just finished optimization')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the info and distortion values for the transformation matrix that
% ultimately was chosen after the optimization
[~, info, infoNormalized, distortion, distortionNormalized,LMSCalFormat_pred] = lossFunction('lambda', 0.0, transformRGB_opt, triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, Disp,paramsStruct);


% rgbCalFormat_pred = LMS2rgbLinCalFormat(LMSCalFormat_pred,Disp);
% rgbCalFormat_pred(rgbCalFormat_pred<0) = 0;
% rgbCalFormat_pred(rgbCalFormat_pred>1) = 1;
% 
% RGB_pred = rgbLin2RGB(rgbCalFormat_pred,Disp);
% RGBimg_pred = CalFormatToImage(RGB_pred,imgParams.m,imgParams.n);
% 
% figure();
% imagesc(RGBimg_pred); axis square

% the idea here was that maybe the constraint for distortion worked but
% keeping it in gamut didn't? or maybe the attempt to keep it in gamut is
% what's making it fail the search? 

% Reshape optimal solution into matrix
transformRGBmatrix = reshape(transformRGB_opt, 3, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEP 7: TRANSFORM CONTRAST IMAGE W OPTIMAL %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert LMS to RGB
triRGBCalFormat = Disp.M_cones2rgb * triLMSCalFormat;

% Convert RGB to contrast RGB
triRGBContrastCalFormat = (triRGBCalFormat - Disp.grayRGB) ./ Disp.grayRGB;

% Transform RGB contrast image
triRGBContrastCalFormat_T = transformRGBmatrix * triRGBContrastCalFormat;
diRGBContrastCalFormat_T  = M_triRGBc2diRGBc * triRGBContrastCalFormat_T;

% Add back in gray before outputting the image
diRGBCalFormatOpt = (diRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;
rgbLinDaltonizedCalFormat = (triRGBContrastCalFormat_T.*Disp.grayRGB) + Disp.grayRGB;

% Cut off values outside of gamut, when there is some weird numerical out
% of bounds 
if (max(rgbLinDaltonizedCalFormat(:))>1) && (max(rgbLinDaltonizedCalFormat(:)) <1+1e-2)
    rgbLinDaltonizedCalFormat(rgbLinDaltonizedCalFormat>1)=1;
end

if (min(rgbLinDaltonizedCalFormat(:))<0)  && (min(rgbLinDaltonizedCalFormat(:))>0-1e-2)
    rgbLinDaltonizedCalFormat(rgbLinDaltonizedCalFormat<0)=0;
end

triRGBCalFormatOpt = rgbLin2RGB(rgbLinDaltonizedCalFormat,Disp);

RGBimg = CalFormatToImage(triRGBCalFormatOpt,imgParams.m,imgParams.n);

% figure();
% imagesc(RGBimg); axis square

% Get LMS values to output
LMSDaltonizedCalFormat = Disp.M_rgb2cones * rgbLinDaltonizedCalFormat;

end

function reason = exitflagReason(exitflag)
% Short description of an fmincon exitflag
switch exitflag
    case 1,  reason = 'converged: first-order optimality met';
    case 2,  reason = 'stalled: step smaller than StepTolerance';
    case 0,  reason = 'hit MaxIterations or MaxFunctionEvaluations';
    case -1, reason = 'stopped by output or plot function';
    case -2, reason = 'no feasible point found';
    case -3, reason = 'objective unbounded below';
    otherwise, reason = sprintf('exitflag %d', exitflag);
end
end

function [c, ceq] = constraintBand(t_vec, triLMSCalFormat, imgParams, dichromatType,infoFnc, distortionFcn, infoNormalizer, distortionNormalizer, paramsStruct, Disp, whichMetric, targetVal, epsBand)

% Doesn't really matter what the loss is here, I just want to grab the
% info and distortion norm vals from the loss function
[~, ~, infoNorm, ~, distNorm] = lossFunction('lambda', 0.0, t_vec, ...
    triLMSCalFormat, imgParams, dichromatType, infoFnc, distortionFcn, ...
    infoNormalizer, distortionNormalizer, Disp,paramsStruct);

switch lower(whichMetric)
    case {'info','information'}
        metricVal = infoNorm;
    case {'dist','distortion'}
        metricVal = distNorm;
    otherwise
        error('constraintBand: whichMetric must be ''info'' or ''distortion''.');
end

% Try to keep the distortion or info value as close as possible to the
% target value (epsBand gives some wiggle room)
% Old band inequality: (metric - target)^2 - eps^2 <= 0
% Its gradient 2*(metric-target)*grad(metric) is ~0 near the target, which
% gives fmincon (interior-point) almost no guidance.
% c   = (metricVal - targetVal).^2 - (epsBand.^2);

% Same band as two one-sided constraints, each linear in the metric:
%   metric <= target + eps   and   metric >= target - eps
c   = [metricVal - targetVal - epsBand; ...
       targetVal - metricVal - epsBand];

ceq = [];
end
