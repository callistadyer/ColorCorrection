# CLAUDE.md

Context for working in this repository.

## What this project is

MATLAB project that evaluates image transformations designed to improve the
color information available to dichromatic observers while minimizing
distortion of the original image. Transformations are found by constrained
optimization: minimize distortion subject to a target information score
('info' sweep) or maximize information subject to a target distortion level
('distortion' sweep). Top-level driver scripts are in `code/` (e.g.
`code/runAllTestImages.m`, `code/runAllInfoTypes.m`); core classes and
helpers are in `toolboxes/` (e.g. the `@daltonize` class,
`toolboxes/Helpers/`); `tutorials/` and `validation/` are what they sound
like. `reconstruction/` holds a separate reconstruction subproject
(Kadkhodaie / Zhang code) that is not part of the main pipeline.

## Setup: ToolboxToolbox (TbTb)

This repo is managed with ToolboxToolbox. `tbUseProject('ColorCorrection')`
reads `configuration/ColorCorrection.json`, deploys dependencies, sets the
path, and runs the local hook (deployed from
`configuration/ColorCorrectionLocalHookTemplate.m` into the user's local
hooks folder as `ColorCorrectionLocalHook.m` — edit the local copy, never
the template, for machine-specific config).

Full resolved dependency tree (traced through the local ToolboxRegistry,
current as of 2026-10):

- Psychtoolbox-3 — https://github.com/Psychtoolbox-3/Psychtoolbox-3
  (the `Psychtoolbox` subfolder goes on the path)
- ISETBioCSFGenerator — https://github.com/isetbio/ISETBioCSFGenerator
  - mQUESTPlus — https://github.com/brainardlab/mQUESTPlus
  - isetbio — https://github.com/isetbio/isetbio
    - isetcam — https://github.com/iset/isetcam

## MATLAB preferences and data layout

The only preference currently used is `outputDir` in the `ColorCorrection`
preference group (set by the local hook; on lab machines it resolves to
`<Dropbox for Business>/DALT_analysis/ColorCorrection`). Important:
`outputDir` is not output-only —

- Input test images (`flower1.png`, `fruit.png`, `Gaugin.png`,
  `ishi45.png`, ...) are read from `<outputDir>/ColorCorrectionImages/`
  (see `toolboxes/Helpers/loadLMSvalues.m`, ~line 100). These images are
  deliberately NOT in the repo (code/data separation). Only the
  programmatically generated `'ishihara'` images work without them.
- `<outputDir>/testImages/...` is a cache: `handleLMSFileLoad.m` loads
  previously generated .mat files back if the saved `imgParams` match.
- `<outputDir>/testImagesTransformed/...` holds sweep results written by
  `computeSweep` and read back by `tutorials/t_applyTransforms.m`.

- Finished sweeps are labeled by "update" (a numbered set of optimizer
  settings): `.../<metric folder>/info_10steps_update<N>`. New runs write
  to plain `info_10steps` and are renamed to `_update<N>` once finished.
  `updateLog.txt` (repo root) records what each update number means, its
  settings table, which runs used it, and results. `runAllTestImages.m`
  has `updateNumber`: `[]` runs a new sweep; a number loads that update's
  finished sweep (look-up only; errors if missing, never computes).

So don't treat `outputDir` as disposable between pipeline stages, and
don't suggest committing the images to the repo — the plan is to post
them on an open science repository.

## Git notes

- Remote: git@github.com:callistadyer/ColorCorrection.git. Default branch
  is `main`; local `main` tracks `origin/main` (fixed 2026-10-15 — it
  previously tracked `origin/main-old`).
- `origin/main-old` is a leftover stale branch the owner chose to leave
  as-is on GitHub. Don't delete it, don't push to it.
- README.md documents dependencies and the preference setup for users who
  don't know TbTb; keep it in sync when prefs or dependencies change.

## TODO

- [ ] Post the input test images (the `ColorCorrectionImages` folder
      contents) to an open science repository so users can grab them.
      Until then they are available from the authors upon request (the
      README says so). When posted, update the README with the download
      location and instructions.
- [ ] At some point, run a code review of the codebase when asked
      (the user mentioned they may request this).
- [ ] Work through the optimizer-robustness TODO list below (from a
      2026-10-20 read-through of the search code; DHB will ask Callista
      to tackle these with Claude's help). In progress: items 1, 3, 5, 9
      done; see the STATUS lines and "Current state and findings" below.

## TODO: improving the parameter search (fmincon robustness)

Context: the sweeps often get stuck in local minima or terminate before
reaching the loss-minimizing transform. These suggestions came from reading
`colorCorrectionOptimize.m`, `lossFunction.m`, `computeSweep.m`,
`chooseStartPoints.m`, `findStartPoints.m`, `buildGamutConstraints.m`, and
the info/distortion metrics. Ordered by expected impact; the first three are
a few lines each. After each change, rerun a known-problematic case and
compare achieved loss / constraint satisfaction before vs. after (note that
sweep results and per-candidate fmincon solutions are cached under
`<outputDir>/testImagesTransformed/...` — delete or redirect the relevant
cache folders when re-testing, or stale results will be loaded).

1. [x] Reformulate the band constraint. In
       `colorCorrectionOptimize.m` (`constraintBand`, ~line 318) the
       constraint is `(metric - target)^2 - eps^2 <= 0`. Its gradient
       `2*(metric-target)*grad(metric)` vanishes at the band center and is
       only ~2*eps near the boundary — nearly useless for interior-point.
       Replace with two one-sided constraints, linear in the metric:
       `c = [metric - target - eps; target - metric - eps]`. Same feasible
       set, well-scaled gradients.
       STATUS (10/22): done in `constraintBand`. Every sweep step now lands
       exactly on the band edge.

2. [ ] Remove `fminconFactor = 1e8` in `lossFunction.m` (~line 77). The
       metrics are already normalized to O(1); scaling by 1e8 makes the
       default absolute `OptimalityTolerance` (1e-6) unreachable, so runs
       stop on StepTolerance/MaxIterations instead of optimality, and it
       skews the interior-point merit function. Set it to 1 and tune
       tolerances deliberately if needed.
       STATUS (10/28): tested, kept 1e8. Factor 1 (40x40 flower1 + Gaugin)
       cut stalls 22 -> 1 but caused 10 "no feasible point" failures, ran
       12% slower, and gave the same final answers. Comment in lossFunction.m.

3. [x] Stop ignoring `exitflag`. In `colorCorrectionOptimize.m` (~line 190)
       `exitflag`/`output` are discarded. Log them into the saved per-
       candidate .mat files (chooseStartPoints caching), and when
       exitflag == 0 (hit MaxIterations = 200) warm-restart fmincon from
       the returned point (e.g. `while exitflag == 0 && nRestarts < 3`).
       exitflag == -2 (no feasible point) should also be surfaced.
       STATUS (10/22): done. colorCorrectionOptimize returns `fminconInfo`
       (exitflag, reason, iterations, funcCount, restarts), warm-restarts up
       to 3 times on exitflag 0, warns on -2; chooseStartPoints saves
       `fminconInfo` in each candidate .mat file.

4. [ ] Loosen the post-hoc feasibility slack. `colorCorrectionOptimize.m`
       (~line 221) accepts solutions only if maxC <= 1e-12, stricter than
       fmincon's own ConstraintTolerance (1e-10), so legitimately feasible
       solutions get rejected and replaced by T_init, discarding the whole
       run. Use the solver's ConstraintTolerance, and on violation return
       the optimized T with its violation recorded (chooseStartPoints
       already has a smallest-violation fallback that can use it) instead
       of reverting to T_init.
       STATUS: first half done (slack = options.ConstraintTolerance; with the
       linear band from #1, 1e-12 rejected a valid solution in testing).
       Second half (return the optimized T on failure) NOT done.

5. [x] Fix the band-width inconsistency for negative info targets. The
       optimizer uses `pctTol * abs(targetInfo)`
       (`colorCorrectionOptimize.m` ~line 147) but the violation check in
       `chooseStartPoints.m` (~line 172) uses `pctTol * thisX` with no
       abs(). With `computeInfo_regressCIE`, normalized info is negative,
       so the checker's band collapses to absTolFloor while the optimizer
       worked to a much wider band — feasible solutions get classified as
       violations. Factor the band computation into one shared helper used
       by both. (Also note the two call sites default different
       pctTol/absTolFloor values; unify.)
       STATUS (10/28): done. pctTol/absTolFloor removed; computeSweep sets
       one absolute band (1% of the endpoint info range) and passes it as
       `bandTol` to both the optimizer and the violation check.

6. [ ] Address nonsmoothness (interior-point + finite differences assume
       smooth functions):
       - Remove the hard threshold `X(abs(X) < 1e-5) = 0` in
         `computeInfo_regressCIE.m` / `computeInfo_regressSquared.m`
         (~line 54) or replace with smooth shrinkage.
       - CIEDE2000 (default distortion, and inside regressCIE) has
         abs()/hue-branch discontinuities. Diagnostic: rerun a problem
         case with `computeDistortion_squared` +
         `computeInfo_regressSquared` (both smooth; distortion exactly
         quadratic in T). If that converges cleanly, DE2000 nonsmoothness
         is the culprit.
       - One-liner either way: add `'FiniteDifferenceType','central'` to
         the fmincon options.
       STATUS (10/28): central differences adopted (update 1). Gaugin
       100x100 L,M,S: stalls 35 -> 9, converged 14 -> 42, failures 2 -> 0,
       identical images, +13% run time. Threshold and DE2000 not changed.
       (computeInfo_regressSquared.m, the smooth diagnostic metric, was
       deleted 10/22; recover from git history if needed.)

7. [ ] (More work, biggest robustness win) Supply analytic gradients via
       `'SpecifyObjectiveGradient', true`. Squared distortion is exactly
       quadratic in the 9 transform parameters; the regression info has a
       closed-form gradient by the envelope theorem (at the least-squares
       beta, only dX/dT contracted with the residual is needed).
       Eliminates finite-difference noise (a common cause of StepTolerance
       stalls) and cuts function evaluations ~10x.
       STATUS: not done; lower priority. Central differences already gave
       the reliability gain, and +84% function evaluations cost only +13%
       run time, so evaluations are not the main cost.

8. [ ] Add random restarts to the existing multistart. The candidate-start
       machinery in `chooseStartPoints.m` is deterministic (identity,
       previous step, info-sweep seeds, presolve seeds). Add 3-5 jittered
       starts per step (`T_candidate + 0.05*randn(3)`) to the
       candNames/candTinit lists. The `rng(1)` in
       `colorCorrectionOptimize.m` (~line 103) currently does nothing
       since nothing downstream is stochastic. MATLAB's MultiStart (Global
       Optimization Toolbox) is a drop-in alternative.
       STATUS: not done; now the main open problem (see Findings below).
       A neighbour-blend start ((T(i-1)+T(i+1))/2) was tested on flower1 M
       100x100 step 9 and did NOT help.

9. [x] Prune gamut constraints to the convex hull.
       `buildGamutConstraints.m` builds 12N x 9 rows (10,800 for a 30x30
       image), but only extreme colors can bind. `convhulln` of the 3D
       RGB-contrast point cloud (and separately of its dichromat-mapped
       image) gives an exactly equivalent constraint set, typically tens
       of rows — faster iterations, better conditioning, and budget for
       more restarts.
       STATUS (10/29): done (update 2; `useConvexHull` in
       buildGamutConstraints.m). One hull of the RGB-contrast vectors covers
       the dichromat rows too (they are linear in the same vectors).
       100x100: 120,000 -> ~1,000-1,500 rows; verified identical constraint
       values. Sweep 8.4x faster, BUT with interior-point, Gaugin 100x100
       steps 3-6 converged to worse local minima (distortion 1.3-2.5x):
       interior-point's barrier sums over all constraints, so removing
       redundant ones changes its path. Also makes 'sqp' usable (with all
       rows sqp needs a ~107 GB matrix).

10. [ ] When fmincon reports infeasibility, enable
        `'EnableFeasibilityMode', true` with
        `'SubproblemAlgorithm', 'cg'` — the documented interior-point
        remedy. Relevant because for band-constrained steps the identity
        start is usually infeasible by construction.
        STATUS: not done. "No feasible point" (exitflag -2) still occurs
        occasionally (e.g. 3 of 51 candidates in Gaugin 100x100 update 2).

11. [ ] Harden the sweep endpoints. Targets are `linspace` between the
        lambda=0 and lambda=1 solutions (`computeSweep.m` ~line 124); the
        lambda=1 endpoint is a single fmincon run, and if it
        under-converges the whole target ladder is misplaced (and then
        cached in `endpoints.mat`). Give that run extra iterations /
        multistart and check its exitflag.
        STATUS (10/29): confirmed as a real risk. With 'sqp' (update 3),
        Gaugin 128x128's lambda=1 run reached only info -0.680 vs -0.543
        with interior-point, so the whole target ladder covered ~47% of the
        real range (sweep looks almost unchanged end to end). Other three
        images reached the same max with both methods. Earlier, a 2D
        plane search confirmed the interior-point lambda=1 endpoints are the
        true max (fruit, flower1, Gaugin at 64x64). Fix not done: e.g. run
        lambda=1 with both algorithms / several starts and keep the best.

12. [ ] Housekeeping (not search-related but confusing):
        `code/Info_functions/computeInfo_regress.m` contains a function
        named `computeInfo_regressSquared` (MATLAB dispatches on filename
        so it silently works); and `lossFunction.m` special-cases
        `isequal(infoFcn, @computeInfo_regressCIE)` by handle identity,
        which is fragile. Rename/clean up.
        STATUS: mostly done. Unused info functions deleted (10/22);
        lossFunction cleaned up (10/28: handle check -> nargout, unused
        'distortion' mode removed, callers documented). Remaining unused
        files listed in Findings below.

## Current state and findings (as of 2026-10-29)

### Code settings right now
- Info: `computeInfo_regressCIE`, plain regression (ridge was tried and
  removed), `abs()` info normalizer in computeSweep. Runner uses
  predictingWhat `'L,M,S'` from `'L and S'`.
- fmincon (colorCorrectionOptimize.m): interior-point, central finite
  differences, fminconFactor 1e8 (lossFunction.m), tolerances 1e-10,
  MaxIterations 200. The 'sqp' + MaxFunctionEvaluations 3000 lines
  (update 3) are present but commented out; comments say how to switch.
- Gamut (buildGamutConstraints.m): `useConvexHull = true`.
- i.e. the code is currently set to UPDATE 2 (see updateLog.txt for
  updates 0-3). Switching updates is done by editing these lines by hand.
- TEMPORARY: `updateNumber` in runAllTestImages + the 'updateNumber'
  option of computeSweep exist only for comparing updates. Once a
  preferred update is chosen, remove them and keep only that update's
  settings.

### What we learned
- Why images barely change as info increases: the original image already
  scores high (M is 96-99% predictable from L and S; the missing
  red-green part is only ~0.03-0.08 cone contrast). The regression info
  score rewards an ACCURATE copy of that signal: injecting it into L
  helps rebuild M but corrupts the rebuild of L, so info peaks near the
  signal's natural strength (k ~ 1-1.5). Larger, more visible changes are
  never rewarded. This is a property of the metric, not the optimizer.
- Ridge regression (noise sigma 0.03 on the dichromat's channels) did not
  fix this: the injected signal rides on the large L/S channels, so ridge
  barely penalizes it; instead ridge rewarded remixing L and S (large
  image changes, no extra red-green). Removed.
- 'M' vs 'L,M,S': 'M' widens the info range for fruit/Gaugin, same for
  flower1, fewer stalls. Not yet decided which to use. 'L-M' as written
  penalizes exactly the useful transforms; its case was removed.
- Open question for the project: should the transform deliver the
  red-green information accurately (current metric) or exaggerate it to
  make it visible (would need a visibility/discriminability metric)?
- Local minima are now the main optimizer problem: several steps have
  more than one solution and small changes (start point, constraint set,
  algorithm) decide which one fmincon finds. Examples: flower1 M 100x100
  step 9 (0.066 found, 0.041 exists, feasible); Gaugin steps 3-6 with
  hull + interior-point. Starting fmincon even AT the better point still
  ended at 0.066 (both interior-point and sqp), so more start points
  alone may not suffice.
- 'sqp' (update 3, 128x128 all 4 images): faster (~1 h vs ~1 h 45 min),
  same max info for 3 images, but for Gaugin its lambda=1 endpoint run
  under-converged (see #11), compressing that sweep's target range.
- Image size: 40x40 and 60x60 runs show stalls but no stalls ending at
  worse answers; only 100x100+ reproduced that. Compare runs only at the
  same size (downsampling changes the image slightly).

### Remaining unused code (candidates for deletion)
code/dichromatCorrection.m (calls deleted GlobalPCA files),
code/infoDistortionPlots.m, code/Distortion_functions/computeDistortion_luminance.m,
toolboxes/ColorSpaceConversion/RGB2LMSCalFormat.m,
toolboxes/CreateDirchomConfusion/plateO.m, toolboxes/Helpers/loadImage.m,
toolboxes/RenderImageAcrossPhenotype/DichromatSimulateBrettel*.m and
tri2dichromatLMSCalFormat.m, tutorials/t_daltonize.m (syntax error),
tutorials/t_junk.m, the old @daltonize methods compute.m /
computeInfoSweep.m / computeInfoDistortionGrid.m (old
colorCorrectionOptimize signature; also remove their wrappers in
daltonize.m), and *.asv autosaves. Keep computeDistortion_squared.m and
LMS2RGBCalFormat.m (used by generateGrayImage). t_applyTransforms.m and
t_compareFixedMetrics.m still list deleted info functions.
