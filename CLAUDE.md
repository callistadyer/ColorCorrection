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
      to tackle these with Claude's help).

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

1. [ ] Reformulate the band constraint. In
       `colorCorrectionOptimize.m` (`constraintBand`, ~line 318) the
       constraint is `(metric - target)^2 - eps^2 <= 0`. Its gradient
       `2*(metric-target)*grad(metric)` vanishes at the band center and is
       only ~2*eps near the boundary — nearly useless for interior-point.
       Replace with two one-sided constraints, linear in the metric:
       `c = [metric - target - eps; target - metric - eps]`. Same feasible
       set, well-scaled gradients.

2. [ ] Remove `fminconFactor = 1e8` in `lossFunction.m` (~line 77). The
       metrics are already normalized to O(1); scaling by 1e8 makes the
       default absolute `OptimalityTolerance` (1e-6) unreachable, so runs
       stop on StepTolerance/MaxIterations instead of optimality, and it
       skews the interior-point merit function. Set it to 1 and tune
       tolerances deliberately if needed.

3. [ ] Stop ignoring `exitflag`. In `colorCorrectionOptimize.m` (~line 190)
       `exitflag`/`output` are discarded. Log them into the saved per-
       candidate .mat files (chooseStartPoints caching), and when
       exitflag == 0 (hit MaxIterations = 200) warm-restart fmincon from
       the returned point (e.g. `while exitflag == 0 && nRestarts < 3`).
       exitflag == -2 (no feasible point) should also be surfaced.

4. [ ] Loosen the post-hoc feasibility slack. `colorCorrectionOptimize.m`
       (~line 221) accepts solutions only if maxC <= 1e-12, stricter than
       fmincon's own ConstraintTolerance (1e-10), so legitimately feasible
       solutions get rejected and replaced by T_init, discarding the whole
       run. Use the solver's ConstraintTolerance, and on violation return
       the optimized T with its violation recorded (chooseStartPoints
       already has a smallest-violation fallback that can use it) instead
       of reverting to T_init.

5. [ ] Fix the band-width inconsistency for negative info targets. The
       optimizer uses `pctTol * abs(targetInfo)`
       (`colorCorrectionOptimize.m` ~line 147) but the violation check in
       `chooseStartPoints.m` (~line 172) uses `pctTol * thisX` with no
       abs(). With `computeInfo_regressCIE`, normalized info is negative,
       so the checker's band collapses to absTolFloor while the optimizer
       worked to a much wider band — feasible solutions get classified as
       violations. Factor the band computation into one shared helper used
       by both. (Also note the two call sites default different
       pctTol/absTolFloor values; unify.)

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

7. [ ] (More work, biggest robustness win) Supply analytic gradients via
       `'SpecifyObjectiveGradient', true`. Squared distortion is exactly
       quadratic in the 9 transform parameters; the regression info has a
       closed-form gradient by the envelope theorem (at the least-squares
       beta, only dX/dT contracted with the residual is needed).
       Eliminates finite-difference noise (a common cause of StepTolerance
       stalls) and cuts function evaluations ~10x.

8. [ ] Add random restarts to the existing multistart. The candidate-start
       machinery in `chooseStartPoints.m` is deterministic (identity,
       previous step, info-sweep seeds, presolve seeds). Add 3-5 jittered
       starts per step (`T_candidate + 0.05*randn(3)`) to the
       candNames/candTinit lists. The `rng(1)` in
       `colorCorrectionOptimize.m` (~line 103) currently does nothing
       since nothing downstream is stochastic. MATLAB's MultiStart (Global
       Optimization Toolbox) is a drop-in alternative.

9. [ ] Prune gamut constraints to the convex hull.
       `buildGamutConstraints.m` builds 12N x 9 rows (10,800 for a 30x30
       image), but only extreme colors can bind. `convhulln` of the 3D
       RGB-contrast point cloud (and separately of its dichromat-mapped
       image) gives an exactly equivalent constraint set, typically tens
       of rows — faster iterations, better conditioning, and budget for
       more restarts.

10. [ ] When fmincon reports infeasibility, enable
        `'EnableFeasibilityMode', true` with
        `'SubproblemAlgorithm', 'cg'` — the documented interior-point
        remedy. Relevant because for band-constrained steps the identity
        start is usually infeasible by construction.

11. [ ] Harden the sweep endpoints. Targets are `linspace` between the
        lambda=0 and lambda=1 solutions (`computeSweep.m` ~line 124); the
        lambda=1 endpoint is a single fmincon run, and if it
        under-converges the whole target ladder is misplaced (and then
        cached in `endpoints.mat`). Give that run extra iterations /
        multistart and check its exitflag.

12. [ ] Housekeeping (not search-related but confusing):
        `code/Info_functions/computeInfo_regress.m` contains a function
        named `computeInfo_regressSquared` (MATLAB dispatches on filename
        so it silently works); and `lossFunction.m` special-cases
        `isequal(infoFcn, @computeInfo_regressCIE)` by handle identity,
        which is fragile. Rename/clean up.
