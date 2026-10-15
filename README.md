# ColorCorrection

Evaluates image transformations designed to improve the color information
available to dichromatic observers while minimizing distortion of the
original image. Transformations are found by constrained optimization,
either minimizing distortion subject to a target information score or
maximizing information subject to a target distortion level.

## Dependencies

This project depends on the following MATLAB toolboxes. The list below is
the full dependency tree, resolved transitively, so if you obtain all of
these you have everything you need.

| Toolbox | GitHub URL | Required by |
| --- | --- | --- |
| Psychtoolbox-3 | https://github.com/Psychtoolbox-3/Psychtoolbox-3 | ColorCorrection |
| ISETBioCSFGenerator | https://github.com/isetbio/ISETBioCSFGenerator | ColorCorrection |
| isetbio | https://github.com/isetbio/isetbio | ISETBioCSFGenerator |
| isetcam | https://github.com/iset/isetcam | isetbio |
| mQUESTPlus | https://github.com/brainardlab/mQUESTPlus | ISETBioCSFGenerator |

All of these need to be on your MATLAB path. For Psychtoolbox-3, it is the
`Psychtoolbox` subfolder of the repository that goes on the path.

## Recommended setup: ToolboxToolbox

The easy way to satisfy the dependencies is with
[ToolboxToolbox](https://github.com/ToolboxHub/ToolboxToolbox), a MATLAB
toolbox manager that reads this repository's configuration
([configuration/ColorCorrection.json](configuration/ColorCorrection.json)),
fetches the dependencies listed above, and sets up the MATLAB path for you.

With ToolboxToolbox installed and this repository cloned into your
ToolboxToolbox projects folder, simply run:

```matlab
tbUseProject('ColorCorrection');
```

This also deploys and runs the project's local hook (created from
[configuration/ColorCorrectionLocalHookTemplate.m](configuration/ColorCorrectionLocalHookTemplate.m)),
which sets the MATLAB preferences described in the next section.

## Manual setup

If you prefer not to use ToolboxToolbox, clone this repository and each of
the dependencies in the table above, add them to your MATLAB path, and set
the MATLAB preferences described in the next section by hand.

## Local configuration: MATLAB preferences

The code does not hard-code the locations of input and output data.
Instead, it looks these up at run time from MATLAB preferences in the
`ColorCorrection` preference group, and you need to point these at
directories on your local machine before running the code.

The preference currently used is:

- `outputDir` — directory where the analysis scripts and tutorials write
  their output. Point this anywhere convenient on your machine. Note that
  this directory also contains data read by the code: the test images the
  analyses operate on live in a `ColorCorrectionImages` subfolder of
  `outputDir`, and intermediate results written by earlier stages of the
  pipeline are read back by later ones.

The test images are not included in this repository — we keep code and
data separate. We plan to post them on an open science repository;
until then, they are available from the authors upon request.

How you set the preferences depends on your setup:

- **ToolboxToolbox users.** The first run of
  `tbUseProject('ColorCorrection')` copies the local hook template into
  your ToolboxToolbox local hooks folder as `ColorCorrectionLocalHook.m`.
  Edit that local copy (not the template in this repository) so that the
  preferences point to your directories; the hook then sets them each time
  the project is deployed.

- **Everyone else.** Set the preferences directly, e.g.:

  ```matlab
  setpref('ColorCorrection', 'outputDir', '/path/where/you/want/output');
  ```

  MATLAB preferences persist across sessions, so this only needs to be
  done once per machine.

## Repository organization

- [code/](code/) — core routines and top-level analysis scripts (e.g.
  [code/runAllTestImages.m](code/runAllTestImages.m))
- [tutorials/](tutorials/) — tutorial scripts
- [validation/](validation/) — validation scripts
- [configuration/](configuration/) — ToolboxToolbox configuration and
  local hook template
