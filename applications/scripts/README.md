---
sort: 5
---

# Bash functions: `solids4FoamScripts.sh`

---

Prepared by Ivan Batistić

---

## Section Aims

- This document describes the bash functions in the
  `solids4foam/applications/scripts/solids4FoamScripts.sh` script;
- These bash functions within `solids4FoamScripts.sh` are used in `Allrun` and
  `Allclean` scripts in `solids4foam` tutorial cases;
- The primary purpose of these functions is to make a case compatible with the
  current version of `OpenFOAM`/`foam-extend`, e.g. convert the case from its
  stored `OpenFOAM.com` format to the `foam-extend-4.1` format.

---

## Stored case format

Tutorial cases are stored in the `OpenFOAM.com` format, as this is the variant
used by most users. A case is converted to the format required by the loaded
`OpenFOAM`/`foam-extend` version when it is run, and converted back when it is
cleaned, so that the stored files are left unchanged.

Files that cannot be converted with a text substitution are stored twice. The
unsuffixed name always holds the `OpenFOAM.com` version, and the `foam-extend`
version is stored alongside it with a `.foamextend` suffix, e.g.
`system/functions` and `system/functions.foamextend`. During conversion the
unsuffixed file is backed up with a `.openfoam` suffix and the `foam-extend`
file is copied into its place; restoring moves the backup back.

---

## `solids4foam::convertCaseFormat()` and `solids4foam::restoreCaseFormat()`

`solids4foam::convertCaseFormat()` converts a case from the stored
`OpenFOAM.com` format to the format required by the loaded version. No changes
are made when `OpenFOAM.com` is loaded.
`solids4foam::restoreCaseFormat()` converts a case back to the stored
`OpenFOAM.com` format, regardless of which version is loaded. It is idempotent,
so it is safe to run it more than once.

- **Function purpose**
  Converts a case between the stored `OpenFOAM.com` format and the format of
  the loaded `OpenFOAM`/`foam-extend` version.

- **Function arguments**
  Path to case directory (most often `.` is used, referring to the current
  directory)

- **Example of usage**

  ```bash
  #!/bin/bash

  # Source solids4Foam scripts
  source solids4FoamScripts.sh

  # Convert the case to the format of the loaded version
  solids4Foam::convertCaseFormat .

  # Convert the case back to the stored OpenFOAM.com format
  solids4Foam::restoreCaseFormat .
  ```

```note
`solids4foam::convertCaseFormatFoamExtend()` is retained as a deprecated alias
for `solids4foam::restoreCaseFormat()`, so that existing user case scripts keep
working. It should not be used in new cases.
```

Changes are performed in different ways:

- Relocating files within case structure (see points 1 and 2);
- Using [sed](https://www.gnu.org/software/sed/manual/sed.html) command to
  perform insertion, deletion, and substitution (see points 3 and 5);
- Having different versions of the same file. The right one is chosen depending
  on the `OpenFOAM` version used. The switch between files is performed simply
  by renaming it (see points 4 or 7).

```note
The following description lists the differences between the formats. The
`foam-extend` conversion applies them, and restoring a case applies them in
reverse. Points 5 and 10 also apply to `OpenFOAM.org`, which differs from
`OpenFOAM.com` in the sampled set type names and in the forces output file name.
```

**1.** `symmetry` in `OpenFOAM` becomes `symmetryPlane` in `foam-extend`.

`blockMeshDict` is located, and every occurrence of the `symmetry` keyword is
updated:

```c++
patches
(
    symmetry left
    (
        (8 9 20 19)
    )
...
```

is transformed into:

```c++
patches
(
    symmetryPlane left
    (
        (8 9 20 19)
    )
...
```

The same is done for the `constant/polyMesh/boundary` file:

```c++
left
{
    type            symmetryPlane;
    inGroups        1(symmetryPlane);
    nFaces          30;
    startFace       1930;
}
```

is transformed into:

```c++
left
{
    type            symmetry;
    inGroups        1(symmetry);
    nFaces          30;
    startFace       1930;
}
```

**2.** If it is found, `blockMeshDict` is moved to the `constant/polyMesh/`
directory, as this is where `foam-extend` reads it from

For solid and fluid simulations, `blockMeshDict` is stored in `system/`:

```plaintext
├── 0
├── constant
└── system
 └── blockMeshDict
```

and it is moved to the `constant/polyMesh` directory:

```plaintext
├── 0
├── constant
│   └── polyMesh
│       └── blockMeshDict
└── system
```

For fluid-solid interaction simulations, there may be two `blockMeshDict` files,
each located in the corresponding `solid/` and `fluid/` subdirectories:

```plaintext
├── 0
├── constant
└── system
    └── fluid
    │   └── blockMeshDict
    └── solid
        └── blockMeshDict
```

These are also relocated to the `constant/` subdirectories:

```plaintext
├── 0
├── constant
│   └── fluid
│   │   └── polyMesh
│   │       └── blockMeshDict
│   └── solid
│       └── polyMesh
│           └── blockMeshDict
└── system
```

**2.1.** Rename the functions file

The function file is used to specify the list of function objects and is loaded
at the bottom of the `controlDict`:

```c++
#include "./system/functions"
```

The stored `functions` file is backed up as `functions.openfoam`, and the
`foam-extend` version, `functions.foamextend`, is copied into its place:

```plaintext
└── system
    ├── controlDict
    ├── fvSchemes
    ├── fvSolution
    ├── functions
    └── functions.foamextend
```

is transformed into:

```plaintext
└── system
    ├── controlDict
    ├── fvSchemes
    ├── fvSolution
    ├── functions
    ├── functions.foamextend
    └── functions.openfoam
```

**3.** Find the `turbulenceProperties` file and rename the value for the
`simulationType` keyword:

```c++
simulationType  RAS;
```

is renamed to:

```c++
simulationType  RASModel;
```

Note that in a fluid-solid interaction case, `turbulenceProperties` file is
located in the `constant/fluid/` directory, while for fluid simulations, it is
in the `constant/`.

**4.** If found,`boundaryData` directory is renamed:

```plaintext
├── 0
├── constant
│   ├── boundaryData
│   ├── boundaryData.openfoam
│   └── polyMesh/
└── system
```

is transformed into:

```plaintext
├── 0
├── constant
│   ├── boundaryData.foamextend
│   ├── boundaryData
│   └── polyMesh/
└── system
```

**5.** If `sample` file is found in the `system/` directory and if the
[OpenFOAM.org](OpenFOAM.org) version is used, `uniform` is replaced with
`lineUniform`:

```c++
sets
(
    lineXX
    {
        type       uniform;
        axis                x;
        nPoints            50;
        start (0.05 1e-6 0.0005);
        end    (0.1 1e-6 0.0005);
    }
...
```

is transformed into:

```c++
sets
(
    lineXX
    {
        type   lineUniform;
        axis               x;
        nPoints             50;
        start (0.05 1e-6 0.0005);
        end    (0.1 1e-6 0.0005);
    }
...
```

**6.** If `p` file is found, and if the `type` keyword is set to
`timeVaryingUniformFixedValue,` its value is changed to `uniformFixedValue` by
commenting out the appropriate lines:

```c++
inlet
{
    //type          uniformFixedValue;
    type        timeVaryingUniformFixedValue;

    uniformValue  tableFile;
    "file|fileName"     "$FOAM_CASE/0/fluid/time-series";
    outOfBounds         clamp;
}

```

is transformed into:

```c++
inlet
{
    type          uniformFixedValue;
    //type        timeVaryingUniformFixedValue;

    uniformValue  tableFile;
    "file|fileName"     "$FOAM_CASE/0/fluid/time-series";
    outOfBounds         clamp;
}

```

**7.** If found, the version of the`changeDictionaryDict` file is updated:

```plaintext
├── 0
├── constant
└── system
    ├── changeDictionaryDict
    ├── changeDictionaryDict.foamextend
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

is transformed into:

```plaintext
├── 0
├── constant
└── system
    ├── changeDictionaryDict
    ├── changeDictionaryDict.openfoam
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

**8.** If found, the version of the`createPatchDict` file is updated:

```plaintext
├── 0
├── constant
└── system
    ├── createPatchDict
    ├── createPatchDict.foamextend
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

is transformed into:

```plaintext
├── 0
├── constant
└── system
    ├── createPatchDict
    ├── createPatchDict.openfoam
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

**9.** Solid cases are stored with the `pointCellsLeastSquares` gradient
method, as `OpenFOAM` requires it to account for boundary non-orthogonal
corrections. For `foam-extend`, it is replaced with `leastSquares`, which is
the equivalent scheme there:

```c++
gradSchemes
{
    default            pointCellsLeastSquares;
}
```

is transformed into:

```c++
gradSchemes
{
    default            leastSquares;
}
```

Note: In fluid-solid interaction cases, this change is performed only on
`fvSchemes,` which refers to solid and is located in `system/solid/fvSchemes`.

**10.** In case the `force.gnuplot` script is found, the path to the forces
output file is changed. This file is generated by the `forces` function object.
`OpenFOAM.com` writes `force.dat` in `postProcessing/fluid/forces/0/`, whereas
`foam-extend` writes `forces.dat` in `forces/0/`:

```bash
plot [0.1:] "< sed s/[\\(\\)]//g ./postProcessing/fluid/forces/0/force.dat" \
    u 1:2 w l
```

is transformed into:

```bash
plot [0.1:] "< sed s/[\\(\\)]//g forces/0/forces.dat" u 1:2 w l
```

Note that `OpenFOAM.org` also writes `forces.dat`, so this change is applied
for `OpenFOAM.org` as well as for `foam-extend`.

**11.** In case the `plot.gnuplot` script is found, the path to the
`sigma_surface.raw` is changed. `sigma_surface.raw` is an output file generated
after using the `sample` function object for post-processing results. In
`OpenFOAM` it is located in `postProcessing/sample.surfaces/1/`, whereas in
`foam-extend` it is located in `postProcessing/surfaces/1/`:

```bash
path = "postProcessing/sample.surfaces/1/sigma_surface.raw"
```

is transformed into:

```bash
path = "postProcessing/surfaces/1/sigma_surface.raw"
```

Note that the `OpenFOAM` path is derived from the name of the function object
entry in the `system/` directory, here `sample`.

---

## `solids4foam::foamFlavour()`

- **Function purpose**
  Echoes the loaded variant: `foamextend`, `com` or `org`.

- **Function arguments**
  None

- **Example of usage**

  ```bash
  if [[ $(solids4Foam::foamFlavour) == "foamextend" ]]
  then
      echo "foam-extend is loaded"
  fi
  ```

---

## `solids4foam::blockMeshDictDir()`

- **Function purpose**
  Echoes the directory `blockMesh` reads `blockMeshDict` from for the loaded
  variant, i.e. `system` for `OpenFOAM` and `constant/polyMesh` for
  `foam-extend`. It is intended for cases that generate their `blockMeshDict`,
  for example with `m4`.

- **Function arguments**
  Optional region name, for multi-region cases

- **Example of usage**

  ```bash
  BLOCK_MESH_DICT_DIR=$(solids4Foam::blockMeshDictDir)
  mkdir -p "${BLOCK_MESH_DICT_DIR}"
  m4 -P system/blockMeshDict.m4 > "${BLOCK_MESH_DICT_DIR}"/blockMeshDict
  ```

---

## `solids4foam::caseOnlyRunsWithFoamExtend()`

- **Function purpose**
  This function gives an error if the `foam-extend` version is not
  sourced/loaded.

- **Function arguments**
  None

- **Example of usage**

  ```bash
  #!/bin/bash

  # Source solids4Foam scripts
  source solids4FoamScripts.sh

  solids4Foam::caseOnlyRunsWithFoamExtend
  ```

---

## `solids4foam::caseDoesNotRunWithFoamExtend()`

- **Function purpose**
  This function gives an error if the [OpenFOAM.com](OpenFOAM.com) or
  [OpenFOAM.org](OpenFOAM.org) version is not sourced/loaded.

- **Function arguments**
  None

- **Example of usage**

  ```bash
  #!/bin/bash

  # Source solids4Foam scripts
  source solids4FoamScripts.sh

  solids4Foam::caseDoesNotRunWithFoamExtend
  ```

---

## `solids4foam::removeEmptyDirs()`

Function loops over time directories (whose name is composed of digits or digits
and the dot) and removes them if there are no results. Checked resulting fields
are `U`, `T`, `D`, `pointD`, `DD`, `pointDD`. The compressed version of these
files with `.gz` extension is also checked. The function also works when the
case is decomposed.

- **Function purpose**
  Remove empty time directories that are inadvertently created when running FSI
  cases with preCICE.

- **Function arguments**
  None

- **Example of usage**

```bash
#!/bin/bash

# Source solids4Foam scripts
source solids4FoamScripts.sh

# Remove empty time directories created by preCICE
solids4Foam::removeEmptyDirs
```

---

## `solids4foam::err()`

It will construct a message string with the current date time and timezone
offset. The message passed as an argument to the `err()` function will be
appended to this message string. The string is written to a file named
`error.txt` and is also displayed on the console as an error output. In case an
optional argument (file name) is prescribed, the context of the file name is
written to `errorCommandLog.txt` file.

- **Function purpose** This function is designed to handle and report errors in
  a script.

- **Function arguments** "error message" - stored to `error.txt`
  optional parameter - name of the log file which will be stored to
  `errorCommandLog.txt`

- **Example of usage**

  ```bash
  #!/bin/bash

  # Source solids4Foam scripts
  source solids4FoamScripts.sh

  # Check if a 0 directory already exists
  if [[ -d "postProcessing" ]]
  then
   solids4foam::err "The postProcessing directory already exists: run Allclean"
  fi
  ```

  In the above example, the error message is stored in the `error.txt` file and
  printed in the console as:

  ```plaintext
  ERROR: see error.txt
  [2023-09-05T10:19:18+0200]: The postProcessing directory already exists: run Allclean
  ```

  If an optional argument wants to be used, the command should look like this:

  ```bash
  #!/bin/bash

  # Source solids4Foam scripts
  source solids4FoamScripts.sh

  # Check if a 0 directory already exists
  if [[ -f "postProcessing/data.dat" ]]
  then
   solids4foam::err "The data.dat file already exists: run Allclean to delete" postProcessing/data.dat
  fi
  ```

  The console output of this will be :

  ```plaintext
  ERROR: see error.txt
  [2023-09-05T10:19:18+0200]: The postProcessing directory already exists: run
  Allclean to delete postProcessing/data.dat
         see errorCommandLog.txt
  ```

  The context of the `data.dat` is stored in the errorCommandLog.txt file.

---

## `solids4foam::runApplication()`

- **Function purpose** This function is designed to run `OpenFOAM` and
  `solids4foam` applications with additional logging and error handling.

- **Function options**
  `-a` = append (append to an existing log file)
  `-o` = overwrite (overwrite existing log file)
  `-s` = suffix (adding suffix to the log file)
  `-decomposeParDict <locationOfAlternativeDecomposeParDict>` Alternative option
  for `decomposeParDict` dictionary location
- **Function arguments**
  `<appName>` Name of the executable (command to run)

  _Note_: Parameters which are added **after** the executable will be passed on
  to it!

- **Example of usage**

```bash
#!/bin/bash

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Source solids4Foam scripts
source solids4FoamScripts.sh

# Create meshes
solids4Foam::runApplication -s solid blockMesh -region solid
```

In the example above, the word `solid` is the suffix. `blockMesh` is the name of
the executable with `-region solid` as an parameter for the`blockMesh`.

```bash
#!/bin/bash

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Source solids4Foam scripts
source solids4FoamScripts.sh

# Create cellZones for materials
solids4Foam::runApplication setSet -batch batch.setSet
solids4Foam::runApplication setsToZones

# Run the solver
solids4Foam::runApplication solids4Foam
```

In this example, options are not used. `setSet` is executable with
`-batch batch.setSet` parameter whereas `setsToZones` and `solids4Foam` are
executables without arguments.

In case the log file already exists and the `-a` option is not set, it will exit
and print that the executable is already running at this location.

Both standard error `sterr` and standard output `stdout` are redirected to the
log file. If the application returns an error (non-zero exit code), it appears
as “ERROR” in the log file. The name of the log file is `log.<executable name>`.
In case the option `-s` is activated, the name of the log file is
`log.<executable name>.<suffix>`.

---

## `solids4foam::runParallel()`

This function has the same functionalities and options as
`solids4Foam::runApplication` function described previously.
The difference is that the executable is run in parallel using MPI.

`mpirun -n nProcs executable -parallel >> log.executable`

is simply replaced with:

`solids4Foam::runParallel executable`

`decomposeParDict` located in the `system` is automatically checked to get the
number of processors `nProcs`.

In case that `$FOAM_MPI` is set to `msmpi`, `mpirun` is replaced with `msmpi`.
