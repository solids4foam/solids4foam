# solids4foam bash functions:  `solids4FoamScripts.sh`

Prepared by Ivan Batistić

---

- These functions are used in `Allrun` and `Allclean` scripts to run and clean `solids4foam` tutorial cases;
- The primarily purpose of these functions is to efficiently change between `OpenFOAM` versions by making corresponding changes in `OpenFOAM` case files;
- In this document, corresponding functions in the `/applications/scripts/solids4FoamScripts.sh` bash script are documented by showing changes in corresponding files.

---

## Convert case format functions: `solids4foam::convertCaseFormat()` and `solids4foam::convertCaseFormatFoamExtend()`

`solids4foam::convertCaseFormat()` is used to convert a case from `foam-extend` format to `OpenFOAM` format, whereas `solids4foam::convertCaseFormatFoamExtend()` is used to convert any case version to the `foam-extend` format.  No changes are applied if the case format matches the format of the `OpenFOAM` version.

- __Function purpose__
  Converts a case from `foam-extend` format to `OpenFOAM` format and vice versa. 

- __Function arguments__
  Path to case directory (most often `.` is used, referring to the current directory)

- __Example of usage__
  ```bash
  #!/bin/bash
  
  # Source solids4Foam scripts
  source solids4FoamScripts.sh
  
  # Convert case format to OpenFOAM
  solids4Foam::convertCaseFormat .
  
  # Convert case format back to foam-extend
  solids4Foam::convertCaseFormatFoamExtend .
  ```


Changes are performed in different ways: 

- Relocating files within case structure (see points 1 and 2);
- Using [sed](https://www.gnu.org/software/sed/manual/sed.html) command to perform insertion, deletion, and substitution (see points 3 and 5);
- Having different versions of the same file. The right one is chosen depending on the `OpenFOAM` version used. The switch between files is performed simply by renaming it (see points 4 or 7).  

```note
The following description refers to the function `solids4foam::convertCaseFormat()`. The description of `solids4foam::convertCaseFormatFoamExtend()` is not provided because it performs the same changes only in reverse order.
```


__1.__ `symmetryPlane` in `foam-extend` becomes `symmetry` in `OpenFOAM`.  

`blockMeshDict` is located and every occurrence of `symmetryPlane` keyword is updated:

```c++
patches
(
    symmetryPlane left
    (
        (8 9 20 19)
    )
...
```

is transformed into:

```c++
patches
(
    symmetry left
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

__2.__ If it is found, `blockMeshDict` is moved to the `/system` directory

For solid and fluid simulations, `blockMeshDict` is found in `constant/polyMesh/`:

```bash
├── 0
├── constant
│   └── polyMesh
│       └── blockMeshDict
└── system
```

and it is moved to the `system` directory:

```
├── 0
├── constant
│   └── polyMesh
└── system
	└── blockMeshDict
```

For FSI simulations, `blockMeshDict` can be used to generate mesh for both solid and fluid. In that case, `blockMeshDict` is located in the corresponding `solid` and `fluid` subdirectories:

```
├── 0
├── constant
│   └── fluid
│   │   └── polyMesh
│   │       └── blockMeshDict
│   └── solid
│       └── polyMesh
│           └── blockMeshDict
└── system
```

and is also relocated to the system:

```
├── 0
├── constant
│   └── fluid
│   │   └── polyMesh
│   └── solid
│       └── polyMesh
└── system
    └── fluid
    │   └── blockMeshDict
    └── solid
        └── blockMeshDict
```

__2.1.__ Rename the functions file

The function file is used to specify the list of function objects and is loaded at the bottom of the  `controlDict`:

```
#include "./system/functions"
```

In `foam-extend` case structure `functions` file is renamed to `functions.foam-extend` whereas `functions.openfoam` is renamed to functions:

```
└── system
    ├── controlDict
    ├── fvSchemes
    ├── fvSolution
    ├── functions
    └── functions.openfoam
```

is transformed into:

```
└── system
    ├── controlDict
    ├── fvSchemes
    ├── fvSolution
    ├── functions.foam-extend
    └── functions
```

__3.__ Find `turbulenceProperties` file and rename value for the  `simulationType` keyword:

```c++
simulationType  RASModel;
```

is renamed to:

```c++
simulationType  RAS;
```

Note that in the case of FSI simulation `turbulenceProperties` file is located in the `constant/fluid/` directory while for fluid simulations it is located in the `constant/`.

__4.__ If found,`boundaryData` directory is renamed:

```bash
├── 0
├── constant
│   ├── boundaryData
│   ├── boundaryData.openfoam
│   └── polyMesh/
└── system
```

is transformed into:

```bash
├── 0
├── constant
│   ├── boundaryData.foamextend
│   ├── boundaryData
│   └── polyMesh/
└── system
```

__5.__ If `sample` file is found in `system/` directory and if the [OpenFOAM.org](OpenFOAM.org) version is used,  `uniform` is replaced with`lineUniform`:

``` 
sets
(
    lineXX
    {
        type 	    	uniform;
        axis             	  x;
        nPoints          	 50;
        start (0.05 1e-6 0.0005);
        end    (0.1 1e-6 0.0005);
    }
...
```

is transformed into:

```
sets
(
    lineXX
    {
        type 		lineUniform;
        axis        	      x;
        nPoints       	     50;
        start (0.05 1e-6 0.0005);
        end    (0.1 1e-6 0.0005);
    }
...
```

__6.__ If `p` file is found, and if `type` keyword is set to `timeVaryingUniformFixedValue` its value is changed to `uniformFixedValue` by commenting out the appropriate lines:

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

__7.__ If found, version of the`changeDictionaryDict` file is updated:

```
├── 0
├── constant
└── system
    ├── changeDictionaryDict
    ├── changeDictionaryDict.openfoam
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

is transformed into:

```bash
├── 0
├── constant
└── system
    ├── changeDictionaryDict.foamextend
    ├── changeDictionaryDict
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

__8.__  If found, version of the`createPatchDict` file is updated:

```bash
├── 0
├── constant
└── system
    ├── createPatchDict
    ├── createPatchDict.openfoam
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

is transformed into:

```bash
├── 0
├── constant
└── system
    ├── createPatchDict.foamextend
    ├── createPatchDict
    ├── controlDict
    ├── fvSchemes
    └── fvSolution
```

__9.__ In case the [OpenFOAM.com](OpenFOAM.com) version is used to solve solid mechanics or FSI problems,  `leastSquare` gradient method in `fvSchemes` file is replaced with `pointCellsLeastSquares` to account for non-orthogonal correction:

```c++
gradSchemes
{
    default            leastSquares;
}
```

is transformed into:

```c++
gradSchemes
{
    default            pointCellsLeastSquares;
}
```

Note: In case of FSI simulation, this change is performed only on `fvSchemes` which refers to solid and is located in `system/solid/fvSchemes`.

__10.__ In case the `force.gnuplot` script is found, path to the `force.dat` file is changed. `force.dat` is an output file generated by the `forces` function object. When `foam-extend` is used, it is located in the `forces/0/` directory; otherwise it is located in `postProcessing/fluid/forces/0/`:

```bash
plot [0.1:] "< sed s/[\\(\\)]//g forces/0/forces.dat" u 1:2 w l
```

is transformed into:

```bash
plot [0.1:] "< sed s/[\\(\\)]//g ./postProcessing/fluid/forces/0/force.dat" u 1:2 w l
```

__11.__ In case the `plot.gnuplot` script is found, the path to the `sigma_surface.raw` is changed.  `sigma_surface.raw` is an output file generated after using `sample` utility for post-processing results. When `foam-extend` is used, it is located in the `"postProcessing/surfaces/1/` directory; otherwise, it is located in `postProcessing/sampleDict/1/`:

```bash
path = "postProcessing/surfaces/1/sigma_surface.raw"
```

is transformed into:

```bash
path = "postProcessing/sampleDict.v2012/1/sigma_surface.raw"
```

Note that the updated path has `sampleDict.v2012` in it, and this is because it is called the same in the `system/` directory where `sampleDIct.v2012` is located.

---

## Case only runs with `foam-extend`: `solids4foam::caseOnlyRunsWithFoamExtend()`

- __Function purpose__
  This function gives an error if the `foam-extend` version is not sourced.

- __Function arguments__
  None

- __Example of usage__

  ```bash
  #!/bin/bash
  
  # Source solids4Foam scripts
  source solids4FoamScripts.sh
  
  solids4Foam::caseOnlyRunsWithFoamExtend
  ```

---

